import warnings
from dataclasses import dataclass, asdict
from functools import cache
from pathlib import Path
from typing import Tuple

import numpy as np
import yaml
from arpeggio.core import InteractionComplex, SelectionError
from biotite.structure import AtomArray, Atom, filter_nucleotides
from biotite.structure.io.pdbx import CIFFile, get_structure, CIFBlock
from more_itertools import one
from pydash import py_, spread, is_none, filter_, map_, is_empty, keys
from tqdm.contrib.concurrent import process_map

from pdb_common import get_all_ligand_ref_atom_mappings, get_all_ligand_hetcodes, \
    get_hbonds_dir_path, NucleobaseHbonds, HydrogenBond, LigandHydrogenBonds, ComplexLigandsHydrogenBonds, \
    get_mmcif_dir_path, get_pdb_ids_file_path, get_protonated_mmcif_dir_path, get_workers


@dataclass(frozen=True)
class LigandDetails:
    chem_comp_id: str
    chain_id: str
    res_id: int


def get_ligands_in_structure(cif_block: CIFBlock) -> set[LigandDetails]:
    n = cif_block.get('pdbx_nonpoly_scheme')

    if n is None:
        return set()

    return set(
        py_([n['mon_id'].as_array(),
             n['pdb_strand_id'].as_array(),
             n['pdb_seq_num'].as_array().astype(int)])
        .unzip()
        .map(spread(LigandDetails))
        .filter(lambda x: x.chem_comp_id in get_all_ligand_hetcodes())
        .value())


def get_heavy_atom_names(atoms: AtomArray):
    return {atom.atom_name for atom in atoms if atom.element != 'H'}


def sort_hydrogen_bond_partners(x: Tuple[Atom, Atom, Atom], chain_id: str, res_id: int):
    if str(x[0].chain_id) == chain_id and int(x[0].res_id) == res_id:
        ligand_atom, protein_atom, hydrogen_atom = x
    else:
        protein_atom, ligand_atom, hydrogen_atom = x
    return ligand_atom, protein_atom, hydrogen_atom


sidechain_connector_atoms = {
    'SER': ['CB'],
    'THR': ['CB'],
    'TYR': ['CZ'],

    'GLN': ['CG', 'CD'],
    'ASN': ['CB', 'CG'],

    'GLU': ['CG', 'CD'],
    'ASP': ['CB', 'CG'],

    'LYS': ['CE'],
    'ARG': ['CZ', 'NE', 'NH1', 'NH2']
}

allowed_protein_atoms = {'OG', 'OG1', 'OH', 'NZ', 'NH1', 'NH2', 'NE', 'OD1', 'OD2', 'OE1', 'OE2', 'ND2', 'NE2'}


@cache
def get_fragment_and_atom_mapping(res_name: str):
    ref_mapping = get_all_ligand_ref_atom_mappings()[res_name]
    fragment_and_atom_mapping = {}
    for fragment, atom_mapping in ref_mapping.items():
        for res_atom, ref_atom in atom_mapping.items():
            fragment_and_atom_mapping[res_atom] = (fragment, ref_atom)
    return fragment_and_atom_mapping


def process_interaction(x, ligand: LigandDetails, s: AtomArray):
    if not any(c in x['contact'] for c in ['hbond', 'weak_hbond']):
        return None

    if (x['bgn']['auth_asym_id'], x['bgn']['auth_seq_id']) == (ligand.chain_id, ligand.res_id):
        ligand_side = x['bgn']
        protein_side = x['end']
    elif (x['end']['auth_asym_id'], x['end']['auth_seq_id']) == (ligand.chain_id, ligand.res_id):
        ligand_side = x['end']
        protein_side = x['bgn']
    else:
        return None

    protein_chain_id = protein_side['auth_asym_id']
    protein_res_id = protein_side['auth_seq_id']
    protein_res_name = protein_side['label_comp_id']
    protein_atom_name = protein_side['auth_atom_id']

    if protein_res_name not in keys(sidechain_connector_atoms) or protein_atom_name in ['N', 'O',
                                                                                        'CA'] or protein_atom_name not in allowed_protein_atoms:
        return None

    fragment_and_atom_mapping = get_fragment_and_atom_mapping(ligand.chem_comp_id)

    ligand_atom_name = ligand_side['auth_atom_id']
    if ligand_atom_name not in fragment_and_atom_mapping:
        return None

    ligand_atoms = s[
        (s.chain_id == ligand.chain_id) & (s.res_id == ligand.res_id) & (s.atom_name == ligand_atom_name) & np.isin(
            s.label_alt_id, ['', '.', '?'])]
    protein_atoms = s[
        (s.chain_id == protein_chain_id) & (s.res_id == protein_res_id) & (s.atom_name == protein_atom_name) & np.isin(
            s.label_alt_id, ['', '.', '?'])]
    if len(protein_atoms) != 1 or len(ligand_atoms) != 1:
        return None

    ligand_atom = ligand_atoms[0]
    protein_atom = protein_atoms[0]

    fragment_name, ref_atom_name = fragment_and_atom_mapping[ligand_atom_name]
    ligand_atom_coord = ligand_atom.coord.tolist()

    protein_atom_coord = protein_atom.coord.tolist()
    if 'pdbx_sifts_xref_db_name' in protein_atom._annot and protein_atom.pdbx_sifts_xref_db_name == 'UNP':
        uniprot_accession = str(protein_atom.pdbx_sifts_xref_db_acc)
        uniprot_res_id = int(protein_atom.pdbx_sifts_xref_db_num)
    else:
        uniprot_accession = None
        uniprot_res_id = None

    sidechain_connector_coords = {
        atom_name: s[(s.chain_id == protein_chain_id) & (s.res_id == protein_res_id) & (s.atom_name == atom_name)][
            0].coord.tolist() for atom_name in sidechain_connector_atoms[protein_res_name] if
        atom_name != protein_atom_name}

    return fragment_name, HydrogenBond(
        ref_atom_name,
        ligand_atom_name,
        ligand_atom_coord,
        protein_atom_coord,
        protein_chain_id,
        protein_res_id,
        protein_res_name,
        protein_atom_name,
        sidechain_connector_coords,
        uniprot_accession,
        uniprot_res_id
    )


def calculate_ligand_hydrogen_bonds(ic: InteractionComplex, s: AtomArray, ligand: LigandDetails, is_nucleic_acid: bool):
    try:
        ic.run_arpeggio([f'/{ligand.chain_id}/{ligand.res_id}/'], 5, 0.1, False)
    except SelectionError:
        return None
    contacts = ic.get_contacts()

    raw_fragment_hbonds: dict[str, list[HydrogenBond]] = (py_(contacts)
                                                          .filter(lambda x: x['interacting_entities'] == 'INTER')
                                                          .map(lambda x: process_interaction(x, ligand, s))
                                                          .reject(is_none)
                                                          .group_by(lambda x: x[0])
                                                          .map_values(lambda x: map_(x, lambda x: x[1]))
                                                          .value())

    if is_empty(raw_fragment_hbonds):
        return None

    ligand_structure = s[
        (s.chain_id == ligand.chain_id) & (s.res_id == ligand.res_id) & (s.res_name == ligand.chem_comp_id) & (
                    s.element != 'H')]
    fragments_with_hbonds: dict[str, NucleobaseHbonds] = {}
    atom_mapping = get_all_ligand_ref_atom_mappings()[ligand.chem_comp_id]

    for fragment, hbonds in raw_fragment_hbonds.items():
        nucleobase_plane_coords = ligand_structure[
            np.isin(ligand_structure.atom_name, filter_(atom_mapping[fragment].keys(), lambda a: a[0] in ('N', 'O')))][
            :3].coord
        if len(nucleobase_plane_coords) < 3:
            continue
        plane = map_(nucleobase_plane_coords, lambda x: map_(x, float))
        fragments_with_hbonds[fragment] = NucleobaseHbonds(plane, hbonds)

    return LigandHydrogenBonds(str(ligand.chain_id), int(ligand.res_id), str(ligand_structure[0].res_name),
                               is_nucleic_acid, fragments_with_hbonds)


def calculate_complex_hydrogen_bonds(protonated_mmcif_path: Path):
    pdb_id = protonated_mmcif_path.stem
    target_path = get_hbonds_dir_path() / f'{pdb_id}.yaml'

    # The unprotonated CIF file is used because some data is lost when ChimeraX exports the protonated structure
    cif_file = CIFFile.read(get_mmcif_dir_path() / f'{pdb_id}.cif')

    cif_block = one(cif_file.values())
    ligands_in_complex = get_ligands_in_structure(cif_block)

    try:
        structure = get_structure(cif_file, model=1, include_bonds=True,
                                  extra_fields=['label_asym_id', 'label_seq_id', 'label_alt_id',
                                                'pdbx_sifts_xref_db_name', 'pdbx_sifts_xref_db_acc',
                                                'pdbx_sifts_xref_db_num'])
    except KeyError:
        structure = get_structure(cif_file, model=1, include_bonds=True,
                                  extra_fields=['label_asym_id', 'label_seq_id', 'label_alt_id'])

    try:
        ic = InteractionComplex(str(protonated_mmcif_path))
    except RuntimeError:
        return

    ic.initialize()

    try:
        ligands_hbonds = py_(ligands_in_complex).map(
            lambda ligand: calculate_ligand_hydrogen_bonds(ic, structure, ligand, False)).reject(is_none).value()
        nucleic_acid_nucleotides = set(py_(structure[filter_nucleotides(structure)]).filter(
            lambda a: str(a.res_name) in get_all_ligand_hetcodes()).map(
            lambda a: LigandDetails(str(a.res_name), str(a.chain_id), int(a.res_id))).value()) - ligands_in_complex
        na_hbonds = py_(nucleic_acid_nucleotides).map(
            lambda ligand: calculate_ligand_hydrogen_bonds(ic, structure, ligand, True)).reject(is_none).value()
    except:
        # A negligible number of structures fail to get processed properly
        return

    nucleobase_residues = ligands_hbonds + na_hbonds
    data = ComplexLigandsHydrogenBonds(pdb_id, nucleobase_residues)

    try:
        yaml.safe_dump(asdict(data), target_path.open('wt'))
    except:
        target_path.unlink(missing_ok=True)
        raise


if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    get_hbonds_dir_path().mkdir(exist_ok=True)
    pdb_ids = set(get_pdb_ids_file_path().read_text().split('\n')) - set(
        map_(get_hbonds_dir_path().glob('*.yaml'), lambda x: x.stem))
    mmcif_paths = filter_(get_protonated_mmcif_dir_path().glob('*.cif'), lambda path: path.stem in pdb_ids)
    process_map(calculate_complex_hydrogen_bonds, mmcif_paths, max_workers=get_workers(), chunksize=10)

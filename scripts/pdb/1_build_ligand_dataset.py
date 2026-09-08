import warnings
from collections import defaultdict
from itertools import product

import numpy as np
import yaml
from biotite.structure import AtomArray
from biotite.structure.info import residue, all_residues
from pydash import reject, is_none, unzip
from tqdm.contrib.concurrent import process_map

from pdb_common import REF_HETCODES_RIBOSE_BOUND_ATOMS, atom_array_to_mol, get_ligands_file_path, get_workers


def get_ref_substructure_of_ligand(ref_hetcode: str, ligand_hetcode: str):
    """Returns a tuple of (nucleobase reference hetcode, ligand hetcode, (ligand atom -> reference atom))"""

    ref_structure: AtomArray = residue(ref_hetcode)
    ref_structure_without_h = ref_structure[ref_structure.element != 'H']
    ref_mol = atom_array_to_mol(ref_structure_without_h)
    ref_atom_names = [atom.GetPDBResidueInfo().GetName() for atom in ref_mol.GetAtoms()]

    try:
        ligand_structure: AtomArray = residue(ligand_hetcode)
        # Discard ligands with deuterium (D) or unknown (X) atoms
        if any(a.element in ('D', 'X') for a in ligand_structure):
            return None
    except:
        return None

    ligand_structure_without_h = ligand_structure[ligand_structure.element != 'H']
    ligand_mol = atom_array_to_mol(ligand_structure_without_h)

    matching_indice = ligand_mol.GetSubstructMatch(ref_mol)
    if len(matching_indice) < len(ref_atom_names):
        return None

    matching_atom_names = [str(ligand_structure_without_h[i].atom_name) for i in matching_indice]
    ligand_ref_mapping = dict(zip(matching_atom_names, ref_atom_names))

    for ligand_atom, ref_atom in ligand_ref_mapping.items():
        ligand_atom_index = np.where(ligand_structure_without_h.atom_name == ligand_atom)[0][0]
        ref_atom_index = np.where(ref_structure_without_h.atom_name == ref_atom)[0][0]

        # For all reference nucleobase atoms except the ribose-bound nitrogen, require identical number of covalent
        # bonds in the equivalent ligand atoms
        if ref_atom != REF_HETCODES_RIBOSE_BOUND_ATOMS[ref_hetcode] and len(
                ligand_structure_without_h.bonds.get_bonds(ligand_atom_index)[0]) != len(
            ref_structure_without_h.bonds.get_bonds(ref_atom_index)[0]):
            return None

        if ref_atom[0] == 'C':
            ligand_atom_index = np.where(ligand_structure.atom_name == ligand_atom)[0][0]
            ref_atom_index = np.where(ref_structure.atom_name == ref_atom)[0][0]
            if len(ligand_structure.bonds.get_bonds(ligand_atom_index)[0]) != len(
                    ref_structure.bonds.get_bonds(ref_atom_index)[0]):
                return None

    return ref_hetcode, ligand_hetcode, ligand_ref_mapping


if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    ligands_dict = defaultdict(dict)
    all_pairs = list(product(REF_HETCODES_RIBOSE_BOUND_ATOMS.keys(), all_residues()))
    ref_ligand_mapping = reject(
        process_map(get_ref_substructure_of_ligand, *unzip(all_pairs), max_workers=get_workers(), chunksize=10),
        is_none)
    for ref_hetcode, ligand_hetcode, atom_mapping in ref_ligand_mapping:
        ligands_dict[ligand_hetcode][ref_hetcode] = atom_mapping
    yaml.safe_dump(dict(ligands_dict), get_ligands_file_path().open('wt'))

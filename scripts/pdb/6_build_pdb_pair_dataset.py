import csv
from dataclasses import dataclass, fields, asdict
from functools import lru_cache
from itertools import combinations
from math import degrees
from pathlib import Path

import yaml
from biotite.structure.info import amino_acid_names
from biotite.structure.io.pdbx import CIFFile
from dacite import from_dict
from more_itertools import one
from pydash import py_, reject, reg_exp_replace, is_none, \
    map_, values
from sympy import Plane, Point3D, Line3D
from tqdm.contrib.concurrent import process_map

from pdb_common import HydrogenBond, ComplexLigandsHydrogenBonds, LigandHydrogenBonds, get_mmcif_dir_path, \
    get_hbonds_dir_path, get_data_dir, get_workers


@dataclass
class Pair:
    nucleic_acid: int
    ligand_ref: str
    ligand_name: str
    uniprot_accession: str
    uniprot_res_id: int
    pdb_id: str
    resolution: float
    ligand_chain_id: str
    ligand_res_id: int
    protein_chain_id: str
    protein_auth_seq_id: int
    protein_res_name: str
    ligand_atom_exact_name_1: str
    ligand_atom_exact_name_2: str
    ligand_atom_1: str
    ligand_atom_2: str
    protein_atom_exact_name_1: str
    protein_atom_exact_name_2: str
    protein_atom_1: str
    protein_atom_2: str
    binding_mode: str
    buckle: float
    propeller: float


class Group:
    OH = 'OH'
    NH3 = 'NH3'
    GUANIDINO_SINGLE_N = 'GUANIDINO_SINGLE_N'
    GUANIDINO_TWO_N = 'GUANIDINO_TWO_N'
    COO_SINGLE_O = 'COO_SINGLE_O'
    COO_TWO_O = 'COO_TWO_O'
    CONH2_BIFURCATED = 'CONH2_BIFURCATED'
    CONH2_BIDENTATE = 'CONH2_BIDENTATE'


def get_protein_group(res_name: str, heavy_atoms: set[str]):
    if res_name in ('SER', 'THR', 'TYR') and heavy_atoms in ({'OG'}, {'OG1'}, {'OH'}):
        return Group.OH

    if res_name == 'LYS' and heavy_atoms == {'NZ'}:
        return Group.NH3

    if res_name == 'ARG':
        if heavy_atoms in ({'NH1'}, {'NH2'}):
            return Group.GUANIDINO_SINGLE_N
        if heavy_atoms in ({'NE', 'NH1'}, {'NE', 'NH2'}):
            return Group.GUANIDINO_TWO_N
        if heavy_atoms == {'NH1', 'NH2'}:
            return Group.GUANIDINO_TWO_N
        return None

    if res_name in ('ASP', 'GLU'):
        if heavy_atoms in ({'OD1', 'OD2'}, {'OE1', 'OE2'}):
            return Group.COO_TWO_O
        if heavy_atoms in ({'OD1'}, {'OD2'}, {'OE1'}, {'OE2'}):
            return Group.COO_SINGLE_O
        return None

    if res_name in ('ASN', 'GLN'):
        if heavy_atoms in ({'OD1', 'ND2'}, {'OE1', 'NE2'}):
            return Group.CONH2_BIDENTATE
        if heavy_atoms in ({'OD1'}, {'OD2'}):
            return Group.CONH2_BIFURCATED
        return None

    return None


def get_binding_mode_name(ligand_ref_name: str, nucleobase_atoms_list: list[str], protein_res_name: str,
                          protein_atoms: list[str]):
    nucleobase_one_letter = ligand_ref_name[0]
    protein_group = get_protein_group(protein_res_name, set(protein_atoms))
    if protein_group is None:
        return None
    nucleobase_atoms = set(nucleobase_atoms_list)
    match nucleobase_one_letter:
        case 'A':
            if nucleobase_atoms == {'N1', 'N6'}:  # WC
                match protein_group:
                    case Group.CONH2_BIDENTATE:
                        return 'A-W(N1,N6)=CONH₂'
                    case Group.OH:
                        return 'A-W(N1,N6)>OH'
            if nucleobase_atoms == {'N6', 'N7'}:  # H
                match protein_group:
                    case Group.CONH2_BIDENTATE:
                        return 'A-H(N6,N7)=CONH₂'
                    case Group.OH:
                        return 'A-H(N6,N7)>OH'
            if nucleobase_atoms == {'C2', 'N3'}:  # S
                if protein_group == Group.CONH2_BIDENTATE:
                    return 'A-S(C2,N3)=CONH₂'
            if nucleobase_atoms == {'N7', 'C8'}:
                if protein_group == Group.CONH2_BIDENTATE:
                    return 'A-H(N7,C8)=CONH₂'

        case 'C':
            if nucleobase_atoms == {'O2', 'N3'}:  # W
                match protein_group:
                    case Group.GUANIDINO_TWO_N:
                        return 'C-W(O2,N3)=NHC(NH₂)₂⁺'
            if nucleobase_atoms == {'N4', 'C5'}:
                if protein_group == Group.COO_TWO_O:
                    return 'C-H(N4,C5)=COO⁻'
                if protein_group == Group.COO_SINGLE_O:
                    return 'C-H(N4,C5)<COO⁻'
            if nucleobase_atoms == {'N3', 'N4'}:
                match protein_group:
                    case Group.CONH2_BIDENTATE:
                        return 'C-W(N3,N4)=CONH₂'
                    case Group.OH:
                        return 'C-W(N3,N4)>OH'
            if nucleobase_atoms == {'N4', 'C5'} and protein_group == Group.OH:
                return 'C-H(N4,C5)>OH'

        case 'G':
            if nucleobase_atoms == {'O6', 'N7'}:  # H
                match protein_group:
                    case Group.GUANIDINO_TWO_N:
                        return 'G-H(O6,N7)=NHC(NH₂)₂⁺'
                    case Group.GUANIDINO_SINGLE_N:
                        return 'G-H(O6,N7)<NHC(NH₂)₂⁺'
                    case Group.NH3:
                        return 'G-H(O6,N7)>NH₃⁺'
            if nucleobase_atoms == {'N7', 'C8'}:
                match protein_group:
                    case Group.CONH2_BIDENTATE:
                        return 'G-H(N7,C8)=CONH₂'
            if nucleobase_atoms == {'N2', 'N3'}:  # S
                match protein_group:
                    case Group.CONH2_BIDENTATE:
                        return 'G-S(N2,N3)=CONH₂'
                    case Group.OH:
                        return 'G-S(N2,N3)>OH'
            if nucleobase_atoms == {'N1', 'O6'}:  # W
                match protein_group:
                    case Group.CONH2_BIDENTATE:
                        return 'G-W(N1,O6)=CONH₂'
                    case Group.OH:
                        return 'G-W(N1,O6)>OH'
            if nucleobase_atoms == {'N1', 'N2'}:
                match protein_group:
                    case Group.COO_TWO_O:
                        return 'G-W(N1,N2)=COO⁻'
                    case Group.OH:
                        return 'G-W(N1,N2)>OH'

        case 'U':
            if nucleobase_atoms == {'O4', 'C5'}:  # H
                match protein_group:
                    case Group.CONH2_BIDENTATE:
                        return 'U-H(O4,C5)=CONH₂'
                    case Group.OH:
                        return 'U-H(O4,C5)>OH'
            if nucleobase_atoms == {'O2', 'N3'}:
                match protein_group:
                    case Group.CONH2_BIDENTATE:
                        return 'U-W(O2,N3)=CONH₂'
                    case Group.OH:
                        return 'U-W(O2,N3)>OH'
            if nucleobase_atoms == {'N3', 'O4'}:
                match protein_group:
                    case Group.CONH2_BIDENTATE:
                        return 'U-W(N3,O4)=CONH₂'
                    case Group.OH:
                        return 'U-W(N3,O4)>OH'
            if protein_group == Group.GUANIDINO_TWO_N:
                if nucleobase_atoms == {'O2'}:
                    return 'U-W(O2)<NHC(NH₂)₂⁺'
                if nucleobase_atoms == {'O4'}:
                    return 'U-W(O4)<NHC(NH₂)₂⁺'
                if nucleobase_atoms == {'O2', ')4'}:
                    return 'U(W)[3]-Arg'

        case 'T':
            if nucleobase_atoms == {'O2', 'N3'}:
                match protein_group:
                    case Group.CONH2_BIDENTATE:
                        return 'T-W(O2,N3)=CONH₂'
                    case Group.OH:
                        return 'T-W(O2,N3)>OH'
            if nucleobase_atoms == {'N3', 'O4'}:
                match protein_group:
                    case Group.CONH2_BIDENTATE:
                        return 'T-W(N3,O4)=CONH₂'
                    case Group.OH:
                        return 'T-W(N3,O4)>OH'
            if protein_group == Group.GUANIDINO_TWO_N:
                if nucleobase_atoms == {'O2'}:
                    return 'T-W(O2)<NHC(NH₂)₂⁺'
                if nucleobase_atoms == {'O4'}:
                    return 'T-W(O4)<NHC(NH₂)₂⁺'


@dataclass
class ChainResidueRange:
    chain_id: str
    start_residue: int
    end_residue: int


def get_pseudopair(pdb_id: str, resolution, l: LigandHydrogenBonds, ligand_ref_name: str, nucleobase_plane: Plane,
                   i1: HydrogenBond, i2: HydrogenBond):
    if not (i1.protein_chain_id == i2.protein_chain_id
            and i1.protein_auth_seq_id == i2.protein_auth_seq_id
            and i1.protein_res_name in amino_acid_names()
            and (i1.protein_atom, i1.ref_ligand_atom) != (i2.protein_atom, i2.ref_ligand_atom)
            and (i1.protein_atom != i2.protein_atom or i1.ligand_atom != i2.ligand_atom)
    ):
        return None

    protein_res_name = i1.protein_res_name if any(i.protein_atom not in ('N', 'O') for i in (i1, i2)) else ''
    if not protein_res_name:
        return None
    binding_mode = get_binding_mode_name(ligand_ref_name, [i1.ref_ligand_atom, i2.ref_ligand_atom], protein_res_name,
                                         [i1.protein_atom, i2.protein_atom])
    if not binding_mode:
        return None

    i1, i2 = sorted((i1, i2), key=lambda x: x.ref_ligand_atom)

    sidechain_hbonded_points: list[Point3D] = [Point3D(*p) for p in [i1.protein_atom_coord, i2.protein_atom_coord]]

    sidechain_connector_coords = None
    match len(i1.sidechain_connector_coords):
        case 1:
            sidechain_connector_coords = [py_(i1.sidechain_connector_coords).values().head().value(),
                                          i1.protein_atom_coord]
        case 2:
            sidechain_connector_coords = values(i1.sidechain_connector_coords)
        case _:  # Guanidinium (Arginine)
            assert i1.protein_res_name == 'ARG', f'{pdb_id} {i1.protein_res_name} {len(i1.sidechain_connector_coords)}'
            if i1.protein_atom == i2.protein_atom:  # bifurcated
                sidechain_connector_coords = [i1.sidechain_connector_coords['CZ'], i1.protein_atom_coord]
            else:  # bidentate
                sidechain_connector_coords = py_(i1.sidechain_connector_coords).keys().without(i1.protein_atom,
                                                                                               i2.protein_atom).map(
                    lambda x: i1.sidechain_connector_coords[x]).value()
                assert len(sidechain_connector_coords) == 2

    sidechain_connector_points = [Point3D(c) for c in sidechain_connector_coords]
    buckle_angle = abs(degrees(float(nucleobase_plane.angle_between(Line3D(*sidechain_connector_points)))))

    sidechain_midpoint = sidechain_hbonded_points[0].midpoint(sidechain_hbonded_points[1])
    sidechain_midpoint_projection = nucleobase_plane.projection(sidechain_midpoint)

    propeller_angle = None

    if i1.protein_atom != i2.protein_atom and i1.ligand_atom != i2.ligand_atom:
        sidechain_line: Line3D = Line3D(*sidechain_hbonded_points)
        propeller_angle = abs(degrees(float(nucleobase_plane.angle_between(sidechain_line))))

    return Pair(1 if l.is_nucleic_acid else 0,
                ligand_ref_name,
                l.ligand_name,
                i1.uniprot_accession,
                i1.uniprot_res_id,
                pdb_id,
                resolution,
                l.ligand_chain_id,
                l.ligand_res_id,
                i1.protein_chain_id,
                i1.protein_auth_seq_id,
                protein_res_name,
                i1.ligand_atom,
                i2.ligand_atom,
                i1.ref_ligand_atom,
                i2.ref_ligand_atom,
                i1.protein_atom,
                i2.protein_atom,
                reg_exp_replace(i1.protein_atom, r'\d', ''),
                reg_exp_replace(i2.protein_atom, r'\d', ''),
                binding_mode,
                buckle_angle,
                propeller_angle
                )


@lru_cache(None)
def get_resolution(pdb_id: str) -> float:
    cif_file = CIFFile.read(get_mmcif_dir_path() / f'{pdb_id}.cif')
    cif_block = one(cif_file.values())
    if 'refine' in cif_block:
        return float(max(cif_block['refine']['ls_d_res_high'].as_array()))
    if 'reflns' in cif_block:
        return float(max(cif_block['reflns']['d_resolution_high'].as_array()))
    raise ValueError(f'No resolution found for {pdb_id}')


def find_pseudopairs(path: Path):
    ci = from_dict(ComplexLigandsHydrogenBonds, yaml.safe_load(path.open('rt')))
    resolution = get_resolution(path.stem)

    pairs: list[Pair] = []

    for ligand in ci.nucleobase_residues:
        for ref_hetcode, fragment in ligand.fragments.items():
            hbonds = [i for i in fragment.hbonds if
                      i.protein_res_name in amino_acid_names() and i.protein_atom[0] in ('N', 'O', 'S')]
            plane = Plane(*map_(fragment.plane, Point3D))
            pairs.extend(reject([get_pseudopair(ci.pdb_id, resolution, ligand, ref_hetcode, plane, i1, i2) for i1, i2 in
                                 combinations(hbonds, 2)], is_none))
    return pairs


if __name__ == '__main__':

    fieldnames = map_(fields(Pair), 'name')
    pairs: list[Pair] = reject(
        process_map(find_pseudopairs, list(get_hbonds_dir_path().glob('*.yaml')), max_workers=get_workers(),
                    chunksize=10), is_none)

    with (get_data_dir() / 'pdb_pair_dataset.csv').open('wt') as file:
        dict_writer = csv.DictWriter(file, fieldnames)
        dict_writer.writeheader()
        for pseudopair in pairs:
            dict_writer.writerow(asdict(pseudopair))

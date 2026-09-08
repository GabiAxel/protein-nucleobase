import os
from argparse import ArgumentParser
from dataclasses import dataclass
from functools import cache
from pathlib import Path

import yaml
from biotite.structure import AtomArray
from pydash import capitalize, keys
from rdkit.Chem import Mol, RWMol, Conformer, Atom, AtomPDBResidueInfo, BondType
from rdkit.Geometry import Point3D

@dataclass
class HydrogenBond:
    ref_ligand_atom: str
    ligand_atom: str
    ligand_atom_coord: list[float]
    protein_atom_coord: list[float]
    protein_chain_id: str
    protein_auth_seq_id: int
    protein_res_name: str
    protein_atom: str
    sidechain_connector_coords: dict[str, list[float]]
    uniprot_accession: str | None
    uniprot_res_id: int | None

@dataclass
class NucleobaseHbonds:
    plane: list[list[float]]
    hbonds: list[HydrogenBond]

@dataclass
class LigandHydrogenBonds:
    ligand_chain_id: str
    ligand_res_id: int
    ligand_name: str
    is_nucleic_acid: bool
    fragments: dict[str, NucleobaseHbonds]


@dataclass
class ComplexLigandsHydrogenBonds:
    pdb_id: str
    nucleobase_residues: list[LigandHydrogenBonds]

REF_HETCODES_RIBOSE_BOUND_ATOMS = {
    'ADE': 'N9',
    'CYT': 'N1',
    'GUN': 'N9',
    'TDR': 'N1',
    'URA': 'N1'
}

ARG_WORKERS = 'workers'
ARG_DATA_DIR = 'data'
ARG_CHIMERAX_EXECUTABLE = 'chimerax'


argparser = ArgumentParser()
argparser.add_argument('-w', f'--{ARG_WORKERS}', type=int, default=os.cpu_count())
argparser.add_argument('-d', f'--{ARG_DATA_DIR}', type=Path, default=Path.home() / 'protein_nucleobase_data')
argparser.add_argument('-x', f'--{ARG_CHIMERAX_EXECUTABLE}', type=Path, required=False)


@cache
def get_args():
    return argparser.parse_args()


def get_path(arg_name: str) -> Path:
    return get_args()[arg_name]


def get_workers() -> int:
    return get_args()[ARG_WORKERS]


def get_data_dir():
    return get_args()[ARG_DATA_DIR]


def get_ligands_file_path():
    return get_data_dir() / 'ligands.yaml'


def get_pdb_ids_file_path():
    return get_data_dir() / 'pdb_ids.txt'


def get_mmcif_dir_path():
    return get_data_dir() / 'mmcif'


def get_protonated_mmcif_dir_path():
    return get_data_dir() / 'protonated_mmcif'


def get_hbonds_dir_path():
    return get_data_dir() / 'hbonds'


def atom_array_to_mol(s: AtomArray) -> Mol:
    """Converts a Biotite AtomArray object to an RDKit Mol object."""

    s = s[s.element != 'H']
    mol = RWMol()
    atom_indices = []
    conformer = Conformer(len(s))
    for i, a in enumerate(s):
        try:
            mol_atom = Atom(capitalize(a.element))
        except:
            mol_atom = Atom(0)
        info = AtomPDBResidueInfo()
        info.SetName(a.atom_name)
        mol_atom.SetMonomerInfo(info)
        atom_indices.append(mol.AddAtom(mol_atom))
        conformer.SetAtomPosition(i, Point3D(*a.coord.astype(float)))
    mol.AddConformer(conformer)
    for atom_index_1, atom_index_2, _ in s.bonds.as_array():
        mol.AddBond(int(atom_index_1), int(atom_index_2), BondType.SINGLE)
    return mol.GetMol()


@cache
def get_all_ligand_ref_atom_mappings() -> dict[str, dict[str, dict[str, str]]]:
    """Returns a nested dictionary of (ligand hetcode -> (fragment hetcode -> (ligand atom -> reference nucleobase atom)))"""

    return yaml.safe_load(get_path(get_ligands_file_path()).open('rt'))


@cache
def get_all_ligand_hetcodes() -> list[str]:
    """returns a list of all ligand hetcodes"""

    return keys(get_all_ligand_ref_atom_mappings())
from argparse import ArgumentParser
from functools import cache
from io import StringIO
from pathlib import Path
from tempfile import gettempdir

import numpy as np
from biotite.structure import molecule_iter
from biotite.structure.io.pdb import PDBFile
from pydash import py_
from rdkit.Chem import MolFromXYZFile, MolToPDBBlock
from rdkit.Chem.rdDetermineBonds import DetermineBonds

ARG_INPUT_PATH = 'input'
ARG_CHARGE = 'charge'
ARG_OUTPUT_PATH = 'output'
ARG_TEMP_DIR_PATH = 'tempdir'
ARG_WORKERS = 'workers'
ARG_MEMORY = 'memory'

argparser = ArgumentParser()
argparser.add_argument('-i', f'--{ARG_INPUT_PATH}', type=Path)
argparser.add_argument('-c', f'--{ARG_CHARGE}', type=int, required=False)
argparser.add_argument('-o', f'--{ARG_OUTPUT_PATH}', type=Path)
argparser.add_argument('-t', f'--{ARG_TEMP_DIR_PATH}', type=Path, default=gettempdir())
argparser.add_argument('-w', f'--{ARG_WORKERS}', type=int)
argparser.add_argument('-m', f'--{ARG_MEMORY}', type=str)


@cache
def _get_args():
    return argparser.parse_args()


def get_input_path() -> Path:
    return _get_args()[ARG_INPUT_PATH]


def get_output_path() -> Path:
    return _get_args()[ARG_OUTPUT_PATH]


def get_temp_dir_path() -> Path:
    return _get_args()[ARG_TEMP_DIR_PATH]


def get_workers() -> int:
    return _get_args()[ARG_WORKERS]


def get_memory() -> str:
    return _get_args()[ARG_MEMORY]


def get_charge() -> int:
    if ARG_CHARGE in _get_args():
        return _get_args()[ARG_CHARGE]
    stem = get_input_path().stem
    if '⁻' in stem:
        return -1
    if '⁺' in stem:
        return +1
    return 0


def get_nucleobase_and_sidechain(input_path: Path, charge: int) -> tuple:
    """Returns a tuple of nucleobase and sidechain structures. Assumes only the nucleobase has 4+ oxygen/nitrogen atoms."""

    mol = MolFromXYZFile(str(input_path))
    DetermineBonds(mol, charge=charge)
    io = StringIO(MolToPDBBlock(mol))
    pdb_file = PDBFile.read(io)
    s = pdb_file.get_structure(model=1, include_bonds=True)
    molecules = py_(molecule_iter(s)).map(
        lambda s: [len(s[np.isin(s.element, ['N', 'O'])]) >= 4, s]).from_pairs().value()
    nucleobase = molecules[True]
    sidechain = molecules[False]
    return nucleobase, sidechain

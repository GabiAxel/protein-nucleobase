from pathlib import Path

import psi4
from pydash import py_, join

from qc.qc_common import get_nucleobase_and_sidechain


def calculate_sapt(input_path: Path, sidechain_charge: int, output_path: Path, workers: int, memory: str,
                   sapt_method: str):
    psi4.set_num_threads(workers)
    psi4.set_memory(memory)

    nucleobase, sidechain = get_nucleobase_and_sidechain(input_path, sidechain_charge)
    nucleobase_xyz, sidechain_xyz = [py_(s).map(lambda a: join([a.element, *a.coord.tolist()], '\t')).join('\n').value()
                                     for s in [nucleobase, sidechain]]

    psi4_dimer = join(
        ['units Angstrom', f'{sidechain_charge} 1', '--', f'{sidechain_charge} 1', sidechain_xyz, '--', '0 1',
         nucleobase_xyz], '\n')

    psi4.set_output_file(str(output_path))

    dimer = psi4.geometry(psi4_dimer)
    psi4.set_options({
        'freeze_core': True,
        'occ_tolerance': 1.0e-6,
        'do_mbpt_disp': True
    })
    psi4.energy(sapt_method, molecule=dimer)

    print(f'{sapt_method} results written to {output_path}')

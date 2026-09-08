from pathlib import Path

from opi.input.blocks import BlockGeom
from opi.input.simple_keywords import Scf, Dft, Task

from orca_calculation import calculate_with_orca
from qc_common import get_memory, get_workers, get_input_path, get_charge, get_output_path, get_temp_dir_path


def optimize_geometry(input_path: Path, charge: int, output_path: Path, temp_dir_path: Path, workers: int, memory: int):
    """Performs geometry optimization using ORCA"""

    output = calculate_with_orca(input_path, charge, workdir, workers, memory, [Task.OPT, Dft.WB97X3C, Scf.TIGHTSCF],
                                 [BlockGeom(maxiter=1000)])

    output.parse()

    if not output.scf_converged() or not output.geometry_optimization_converged():
        print(f'Not converged. See output: {output.get_outfile()}')
        return

    optimized = output.get_structure()
    output_path.write_text(optimized.to_xyz_block())

    print(f'Output written to: {output_path.absolute()}')


if __name__ == '__main__':
    input_path = get_input_path()
    charge = get_charge()
    output_path = get_output_path()
    workdir = get_temp_dir_path()
    workers = get_workers()
    memory = int(get_memory())
    optimize_geometry(input_path, charge, output_path, workdir, workers, memory)

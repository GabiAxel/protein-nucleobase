import shutil

from opi.input.simple_keywords import Task, Dft, Scf

from orca_calculation import calculate_with_orca
from qc_common import get_input_path, get_charge, get_output_path, get_temp_dir_path, get_workers, get_memory

if __name__ == '__main__':
    input_path = get_input_path()
    charge = get_charge()
    output_path = get_output_path()
    workdir = get_temp_dir_path()
    workers = get_workers()
    memory = int(get_memory())
    output = calculate_with_orca(input_path, charge, workdir, workers, memory,
                                 [Task.ENGRAD, Dft.WB97X3C, Scf.TIGHTSCF])
    shutil.copy(output.get_outfile(), output_path)
    print(f'Analytical gradient analysis written to: {output_path.absolute()}')

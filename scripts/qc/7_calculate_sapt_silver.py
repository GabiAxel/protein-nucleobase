from psi4_sapt_calculation import calculate_sapt
from qc_common import get_input_path, get_charge, get_output_path, get_temp_dir_path, get_workers, get_memory

if __name__ == '__main__':
    input_path = get_input_path()
    charge = get_charge()
    output_path = get_output_path()
    workdir = get_temp_dir_path()
    workers = get_workers()
    memory: str = get_memory()
    calculate_sapt(input_path, charge, output_path, workers, memory, 'SAPT2+/aug-cc-pVDZ')
from pathlib import Path

from opi.core import Calculator
from opi.input.blocks import Block
from opi.input.structures import Structure


def calculate_with_orca(input_path: Path, charge: int, temp_dir_path: Path, workers: int,
                        memory: int, orca_keywords: list[str], orca_blocks: list[Block] = None):
    calc = Calculator(input_path.stem, temp_dir_path)
    calc.structure = Structure.from_xyz(input_path, charge=charge, multiplicity=1)
    calc.input.ncores = workers
    calc.input.memory = memory
    calc.input.add_simple_keywords(*orca_keywords)
    if orca_blocks is not None:
        calc.input.add_blocks(*orca_blocks)
    calc.write_input()
    calc.run()

    output = calc.get_output()
    if not output.terminated_normally():
        print(f'Not terminated normally. See output: {output.get_outfile()}')
        raise

    return output

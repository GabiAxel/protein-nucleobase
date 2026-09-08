import platform
import subprocess
from pathlib import Path

from pydash import py_
from tqdm import tqdm

from pdb_common import get_path, get_mmcif_dir_path, get_protonated_mmcif_dir_path, ARG_CHIMERAX_EXECUTABLE, get_args


def get_chimerax_executable_path():
    if ARG_CHIMERAX_EXECUTABLE in get_args():
        return get_path(ARG_CHIMERAX_EXECUTABLE)
    match platform.system():
        case 'Darwin':
            return list(Path('/Applications').glob('ChimeraX-*.app'))[-1] / 'Contents' / 'bin' / 'ChimeraX'
        case 'Linux':
            return Path('/usr/bin/chimerax')
        case 'Windows':
            return list(Path('C:\Program Files').glob('ChimeraX*'))[-1] / 'bin' / 'ChimeraX.exe'


if __name__ == '__main__':

    chimerax_executable_path = get_chimerax_executable_path()
    if not chimerax_executable_path.exists():
        print(f'Path {chimerax_executable_path} does not exist. Please set -x/--chimerax flag to ChimeraX executable.')

    get_protonated_mmcif_dir_path().mkdir(exist_ok=True, parents=True)
    paths: list[Path] = py_(get_mmcif_dir_path().glob('*.cif')).reject(
        lambda path: (get_protonated_mmcif_dir_path() / path.name).exists()).value()
    for path in tqdm(paths):
        subprocess.run([chimerax_executable_path, '--nogui', '--silent', '--exit', '--cmd',
                        f'open {path};addh;save {get_protonated_mmcif_dir_path() / path.name}'])

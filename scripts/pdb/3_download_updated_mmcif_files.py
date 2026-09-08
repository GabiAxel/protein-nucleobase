import requests
from tqdm.contrib.concurrent import process_map

from pdb_common import get_mmcif_dir_path, get_pdb_ids_file_path, get_workers


def download_updated_mmcif(pdb_id: str):
    updated_mmcif_path = get_mmcif_dir_path() / f'{pdb_id}.cif'
    if not updated_mmcif_path.exists():
        response = requests.get(f'https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id}_updated.cif')
        response.raise_for_status()
        updated_mmcif_path.write_bytes(response.content)


if __name__ == '__main__':
    get_mmcif_dir_path().mkdir(exist_ok=True, parents=True)
    pdb_ids = get_pdb_ids_file_path().read_text().split('\n')
    process_map(download_updated_mmcif, pdb_ids, max_workers=get_workers(), chunksize=10)

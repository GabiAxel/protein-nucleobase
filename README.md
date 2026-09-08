Python scripts for producing the data in

The script functionalities and arguments are detailed below. Arguments can be passed with either short or long flags.

Example of short flags:

```sh
python 1_build_ligand_dataset.py -d /my_data_dir -w 8
```

Example of short flags:

```sh
python 1_build_ligand_dataset.py --data /my_data_dir --workers 8
```

Example of no flags wth default values applied:

```sh
python 1_build_ligand_dataset.py
```

## Scripts

### Package [pdb](scripts/pdb) - PDB survey of nucleobase-side chain occurrences

⚠️ The scripts in this package rely on ordered execution with the same data directory (-d/--data flag).

#### [1_build_ligand_dataset.py](scripts/pdb/1_build_ligand_dataset.py)
Finds A, G, C, U and T nucleobase-containing molecules in the CCD and produces a YAML file with molecule-fragment-atom mapping nested dictionary.

**Arguments:**

| Flag               | Description       | Default                   |
|--------------------|-------------------|---------------------------|
| `-d`/ `--data`     | Data directory    | ~/protein_nucleobase_data |
| `-w` / `--workers` | Number of workers | Number of CPU cores       |


#### [2_find_pdb_ids_with_ligands_or_nucleic_acid.py](scripts/pdb/2_find_pdb_ids_with_ligands_or_nucleic_acid.py)
Searches RCSB PDB API for structures of proteins in complex with either any small molecule from the previous step or with a nucleic acid, determined with X-ray crystallography at resolution of 2.5Å or better, and saves the PDB IDs to a text file.

**Arguments:**

| Flag               | Description       | Default                   |
|--------------------|-------------------|---------------------------|
| `-d`/ `--data`     | Data directory    | ~/protein_nucleobase_data |
| `-w` / `--workers` | Number of workers | Number of CPU cores       |

#### [3_download_updated_mmcif_files.py](scripts/pdb/3_download_updated_mmcif_files.py)
Downloads the Updated mmCIF files from PDBe for the PDB IDs from the previous step.

**Arguments:**

| Flag               | Description       | Default                   |
|--------------------|-------------------|---------------------------|
| `-d`/ `--data`     | Data directory    | ~/protein_nucleobase_data |
| `-w` / `--workers` | Number of workers | Number of CPU cores       |

#### [4_protonate_mmcif_files.py](scripts/pdb/4_protonate_mmcif_files.py)
Protonates the structures from the previous step using ChimeraX.

**Arguments:**

| Flag                | Description                 | Default                                    |
|---------------------|-----------------------------|--------------------------------------------|
| `-d`/ `--data`      | Data directory              | ~/protein_nucleobase_data                  |
| `-w` / `--workers`  | Number of workers           | Number of CPU cores                        |
| `-c` / `--chimerax` | Path of ChimeraX executable | OS-dependent default installation location |

#### [5_calculate_hydrogen_bonds.py](scripts/pdb/5_calculate_hydrogen_bonds.py)
Calculates hydrogen bonds between nucleobases and side chain functional groups in the protonated structures from the previous step.

**Arguments:**

| Flag               | Description       | Default                   |
|--------------------|-------------------|---------------------------|
| `-d`/ `--data`     | Data directory    | ~/protein_nucleobase_data |
| `-w` / `--workers` | Number of workers | Number of CPU cores       |

#### [6_build_pdb_pair_dataset.py](scripts/pdb/6_build_pdb_pair_dataset.py)
Composes a dataset of bidentate and bifurcated hydrogen-bonded pair occurrences in the PDB based on the hydrogen bonds from the previous step.

**Arguments:**

| Flag               | Description       | Default                   |
|--------------------|-------------------|---------------------------|
| `-d`/ `--data`     | Data directory    | ~/protein_nucleobase_data |
| `-w` / `--workers` | Number of workers | Number of CPU cores       |

### Package [qc](scripts/qc) - Quantum-chemical calculations

Each script in this package contains a minimal "__main__" block that bridges the CLI usage to a main function. The functions can also be invoked in a any pipeline. 

#### [1_optimize_geometry.py](scripts/qc/1_optimize_geometry.py)
Performs geometry optimization using ORCA.

**Arguments:**

| Flag               | Description              | Default                |
|--------------------|--------------------------|------------------------|
| `-i` / `--input`   | Input XYZ file path      |                        |
| `-c` / `--charge`  | Charge                   | Inferred from filename |
| `-o` / `--output`  | Optimized XYZ file path  |                        |
| `-t` / `--tempdir` | Temporary work directory | OS-dependent           |
| `-w` / `--workers` | Number of workers        | Number of CPU cores    |
| `-m` / `--memory`  | Memory per worker        |                        |

#### [2_calculate_gradients.py](scripts/qc/2_calculate_gradients.py)
Performs analytical gradient analysis using ORCA.

**Arguments:**

| Flag               | Description              | Default                |
|--------------------|--------------------------|------------------------|
| `-i` / `--input`   | Input XYZ file path      |                        |
| `-c` / `--charge`  | Charge                   | Inferred from filename |
| `-o` / `--output`  | ORCA output file path    |                        |
| `-t` / `--tempdir` | Temporary work directory | OS-dependent           |
| `-w` / `--workers` | Number of workers        | Number of CPU cores    |
| `-m` / `--memory`  | Memory per worker        |                        |

#### [3_calculate_frequencies.py](scripts/qc/3_calculate_frequencies.py)
Performs harmonic vibrational frequency analysis using ORCA.

**Arguments:**

| Flag               | Description              | Default                |
|--------------------|--------------------------|------------------------|
| `-i` / `--input`   | Input XYZ file path      |                        |
| `-c` / `--charge`  | Charge                   | Inferred from filename |
| `-o` / `--output`  | ORCA output file path    |                        |
| `-t` / `--tempdir` | Temporary work directory | OS-dependent           |
| `-w` / `--workers` | Number of workers        | Number of CPU cores    |
| `-m` / `--memory`  | Memory per worker        |                        |

#### [4_get_plane_angles.py](scripts/qc/4_get_plane_angles.py)
Calculates the buckle and propeller angles for a nucleobase-sidechain dimer.

**Arguments:**

| Flag              | Description              | Default                |
|-------------------|--------------------------|------------------------|
| `-i` / `--input`  | Input XYZ file path      |                        |
| `-c` / `--charge` | Charge                   | Inferred from filename |

#### [5_calculate_ccsd_t_cbs.py](scripts/qc/5_calculate_ccsd_t_cbs.py)
Calculates interaction energy at the extrapolated CCSD(T)/CBS limit using PySCF.

**Arguments:**

| Flag              | Description                                 | Default                |
|-------------------|---------------------------------------------|------------------------|
| `-i` / `--input`  | Input XYZ file path                         |                        |
| `-c` / `--charge` | Charge                                      | Inferred from filename |
| `-o` / `--output` | Output file path with the energy in Hartree |                        |

#### [6_calculate_sapt_gold.py](scripts/qc/6_calculate_sapt_gold.py)
Calculates SAPT "gold standard" energy components using Psi4.

**Arguments:**

| Flag               | Description           | Default                |
|--------------------|-----------------------|------------------------|
| `-i` / `--input`   | Input XYZ file path   |                        |
| `-c` / `--charge`  | Charge                | Inferred from filename |
| `-o` / `--output`  | Psi4 output file path |                        |
| `-w` / `--workers` | Number of workers     | Number of CPU cores    |
| `-m` / `--memory`  | Total memory          |                        |

#### [7_calculate_sapt_silver.py](scripts/qc/7_calculate_sapt_silver.py)
Calculates SAPT "silver standard" energy components using Psi4.

**Arguments:**

| Flag               | Description           | Default                |
|--------------------|-----------------------|------------------------|
| `-i` / `--input`   | Input XYZ file path   |                        |
| `-c` / `--charge`  | Charge                | Inferred from filename |
| `-o` / `--output`  | Psi4 output file path |                        |
| `-w` / `--workers` | Number of workers     | Number of CPU cores    |
| `-m` / `--memory`  | Total memory          |                        |

## Dependencies

[Biotite](https://github.com/biotite-dev/biotite/)

[Dacite](https://github.com/konradhalas/dacite)

[More Itertools](https://github.com/more-itertools/more-itertools)

[numpy](https://github.com/numpy/numpy)

[OPI](https://github.com/faccts/opi)

[PDBe Arpeggio](https://github.com/PDBeurope/arpeggio)

[psi4](https://github.com/psi4/psi4/)

[pydash](https://github.com/dgilland/pydash)

[pySCF](https://github.com/pyscf/pyscf)

[RDKit](https://github.com/rdkit/rdkit/)

[Requests](https://github.com/psf/requests)

[SymPy](https://github.com/sympy/sympy)

[tqdm](https://github.com/tqdm/tqdm)
from dataclasses import dataclass
from itertools import product
from pathlib import Path

from pydash import join, py_
from pyscf import gto, scf

from qc_common import get_nucleobase_and_sidechain, get_input_path, get_charge, get_output_path


@dataclass
class GeometryDetails:
    name: str
    xyz: str
    charge: int


def calculate_interaction_energy(input_path: Path, sidechain_charge: int, output_path: Path):
    """Calculates interaction energy at the extrapolated CCSD(T)/CBS limit using PySCF."""

    nucleobase, sidechain = get_nucleobase_and_sidechain(input_path, sidechain_charge)
    nucleobase_xyz, nucleobase_ghost_xyz, sidechain_xyz, sidechain_ghost_xyz = [
        py_(s).map(lambda a: join([a.element if not is_ghost else f'X-{a.element}', *a.coord.tolist()], '\t')).join(
            '\n').value() for (s, is_ghost) in product([nucleobase, sidechain], [False, True])]

    dimer_details = GeometryDetails('dimer', '\n'.join([sidechain_xyz, nucleobase_xyz]), sidechain_charge)
    sidechain_details = GeometryDetails('sidechain', '\n'.join([sidechain_xyz, nucleobase_ghost_xyz]), sidechain_charge)
    nucleobase_details = GeometryDetails('nucleobase', '\n'.join([sidechain_ghost_xyz, nucleobase_xyz]), 0)

    e = {}

    for d in [dimer_details, nucleobase_details, sidechain_details]:
        mol_aQZ = gto.M(atom=d.xyz, basis='aug-cc-pvqz', charge=d.charge)
        hf_aQZ = mol_aQZ.HF().apply(scf.addons.remove_linear_dep_).density_fit().run()
        e_hf_aQZ = hf_aQZ.e_tot

        mp2_aQZ = hf_aQZ.MP2().set_frozen()
        mp2_aQZ.kernel(with_t2=False)
        e_mp2_corr_aQZ = mp2_aQZ.e_corr

        mol_aTZ = gto.M(atom=d.xyz, basis='aug-cc-pvtz', charge=d.charge)
        mp2 = mol_aTZ.HF().apply(scf.addons.remove_linear_dep_).density_fit().run().MP2().set_frozen()
        mp2.kernel(with_t2=False)
        e_mp2_corr_aTZ = mp2.e_corr

        e_mp2_corr_cbs = (64 * e_mp2_corr_aQZ - 27 * e_mp2_corr_aTZ) / 37  # beta = 3

        mol_aDZ = gto.M(atom=d.xyz, basis='aug-cc-pvdz', charge=d.charge)
        hf_aDZ = mol_aDZ.HF().apply(scf.addons.remove_linear_dep_).density_fit().run()
        mp2 = hf_aDZ.MP2().set_frozen()
        mp2.kernel(with_t2=False)
        e_mp2_aDZ = mp2.e_tot

        ccsd_aDZ = hf_aDZ.CCSD().set_frozen().run()

        e_ccsd_t_aDZ = ccsd_aDZ.e_tot + ccsd_aDZ.ccsd_t()
        delta_e_ccsd_t_mp2_aDZ = e_ccsd_t_aDZ - e_mp2_aDZ

        e_tot_cbs = e_hf_aQZ + e_mp2_corr_cbs + delta_e_ccsd_t_mp2_aDZ

        e[d.name] = e_tot_cbs

    interaction_energy = e['dimer'] - e['sidechain'] - e['nucleobase']
    output_path.write_text(str(interaction_energy))

    print(f'Energy {interaction_energy} Eh written to {output_path}')


if __name__ == '__main__':
    input_path = get_input_path()
    charge = get_charge()
    output_path = get_output_path()
    calculate_interaction_energy(input_path, charge, output_path)

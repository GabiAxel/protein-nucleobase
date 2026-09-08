from math import degrees
from pathlib import Path

import numpy as np
from biotite.structure import distance, Atom
from pydash import py_, min_by, map_
from sympy import Point3D, Line3D, Plane

from qc_common import get_nucleobase_and_sidechain, get_input_path, get_charge

content = [['binding_mode', 'buckle', 'propeller']]


def get_plane_angles(input_path: Path, charge: int):

    nucleobase, sidechain = get_nucleobase_and_sidechain(input_path, charge)

    sidechain_heavy = sidechain[sidechain.element != 'H']

    nucleobase_bound_atoms = py_(nucleobase[nucleobase.element != 'H']).filter(lambda a: any(distance(a, sa) < 3.5 for sa in sidechain_heavy)).take(2).value()
    if len(nucleobase_bound_atoms) == 2:
        sidechain_bound_atoms = py_(nucleobase_bound_atoms).map(lambda nucleobase_atom: min_by(sidechain_heavy[np.isin(sidechain_heavy.element, ['N', 'O'])], lambda sidechain_atom: distance(sidechain_atom, nucleobase_atom))).uniq_by('atom_name').value()
    else:
        sidechain_bound_atoms = py_(sidechain_heavy[np.isin(sidechain_heavy.element, ['N', 'O'])]).sort_by(lambda a: distance(a, nucleobase_bound_atoms[0])).take(2).value()

    closest_atom: Atom = py_(sidechain_heavy).without(*sidechain_bound_atoms).min_by(lambda a: distance(a, sidechain_bound_atoms[0])).value()
    if len(sidechain_bound_atoms) == 1:
        sidechain_connector_atoms = [closest_atom, sidechain_bound_atoms[0]]
    else:
        far_atom = py_(sidechain_heavy).reject(lambda a: a.atom_name in map_([closest_atom, *sidechain_bound_atoms], 'atom_name')).min_by(lambda a: distance(a, closest_atom)).value()
        sidechain_connector_atoms = [closest_atom, far_atom]

    sidechain_connector_line: Line3D = Line3D(*[Point3D(a.coord) for a in sidechain_connector_atoms])
    nucleobase_plane: Plane = Plane(*[Point3D(a.coord) for a in nucleobase[np.isin(nucleobase.element, ['N', 'O'])][:3]])
    buckle_angle = abs(degrees(float(nucleobase_plane.angle_between(sidechain_connector_line))))

    propeller_angle = None
    if len(sidechain_bound_atoms) == 2:
        sidechain_face_line = Line3D(*[Point3D(a.coord) for a in sidechain_bound_atoms])
        propeller_angle = abs(degrees(float(nucleobase_plane.angle_between(sidechain_face_line).evalf())))

    return buckle_angle, propeller_angle


if __name__ == '__main__':
    input_path = get_input_path()
    charge = get_charge()
    buckle_angle, propeller_angle = get_plane_angles(input_path, charge)
    print(f'Buckle: {buckle_angle}° ; Propeller: {propeller_angle}°')
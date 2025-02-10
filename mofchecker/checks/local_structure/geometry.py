# -*- coding: utf-8 -*-
"""Utilities for geometry operations."""
import math
from typing import Iterable

import numpy as np
from pymatgen.core import Site, Structure
from pymatgen.util.coord import get_angle

from ..utils.get_indices import is_metal, num_neighbor_H, non_metal_neighbors, non_metal_neighbor
from ..data.definitions import COVALENT_RADII


def rotation_matrix(axis, theta):
    """Return the rotation matrix associated with counterclockwise rotation about the given axis by theta radians.

    Stolen from https://stackoverflow.com/questions/6802577/rotation-of-3d-vector

    Args:
        axis (np.array): axis of rotation
        theta (float): angle of rotation in radians

    Returns:
        np.array: rotation matrix
    """
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)  # pylint:disable=invalid-name
    b, c, d = -axis * math.sin(theta / 2.0)  # pylint:disable=invalid-name
    aa, bb, cc, dd = a * a, b * b, c * c, d * d  # pylint:disable=invalid-name
    bc, ad, ac, ab, bd, cd = (  # pylint:disable=invalid-name
        b * c,
        a * d,
        a * c,
        a * b,
        b * d,
        c * d,
    )
    return np.array(
        [
            [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
            [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
            [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc],
        ]
    )


def _maximum_angle(angle):
    diff_to_180 = np.abs(180 - angle)
    return max([angle, diff_to_180])


def get_angle_between_site_and_neighbors(site: Site, neighbors: Iterable[Site]) -> float:
    """Get the minimum angle between a site and its two neighbors."""
    vec_1 = site.coords - neighbors[1].site.coords
    vec_2 = site.coords - neighbors[0].site.coords
    return get_angle(vec_1, vec_2)


def _guess_underbound_nitrogen_cn3(
    structure: Structure, site_index: int, neighbors: list, tolerance: int = 10
) -> bool:
    """Check if there is a nitrogen with three neighbors that likely misses some coordination.

    Args:
        structure (Structure): pymatgen Structure object
        site_index (int): index of the central site that is check
        neighbors (list): list of neighboring sites
        tolerance (int): Tolerance for angle checks in degree.
            Defaults to 10.

    Returns:
        bool: True if the nitrogen is likely missing some coordination partner
    """
    min_angle = get_angle_between_site_and_neighbors(structure[site_index], neighbors)
    any_metal = False
    for neighbor in neighbors:
        if is_metal(neighbor.site):
            any_metal = True

    num_h = 0
    for neighbor in neighbors:
        if str(neighbor.site.specie) == "H":
            num_h += 1

    if min_angle < 110 + tolerance:
        # let's only do this if one of the neighbors is a metal.
        # sometimes the M-N bond is so long that it isn't correctly recognized
        # obviously, this now won't detect missing H on a floating NH3
        # but this is probably a rare situation
        if any_metal and (num_h == 2):
#            print(site.index)
            return True

    return False


def _guess_underbound_nitrogen_cn2(
    structure: Structure,
    site_index: int,
    neighbors: list,
    connected_sites_a: list,
    connected_sites_b: list,
    tolerance: float = 25,
) -> bool:
    """Check if there is a nitrogen with CN 2 that probably misses some coordination.

    Args:
        structure (Structure): pymatgen Structure object
        site_index (int): Index of the site on which the check is performed
        neighbors (list): List of neighboring sites
        connected_sites_a (list): List of neighbor sites for first neighbor
        connected_sites_b (list): List of neighbor sites for second neighbor
        tolerance (float): Tolerance for angle checks in degree.
             Defaults to 25.

    Returns:
        bool: True if there is a nitrogen that likely misses some coordination.
    """
    def vector_angle_to_plane(a, b, c, d):
        a, b, c, d = map(np.array, [a, b, c, d])
        ab = b - a
        ac = c - a
        normal_vector = np.cross(ab, ac)
        dc = c - d
        dot_product = np.dot(dc, normal_vector)
        magnitude_dc = np.linalg.norm(dc)
        magnitude_normal = np.linalg.norm(normal_vector)
        angle_rad = np.arcsin(dot_product / (magnitude_dc * magnitude_normal))
        angle_deg = np.degrees(angle_rad)
        return abs(angle_deg)
    
    angle = get_angle_between_site_and_neighbors(structure[site_index], neighbors)
    num_h = 0
  #  print(angle)
    bond_lengths = np.array(
        [
            structure.get_distance(site_index, neighbors[0].index),
            structure.get_distance(site_index, neighbors[1].index),
        ]
    )
    expected_bond_lengths = np.array(
        [
            float(1.36 + COVALENT_RADII[str(neighbors[0].site.specie)] - 0.76),
            float(1.36 + COVALENT_RADII[str(neighbors[1].site.specie)] - 0.76),
        ]
    )
    if str(neighbors[0].site.specie) == 'H':
        num_h = num_h+1
    if str(neighbors[1].site.specie) == 'H':
        num_h = num_h+1
    
    if (np.abs(180 - angle) < tolerance) or (  # pylint: disable=no-else-return
        np.abs(0 - angle) < tolerance
    ):
        # sp hybridization if the nitrogen is linear
        # this could be a nitride or a nitrosyl
        # usually, there is nothing to worry about if this is the case
        return False
    else:
        # typically angle around 109.5 degree for sp3 hybridization
        # if we only have two neighbors but the nitrogen is likely
        # sp3 this is suspicious
        # to be sure we will check if it is planar (pyridine) or
        # not (piperazine) in the case the two neighbors are carbon
        # if neighbor_species == set(["C", "C"]):

#        dihedral_a = structure.get_dihedral(
#            neighbors[0].index,
#            site_index,
#            neighbors[1].index,
#            connected_sites_a[0].index,
#        )
#        dihedral_b = structure.get_dihedral(
#            neighbors[0].index,
#            site_index,
#            neighbors[1].index,
#            connected_sites_b[0].index,
#        )
#
#        dihedral_c = structure.get_dihedral(
#            connected_sites_b[0].index,
#            neighbors[0].index,
#            site_index,
#            neighbors[1].index,
#        )
#
#        dihedral_d = structure.get_dihedral(
#            connected_sites_a[0].index,
#            neighbors[0].index,
#            site_index,
#            neighbors[1].index,
#        )
        for neighbors_neighbor in connected_sites_a:
            if neighbors_neighbor.index == site_index:
                image_shift = structure[site_index].coords-neighbors_neighbor.site.coords
        for neighbors_neighbor in connected_sites_a:
            dihedral_a = 0
            if neighbors_neighbor.index != site_index:
                dihedral = vector_angle_to_plane(neighbors[1].site.coords, structure[site_index].coords, neighbors[0].site.coords, image_shift+neighbors_neighbor.site.coords)
                if dihedral > dihedral_a:
                    dihedral_a = dihedral
        for neighbors_neighbor in connected_sites_a:
            if neighbors_neighbor.index == site_index:
                image_shift = structure[site_index].coords-neighbors_neighbor.site.coords
        for neighbors_neighbor in connected_sites_b:
            dihedral_b = 0
            if neighbors_neighbor.index != site_index:
                dihedral = vector_angle_to_plane(neighbors[0].site.coords, structure[site_index].coords, neighbors[1].site.coords, image_shift+neighbors_neighbor.site.coords)
                if dihedral > dihedral_b:
                    dihedral_b = dihedral
        mean_dihedral = np.min(np.abs([dihedral_a, dihedral_b]))
        if (np.abs(mean_dihedral - 180) < tolerance) or (np.abs(mean_dihedral - 0) < tolerance):
            #print(bond_lengths, expected_bond_lengths)
            # we have to check if the one bond length is smaller than double bonds or two bonds smaller than aromatic partial bonds
            if num_h == 1 and (bond_lengths[0] < expected_bond_lengths[0] and bond_lengths[1] < expected_bond_lengths[1]):
                return False
            if num_h == 0 and (bond_lengths[0] < expected_bond_lengths[0] or bond_lengths[1] < expected_bond_lengths[1]):
                return False

            if num_h == 0 and (bond_lengths[0] < (expected_bond_lengths[0]+0.1) and bond_lengths[1] <(expected_bond_lengths[1]+0.1)):
                return False         
            return True
        else:
            
            return True
    # # larger angles should indicate sp2 hybridization
    # # one case where MOFs might have an issue with sp2
    # # is an NH2 group planar to the ring where one H is missing
    # # the heuristic we use to catch this is if one of the neighbors
    # # is H
    # if "H" in neighbor_species:
    #     return True
    #return False


def make_vec(start, end, length=None):
    """Create a vector based on a start and end position."""
    vector = end - start
    if length is not None:
        vector = vector / np.linalg.norm(vector) * length
    return vector


def add_sp_hydrogen(site, neighbors, length: float = 1):
    """Turn x#C -> x#C-H."""
    assert len(non_metal_neighbors(neighbors)) == 1
    vector = make_vec(site.coords, neighbors[0].site.coords, length)
    h_coords = site.coords + vector
    return h_coords


def add_sp2_hydrogen(site, neighbors, length: float = 1):
    """Convert x-C=z to x-CH-z."""
    assert len(neighbors) == 2

    vector0 = make_vec(neighbors[0].site.coords, site.coords)
    vector1 = make_vec(neighbors[1].site.coords, site.coords)
    summed = vector0 + vector1
    summed = summed / np.linalg.norm(summed) * length
    h_coords = site.coords + summed
    return h_coords

def add_O_hydrogen(structure: Structure,site, neighbors, length: float = 1):
    """Convert M-O-M to M-OH-M."""
    theta = np.linspace(0, 2 * np.pi, 72)
    phi = np.linspace(0, np.pi, 36)
    max_distance_sum = -np.inf
    best_x = None
    for t in theta:
        for p in phi:
            x = site.coords[0] + length * np.sin(p) * np.cos(t)
            y = site.coords[1] + length * np.sin(p) * np.sin(t)
            z = site.coords[2] + length * np.cos(p)
            atom = np.array([x,y,z])
            fra_coords = structure.lattice.get_fractional_coords(atom)
            dist_sum = 0
            for neighbor in neighbors:
                site_fra_coords = structure.lattice.get_fractional_coords(neighbor.site.coords)
                dist = structure.lattice.get_distance_and_image(fra_coords, site_fra_coords)[0]
                dist_sum += dist
            if dist_sum > max_distance_sum:
                max_distance_sum = dist_sum
                best_x = atom
    return best_x
#def add_x_atom(site, neighbors, length: float = 2):
#    """Convert O-X-O to O-XZ-O."""
#    n = len(neighbors)+1
#    summed = np.array([0,0,0])
#    for i in range(n):
#        vector = make_vec(neighbors[i-1].site.coords, site.coords)
#        vector /= np.linalg.norm(vector)
#        summed = summed + vector
#    summed = summed / np.linalg.norm(summed) * length
#    z_coords = site.coords + summed

#    return z_coords

def add_x_atom(structure: Structure, site, neighbors, length: float = 2.3):
    """Convert O-X-O to O-XZ-O."""
    theta = np.linspace(0, 2 * np.pi, 72)
    phi = np.linspace(0, np.pi, 36)
    max_distance_sum = -np.inf
    best_x = None
    for t in theta:
        for p in phi:
            x = site.coords[0] + length * np.sin(p) * np.cos(t)
            y = site.coords[1] + length * np.sin(p) * np.sin(t)
            z = site.coords[2] + length * np.cos(p)
            atom = np.array([x,y,z])
            fra_coords = structure.lattice.get_fractional_coords(atom)
            dist_sum = 0
            for neighbor in neighbors:
                site_fra_coords = structure.lattice.get_fractional_coords(neighbor.site.coords)
                dist = structure.lattice.get_distance_and_image(fra_coords, site_fra_coords)[0]
                dist_sum += dist
            if dist_sum > max_distance_sum:
                max_distance_sum = dist_sum
                best_x = atom
    return best_x
def add_x2_atom(site, neighbors, length: float = 2):
    """Convert O-X-O to O-XZZ-O."""
    n = len(neighbors)+1
    summed = np.array([0,0,0])
    for i in range(n):
        vector = make_vec(neighbors[i-1].site.coords, site.coords)
        vector /= np.linalg.norm(vector)
        summed = summed + vector
    summed = summed / np.linalg.norm(summed) * length
    z2_coords = site.coords - summed

    return z2_coords

def add_methylene_hydrogens(site, neighbors, length: float = 1):
    """Convert x-C-z to z-CH2-z."""
    assert len(neighbors) == 2
    vector = make_vec(neighbors[0].site.coords, site.coords)
    vector1 = make_vec(neighbors[1].site.coords, site.coords)
    summed = vector + vector1

    normal = np.cross(vector, vector1)

    hydrogen_1 = summed + normal
    hydrogen_1 /= np.linalg.norm(hydrogen_1) * length

    hydrogen_2 = summed - normal
    hydrogen_2 /= np.linalg.norm(hydrogen_2) * length

    hydrogen_1 = site.coords + hydrogen_1
    hydrogen_2 = site.coords + hydrogen_2
    return [hydrogen_1, hydrogen_2]


def get_some_orthorgonal_vector(vector):
    """Generate a orthogonal vector by cross product with a random vector.

    Will fail if the randly chosen vector is parallel to the input vector.

    Args:
        vector (np.array): Input vector

    Returns:
        np.array: Orthogonal vector
    """
    rand_vec = np.array([np.random.rand(), np.random.rand(), np.random.rand()])
    new_vec = np.cross(rand_vec, vector)
    new_vec /= np.linalg.norm(new_vec)
    return new_vec


def add_sp3_hydrogen(site, neighbors, length: float = 1):
    """Turn H2N-M --> H3N-M."""
    vector = make_vec(neighbors[0].site.coords, site.coords)
    vector1 = make_vec(neighbors[1].site.coords, site.coords)

    vector /= np.linalg.norm(vector)
    vector1 /= np.linalg.norm(vector1)
    summed = vector + vector1
    summed = summed / np.linalg.norm(summed) * length
    new_position = site.coords + summed

    return new_position


def add_sp3_hydrogens_on_cn1(site, neighbors, length: float = 1):
    """Add hydrogens to a SP3 carbon with coordination number 1.

    We make a simple geometric construction based on a triangle
    which normal vector is the vector from the current neighbor and the central atom.
    The cross product then gives us the next vector which we then only need to rotate
    twice around 120 degrees.

    Args:
        site (pymatgen.core.structure.Site): Site to add hydrogens to
        neighbors (list): List of neighbors
        length (float): Length of the bond. Defaults to 1.

    Returns:
        List[np.array]: List of new positions
    """
    vector = make_vec(neighbors[0].site.coords, site.coords)

    center = site.coords + vector / np.linalg.norm(vector) * length * np.cos(np.deg2rad(71))
    orthogonal_vector = get_some_orthorgonal_vector(vector) * length * np.sin(np.deg2rad(71))

    second_vec = np.dot(rotation_matrix(vector, np.deg2rad(120)), orthogonal_vector)
    third_vec = np.dot(rotation_matrix(vector, np.deg2rad(240)), orthogonal_vector)

    first_h = center + orthogonal_vector
    second_h = center + second_vec
    third_h = center + third_vec

    return [first_h, second_h, third_h]
def add_2_hydrogens_on_cn1(site, neighbors, length: float = 1):
    vector = make_vec(non_metal_neighbor(neighbors).site.coords, site.coords)

    center = site.coords + vector / np.linalg.norm(vector) * length * np.cos(np.deg2rad(71))
    orthogonal_vector = get_some_orthorgonal_vector(vector) * length * np.sin(np.deg2rad(71))

    second_vec = np.dot(rotation_matrix(vector, np.deg2rad(120)), orthogonal_vector)

    first_h = center + orthogonal_vector
    second_h = center + second_vec

    return [first_h, second_h]
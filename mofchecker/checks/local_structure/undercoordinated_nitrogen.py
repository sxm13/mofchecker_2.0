# -*- coding: utf-8 -*-
"""Check for undercoordinated nitrogens."""
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.util.coord import get_angle
from mofchecker.types import StructureIStructureType
import numpy as np
from .base_missing_check import BaseMissingCheck
from .geometry import (
    _guess_underbound_nitrogen_cn2,
    _guess_underbound_nitrogen_cn3,
    add_sp2_hydrogen,
    add_sp3_hydrogen,
    add_sp_hydrogen,
    get_angle_between_site_and_neighbors,
    add_2_hydrogens_on_cn1
)
from ..utils.get_indices import get_n_indices, num_neighbor_metal, non_metal_neighbors


class UnderCoordinatedNitrogenCheck(BaseMissingCheck):
    """Check for undercoordinated nitrogens."""

    def __init__(self, structure: StructureIStructureType, structure_graph: StructureGraph):
        """Initialize the check.

        Args:
            structure (StructureIStructureType): The structure to check
            structure_graph (StructureGraph): The structure graph of the structure
        """
        self.structure = structure
        self.n_indices = get_n_indices(self.structure)
        self.structure_graph = structure_graph

    @property
    def name(self):
        """Return the name of the check."""
        return "Undercoordinated nitrogen"

    @property
    def description(self):
        """Return a description of the check."""
        return "Checks, using geometric heuristics,\
             if there are any nitrogens that are likely undercoordinated."

    def _run_check(self):
        undercoordinated_nitrogens, positions = self._get_undercoordinated_nitrogens()
        return (
            len(undercoordinated_nitrogens) == 0,
            undercoordinated_nitrogens,
            positions,
        )

    def _get_undercoordinated_nitrogens(self, tolerance: float = 35):
        """Capture missing hydrogens on nitrogen groups using heuristics.

        Args:
            tolerance (float): angle tolerance for the check

        Returns:
            List[int], np.typing.ArrayLike: list of undercoordinated nitrogens and candidate positions
        """       
        undercoordinated_nitrogens = []
        h_positions = []
        for site_index in self.n_indices:
            cn = self.get_cn(site_index)  # pylint:disable=invalid-name
            neighbors = self.get_connected_sites(site_index)
            cm = num_neighbor_metal(neighbors)
            if cn-cm == 1:
                # this is suspicous, but it also might a CN which is perfectly fine.
                # to check this, we first see if the neighbor is carbon
                # and then what its coordination number is. If it is greater than 2
                # then we likely do not have a CN for which the carbon should be a
                # linear sp one
            #    if (self.get_cn(neighbors[0].index) > 2) and not neighbors[0].site.specie.is_metal:
                '''check if this terminal N has triple bond with neighbor atom, normal NN,CN triple bond is smaller than 1.25A and sp1 should be linear'''
                neighbors = non_metal_neighbors(neighbors)
                if (str(neighbors[0].site.specie) == 'C' or str(neighbors[0].site.specie) == 'N') and self.structure.get_distance(neighbors[0].index,site_index) < 1.25:
                    continue
                elif (str(neighbors[0].site.specie) == 'C' or str(neighbors[0].site.specie) == 'N') and self.structure.get_distance(neighbors[0].index,site_index) > 1.35:
                    undercoordinated_nitrogens.append(site_index)
                    for h_site in add_2_hydrogens_on_cn1(self.structure[site_index], neighbors):
                        h_positions.append(h_site)
                else:
                    undercoordinated_nitrogens.append(site_index)
                    h_positions.append(add_sp_hydrogen(self.structure[site_index], neighbors))
            elif cn-cm == 2:
                neighbors = non_metal_neighbors(neighbors)
                undercoordinated_nitrogen = _guess_underbound_nitrogen_cn2(
                    self.structure,
                    site_index,
                    neighbors,
                    self.get_connected_sites(neighbors[0].index),
                    self.get_connected_sites(neighbors[1].index),
                    tolerance,
                )
                if undercoordinated_nitrogen:
                    h_positions.append(add_sp2_hydrogen(self.structure[site_index], neighbors))
 #           elif cn == 3:
 #               undercoordinated_nitrogen = _guess_underbound_nitrogen_cn3(
 #                   self.structure, site_index, neighbors, tolerance
 #               )
 #               if undercoordinated_nitrogen:
 #                   undercoordinated_nitrogens.append(site_index)
  #                  h_positions.append(add_sp3_hydrogen(self.structure[site_index], neighbors))
        return undercoordinated_nitrogens, h_positions

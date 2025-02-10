# -*- coding: utf-8 -*-
from pymatgen.analysis.graphs import StructureGraph
from mofchecker.types import StructureIStructureType
from .base_coordination_check import BaseCoordinationCheck
from .base_missing_check import BaseMissingCheck
from ..utils.get_indices import get_x_indices, get_x2_indices, _get_all_indices,  non_metal_neighbor, get_metal_indices
from .geometry import  add_x_atom, add_x2_atom
from ..oms import MOFOMS

class adding_linker(BaseMissingCheck):

    def __init__(self, structure: StructureIStructureType, structure_graph: StructureGraph):
        """Initialize a function for connecting open metal site and missing charged linker.

        Args:
            structure (StructureIStructureType): Structure to check.
            structure_graph (StructureGraph): StructureGraph of the structure.
        """
        
        self.structure = structure
        self.x_indices = get_x_indices(self.structure)
        self.metal_indices = get_metal_indices(self.structure)
        self.structure_graph = structure_graph
        self.graph = structure_graph
        self.MOFOMS = MOFOMS.from_mofchecker(self)

    @property
    def name(self):
        """Return the name of the check."""
        return "Add missing linkers"

    @property
    def description(self):
        """Return a description of the check."""
        return "Return where to add missing linkers."

    def _run_check(self):
        (number, z_coord,) = self._get_z_position()
        return (len(number) == 0, number, z_coord,)
        

    
    def _get_z_position(self):
        """return the coordinates of first atom z added to X"""
        z_positions = []
        x_sum = []
        oms_indices = self.MOFOMS.flagged_indices
        for site_index in oms_indices:
            neighbors = self.get_connected_sites(site_index)
            z_positions.append(add_x_atom(self.structure, self.structure[site_index], neighbors))
            x_sum.append(site_index)        
#        for site_index in self.x_indices:
#            neighbors = self.get_connected_sites(site_index)
#            z_positions.append(add_x_atom(self.structure, self.structure[site_index], neighbors))
#            x_sum.append(site_index)
#        for site_index in self.x2_indices:
 #           neighbors = self.get_connected_sites(site_index)
 #           z_positions.append(add_x_atom(self.structure[site_index], neighbors))
 #           z_positions.append(add_x2_atom(self.structure[site_index], neighbors))
 #           x_sum.append(site_index)
 #           x_sum.append(site_index)
        return x_sum, z_positions
    

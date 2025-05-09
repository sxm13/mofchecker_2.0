# -*- coding: utf-8 -*-
"""Check the positive charge from linkers like N with CN4 and O with CN3."""
from pymatgen.analysis.graphs import StructureGraph
from structuregraph_helpers.create import construct_clean_graph
from mofchecker.checker_types import StructureIStructureType
import networkx as nx
from .base_coordination_check import BaseCoordinationCheck
from ..utils.get_indices import num_neighbor_metal, get_n_indices, get_o_indices, get_ge_indices, get_sb_indices, is_metal


class Positive_charge_Check(BaseCoordinationCheck):

    def __init__(self, structure: StructureIStructureType, structure_graph: StructureGraph):
        """Initialize a new PositivechargeCheck.

        Args:
            structure (StructureIStructureType): Structure to check.
            structure_graph (StructureGraph): StructureGraph of the structure.
        """
        self.structure = structure
        self.n_indices = get_n_indices(self.structure)
        self.o_indices = get_o_indices(self.structure)
        self.ge_indices = get_ge_indices(self.structure)
        self.sb_indices = get_sb_indices(self.structure)
        self.structure_graph = structure_graph

    @property
    def name(self):
        """Return the name of the check."""
        return "Positive charge from linkers"

    @property
    def description(self):
        """Return a description of the check."""
        return "Check the positive charge from the linkers."

    def _run_check(self):
        nitrogen_oxygen = self._get_overcoordinated_nitrogen()
        return len(nitrogen_oxygen) == 0, nitrogen_oxygen

    def _get_overcoordinated_nitrogen(self):
        """Check for all N if CN=4 and all O if CN = 3, ignore metal bonds."""
        overcoordinated_nitrogen_oxygen = []
        N_sum = []
        N_jump = []
        cycles = []
        nx_graph = construct_clean_graph(self.structure_graph)
        simple_cycles = nx.simple_cycles(nx_graph, length_bound=16)
        for cycle in simple_cycles:
            cycles.append(cycle)
        for site_index in self.n_indices:
            if site_index not in N_jump:
                cn = self.get_cn(site_index)  # pylint:disable=invalid-name
                neighbors = self.get_connected_sites(site_index)             
                cm = num_neighbor_metal(neighbors)
                if cn-cm == 4:
                    overcoordinated_nitrogen_oxygen.append(site_index)
                else:
                    for cycle in cycles:
                        N_possible_jump = []
                        if len(cycle) == 16 and site_index in cycle and all( not is_metal(self.structure[ring]) for ring in cycle):
                            num_n = 0
                            num_c = 0
                            num_n3 = 0
                            for ring in cycle:
                                if str(self.structure[ring].specie) == 'N':
                                    N_possible_jump.append(ring)
                                    num_n += 1
                                    n_cn = self.get_cn(ring)
                                    n_cm = num_neighbor_metal(self.get_connected_sites(ring))
                                    if n_cn-n_cm == 3:
                                        num_n3 += 1
                                if str(self.structure[ring].specie) == 'C':
                                    num_c += 1
                            if (num_n == 4) and (num_n3 == 4) and (num_c == 12):
                                overcoordinated_nitrogen_oxygen.append(site_index)
                                overcoordinated_nitrogen_oxygen.append(site_index)
                            if (num_n == 4) and (num_n3 == 3) and (num_c == 12):
                                overcoordinated_nitrogen_oxygen.append(site_index)
                            if (num_n == 4) and (num_c == 12):
                                for n in N_possible_jump:
                                    N_jump.append(n)
                                break
        
        for site_index in self.o_indices:
            cn = self.get_cn(site_index)
            neighbors = self.get_connected_sites(site_index)             
            cm = num_neighbor_metal(neighbors)  # pylint:disable=invalid-name
            if cn-cm == 3:
                overcoordinated_nitrogen_oxygen.append(site_index)
        for site_index in self.sb_indices:
            overcoordinated_nitrogen_oxygen.append(site_index)
            overcoordinated_nitrogen_oxygen.append(site_index)
            overcoordinated_nitrogen_oxygen.append(site_index)
        for site_index in self.ge_indices:
            overcoordinated_nitrogen_oxygen.append(site_index)
            overcoordinated_nitrogen_oxygen.append(site_index)
            overcoordinated_nitrogen_oxygen.append(site_index)
            overcoordinated_nitrogen_oxygen.append(site_index)
        return overcoordinated_nitrogen_oxygen

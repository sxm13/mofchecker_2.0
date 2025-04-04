# -*- coding: utf-8 -*-
"""Check if there are fused rings with 5-ring N which is possible charged and have to be checked manaully."""
from pymatgen.analysis.graphs import StructureGraph
from structuregraph_helpers.create import construct_clean_graph
from mofchecker.checker_types import StructureIStructureType
from collections import Counter
from .base_coordination_check import BaseCoordinationCheck
from ..utils.get_indices import get_n_indices, get_c_indices,num_neighbor_metal, non_metal_neighbors,num_neighbor_H, non_H_neighbors, is_metal
from .geometry import _guess_underbound_nitrogen_cn2
import networkx as nx

class Fusedring_Check(BaseCoordinationCheck):


    def __init__(self, structure: StructureIStructureType, structure_graph: StructureGraph):
        """Initialize the fusedring check.

        Args:
            structure (StructureIStructureType): The structure to check.
            structure_graph (StructureGraph): The structure graph to use for the check.
        """
        self.structure = structure
        self.c_indices = get_c_indices(self.structure)
        self.n_indices = get_n_indices(self.structure)
        self.structure_graph = structure_graph

    @property
    def name(self):
        """Return the name of the check."""
        return "Possible charged fused ring with N"

    @property
    def description(self):
        """Return a description of the check."""
        return "Checks if there are fused ring with 5-ring which contains N, charge should be checked manually."

    def _run_check(self):
        fused_ring = self._get_fused_ring()
        return len(fused_ring) == 0, fused_ring
    def _get_fused_ring(self):
        N_sum = []
        N_jump = []
        cycles = []
        nx_graph = construct_clean_graph(self.structure_graph)
        simple_cycles = nx.simple_cycles(nx_graph, length_bound=16)
        for cycle in simple_cycles:
            cycles.append(cycle)
        for site_index in self.n_indices:
            if site_index not in N_jump:
                cn = self.get_cn(site_index)
                neighbors = self.get_connected_sites(site_index)                
                cm = num_neighbor_metal(neighbors)
                if cn-cm==2:
                    N_5_ring = False
                    non_metals = non_metal_neighbors(neighbors)
                    neighbors_index = [site_index]
                    for neighbor in non_metals:
                        neighbors_index.append(neighbor.index)
                    for cycle in cycles:
                        if len(cycle) == 5 and site_index in cycle and all( not is_metal(self.structure[ring]) for ring in cycle):
                            N_5_ring = True
                            break
                    if N_5_ring:
                        for cycle in cycles:
                            if (len(cycle) < 10 and len(cycle) > 5)and all(index in cycle for index in neighbors_index) and all( not is_metal(self.structure[ring]) for ring in cycle):
                                N_sum.append(site_index)
        return N_sum

                                   
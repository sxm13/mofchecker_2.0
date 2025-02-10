# -*- coding: utf-8 -*-
from pymatgen.analysis.graphs import StructureGraph
from mofchecker.types import StructureIStructureType
from .base_coordination_check import BaseCoordinationCheck
from .base_missing_check import BaseMissingCheck
from ..utils.get_indices import get_o_indices, num_neighbor_metal, _get_all_indices,  non_metal_neighbor
from .geometry import  add_O_hydrogen

class O_site_adding_hydrogen(BaseMissingCheck):

    def __init__(self, structure: StructureIStructureType, structure_graph: StructureGraph):
        """Initialize a new Check for missing H atom site.

        Args:
            structure (StructureIStructureType): Structure to check.
            structure_graph (StructureGraph): StructureGraph of the structure.
        """
        self.structure = structure
        self.o_indices = get_o_indices(self.structure)
        self.all_indices = _get_all_indices(self.structure)
        self.structure_graph = structure_graph

    @property
    def name(self):
        """Return the name of the check."""
        return "Add missing hydrogen"

    @property
    def description(self):
        """Return a description of the check."""
        return "Return the coordniates of missing hydrogen."

    def _run_check(self):
        (number, hydrogen_coord,) = self._get_O_charge()
        return (len(number) == 0, number, hydrogen_coord,)
        

    
    def _get_O_charge(self):
        """Only check missing H at O site"""
        O2_sum = []
        O3_sum = []
        O1_sum = []
        O_sum = []
        h_positions = []
        h2_positions = []
        h3_positions = []
        h1_positions = []
        O_jump = []
        for site_index in self.o_indices:
            if site_index not in O_jump:
                cn = self.get_cn(site_index)
                neighbors = self.get_connected_sites(site_index)                
                cm = num_neighbor_metal(neighbors)
                """check for single O atom"""
                if cn-cm==0:
                    O2_sum.append(site_index)
                    h2_positions.append(add_O_hydrogen(self.structure, self.structure[site_index], neighbors))
                    """check for netural O"""
                elif cn-cm==2:
                    continue

                elif cn-cm==1:
                    non_metal = non_metal_neighbor(neighbors)
                    """check for single OH group""" 
                    if str(non_metal.site.specie) == 'H':
                        O1_sum.append(site_index)
                        h1_positions.append(add_O_hydrogen(self.structure, self.structure[site_index], neighbors))
                        """check for C-O connection"""
                    elif str(non_metal.site.specie) == 'C':
                        neighbor_C = non_metal.index
                        disC_O = self.structure.get_distance(neighbor_C,site_index)
                        C_neighbors = self.get_connected_sites(neighbor_C)
                        num_O = 0
                        n_O = 0
                        for C_neighbor in C_neighbors:
                                C_neighbor_site = C_neighbor.site
                                neighbor_cn = self.get_cn(C_neighbor.index)
                                neighbor_cm = num_neighbor_metal(self.get_connected_sites(C_neighbor.index))
                                if str(C_neighbor_site.specie) == 'O':
                                    n_O = n_O+1
                                if str(C_neighbor_site.specie) == 'O' and neighbor_cn-neighbor_cm ==1:
                                    num_O = num_O+1
                                    """only have to check one atom out of two O from carboxylic group"""
                                    if C_neighbor.index != site_index:
                                        O_possible_jump = C_neighbor.index
                        """check for charged carboxylic group by finding one C atom coordinated with two O atoms and without H or R on O site"""
                        if num_O == 2:
                            O1_sum.append(site_index)
                            h1_positions.append(add_O_hydrogen(self.structure, self.structure[site_index], neighbors))
                            O_jump.append(O_possible_jump)
                        """check for alcohol or phenol by analysing C-O bond distance to predict if it is a single bond"""
                        if n_O ==1 and disC_O > 1.315:
                            O1_sum.append(site_index)
                            h1_positions.append(add_O_hydrogen(self.structure, self.structure[site_index], neighbors))

                        """check for N-O connection"""
                    elif str(non_metal.site.specie) == 'N':

                        neighbor_N = non_metal.index
                        N_neighbors = self.get_connected_sites(neighbor_N)
                        neighbor_cn = self.get_cn(neighbor_N)
                        num_O = 0
                        O_NOx = []
                        num_metal = num_neighbor_metal(self.get_connected_sites(neighbor_N))
                        for N_neighbor in N_neighbors:
                            if str(N_neighbor.site.specie) == 'O':
                                num_O = num_O+1
                                O_NOx.append(N_neighbor)

                        """check for charged nitrate group rather than nitro group"""

                        negative_O = 1-(neighbor_cn-num_metal-num_O)
                        if negative_O > 0:
                            for O_atom in O_NOx:
                                O_NOx_cn = self.get_cn(O_atom.index)
                                O_NOx_neighbors = self.get_connected_sites(O_atom.index)                
                                O_NOx_cm = num_neighbor_metal(O_NOx_neighbors)
                                if O_atom.index != site_index:
                                    O_jump.append(O_atom.index)
                                if O_NOx_cn-O_NOx_cm != 1:
                                    negative_O = negative_O-1
                            if negative_O==1:
                                O1_sum.append(site_index)
                                h1_positions.append(add_O_hydrogen(self.structure, self.structure[site_index], neighbors))
                            elif negative_O==0:
                                continue

                        """check for S-O connection"""
                    elif str(non_metal.site.specie) == 'S':                   
                        neighbor_S = non_metal.index
                        S_neighbors = self.get_connected_sites(neighbor_S)
                        neighbor_cn = self.get_cn(neighbor_S)
                        num_O = 0
                        O_SOx = []
                        num_metal = num_neighbor_metal(self.get_connected_sites(neighbor_S))
                        for S_neighbor in S_neighbors:
                            if str(S_neighbor.site.specie) == 'O':
                                num_O = num_O+1
                                O_SOx.append(S_neighbor)

                        """check for sulfate,sulfite,sulfonic group and sulfinic group SOx"""

                        negative_O = 2-(neighbor_cn-num_metal-num_O)
                        if negative_O > 0:
                            for O_atom in O_SOx:
                                O_SOx_cn = self.get_cn(O_atom.index)
                                O_SOx_neighbors = self.get_connected_sites(O_atom.index)                
                                O_SOx_cm = num_neighbor_metal(O_SOx_neighbors)
                                if O_atom.index != site_index:
                                    O_jump.append(O_atom.index)
                                if O_SOx_cn-O_SOx_cm != 1:
                                    negative_O = negative_O-1
                            if negative_O==2:
                                O2_sum.append(site_index)
                                h2_positions.append(add_O_hydrogen(self.structure, self.structure[site_index], neighbors))
                            elif negative_O==1:
                                O1_sum.append(site_index)
                                h1_positions.append(add_O_hydrogen(self.structure, self.structure[site_index], neighbors))
                            elif negative_O==0:
                                continue

                        """check for P-O connection"""
                    elif str(non_metal.site.specie) == 'P':
                        neighbor_P = non_metal.index
                        P_neighbors = self.get_connected_sites(neighbor_P)
                        neighbor_cn = self.get_cn(neighbor_P)
                        num_O = 0
                        O_POx = []
                        num_metal = num_neighbor_metal(self.get_connected_sites(neighbor_P))
                        for P_neighbor in P_neighbors:
                            if str(P_neighbor.site.specie) == 'O':
                                num_O = num_O+1
                                O_POx.append(P_neighbor)
                        """check for POx"""

                        negative_O = 3-(neighbor_cn-num_metal-num_O)

                        if negative_O > 0:
                            for O_atom in O_POx:
                                O_POx_cn = self.get_cn(O_atom.index)
                                O_POx_neighbors = self.get_connected_sites(O_atom.index)                
                                O_POx_cm = num_neighbor_metal(O_POx_neighbors)
                                if O_atom.index != site_index:
                                    O_jump.append(O_atom.index)
                                if O_POx_cn-O_POx_cm != 1:
                                    negative_O = negative_O-1

                            if negative_O==2:
                                O3_sum.append(site_index)
                                h3_positions.append(add_O_hydrogen(self.structure, self.structure[site_index], neighbors))
                            elif negative_O==1:
                                O1_sum.append(site_index)
                                h1_positions.append(add_O_hydrogen(self.structure, self.structure[site_index], neighbors))
                            elif negative_O==0:
                                continue
                            elif negative_O == 3:
                                O3_sum.append(site_index)  
                                h3_positions.append(add_O_hydrogen(self.structure, self.structure[site_index], neighbors))                              

                        """check for halogen-O connection like ClO4-, IO3-..."""
                    elif str(non_metal.site.specie) == 'Cl' or str(non_metal.site.specie) == 'Br' or str(non_metal.site.specie) == 'I':
                        neighbor_F = non_metal.index
                        F_neighbors = self.get_connected_sites(neighbor_F)
                        neighbor_cn = self.get_cn(neighbor_F)
                        num_O = 0
                        O_ClOx = []
                        num_metal = num_neighbor_metal(self.get_connected_sites(neighbor_F))
                        for F_neighbor in F_neighbors:
                            if str(F_neighbor.site.specie) == 'O':
                                num_O = num_O+1
                                O_ClOx.append(F_neighbor)

                        """check for charged group ClOx"""

                        negative_O = 1-(neighbor_cn-num_metal-num_O)
                        if negative_O > 0:
                            for O_atom in O_ClOx:
                                O_ClOx_cn = self.get_cn(O_atom.index)
                                O_ClOx_neighbors = self.get_connected_sites(O_atom.index)                
                                O_ClOx_cm = num_neighbor_metal(O_ClOx_neighbors)
                                if O_atom.index != site_index:
                                    O_jump.append(O_atom.index)
                                if O_ClOx_cn-O_ClOx_cm != 1:
                                    negative_O = negative_O-1
                            if negative_O==1:
                                O1_sum.append(site_index)
                                h1_positions.append(add_O_hydrogen(self.structure, self.structure[site_index], neighbors))
                            elif negative_O==0:
                                continue


                        """other non_metal elements connection"""
                    else:
                        O1_sum.append(site_index)
                        h1_positions.append(add_O_hydrogen(self.structure, self.structure[site_index], neighbors))
                else:
                    continue
        O_sum = O2_sum + O3_sum + O1_sum
        h_positions = h2_positions + h3_positions + h1_positions
        return O_sum, h_positions
    

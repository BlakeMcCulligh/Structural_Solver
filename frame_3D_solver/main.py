"""
Holds the 3D frame solver object and handels everything to do with the solving of 3D frames.
"""

import frame_3D_solver.helper_functions as hf
import numpy as np
import scipy.optimize as opt

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

class Frame3D:
    """
    Object that handels the solving of a 3D frame.
    """
    def __init__(self) -> None:
        """
        Initialises the 3D frame solver object.
        """

        # general Lists
        self.casses: np.ndarray | list = []

        # Node Lists
        self.nodes_cord: np.ndarray | list = [] # shape: (# Nodes, 3)
        self.nodes_loads: np.ndarray | list = [] # shape: (# Nodes, # Casses, 6: [Px, Py, Pz, Mx, My, Mz])
        self.nodes_support: np.ndarray | list = [] # shape: (# Nodes, 6)
        self.nodes_dof_unknown: np.ndarray | list = None
        self.nodes_dof_known: np.ndarray | list = None

        # Material Lists
        self.materials: np.ndarray | list = [] # shape: (# Materials, 5: [E, G, nu, rho, fy])

        # Member Lists
        self.members: np.ndarray | list = [] # shape: (# Members, 4: [i_node, j_node, material, set_cross_section])
        self.members_cross_section_props = [] # shape: (# Members, 4: [A, Iy, Iz, J])
        self.members_releases = [] # shape: (# Members, 12)
        self.members_dof = [] # shape: (# Members, 12)
        self.members_L: list or np.ndarray = [] # shape: (# Members)
        self.members_dof_unreleced: list or np.ndarray = []
        self.members_dof_releced: list or np.ndarray = []
        self.members_T: list or np.ndarray = []

        # shape: (# Members, # Casses, 2: [Case Index, 2: [# loads: locations, # loads: [6]]])
        self.members_point_loads = []
        # shape: (# Members, # Casses, 2: [Case Index, # Loads:[2: locations,6:loads]])
        self.members_dist_loads = []

        # preped lists
        self.point_loads = None
        self.dist_loads = None

    def AddNode(self, X: float, Y: float, Z: float) -> None:
        """
        Adds a node to the frame.

        :param X: float. x coordanite
        :param Y: float. y coordanite
        :param Z: float. z coordanite
        """

        self.nodes_cord.append([X, Y, Z])
        self.nodes_loads.append([])
        self.nodes_support.append([False, False, False, False, False, False])

    def AddMaterial(self, E: float, G: float, nu: float, rho: float, fy: float | None = None) -> None:
        """
        Adds a material to the frame.

        :param E: float. Youngs Modulus
        :param G: float. Shear Modulus
        :param nu: float. Poisson's Ratio
        :param rho: float. Density
        :param fy: float. Yeald Strangth
        """

        self.materials.append([E, G, nu, rho, fy])

    def AddMember(self, i_node: int, j_node: int, material_index: int, set_cross_section_props: bool, A: float,
                  Iy: float, Iz: float, J: float) -> None:
        """
        Adds a member to the frame.

        :param i_node: int. Index of the start node of the member.
        :param j_node: int. Index of the end node of the member.
        :param material_index: int. Index of the material of the member.
        :param set_cross_section_props: bool. Whether the member's cross-sections properties should be optimized when
                                        a cross-section optimizer is run.
        :param A: float. Cross-section area.
        :param Iy: float. Moment of Inerta in the y direction.
        :param Iz: float. Moment of Inerta in the z direction.
        :param J: float. Moment of Inerta in the z direction.
        """

        self.members.append([i_node, j_node, material_index, set_cross_section_props])
        self.members_cross_section_props.append([float(A), float(Iy), float(Iz), float(J)])
        self.members_releases.append([False, False, False, False, False, False,
                                      False, False, False, False, False, False])
        self.members_dof.append([])
        self.members_L.append([])
        self.members_dof_unreleced.append([])
        self.members_dof_releced.append([])
        self.members_T.append([])

        self.members_point_loads.append([])
        self.members_dist_loads.append([])

    def AddSupport(self, node_index: int, support_DX: bool = False, support_DY: bool = False, support_DZ: bool = False,
                   support_RX: bool = False, support_RY: bool = False, support_RZ: bool = False) -> None:
        """
        Adds a support to the 3D frame.

        :param node_index: int. Index of node to add supports to.
        :param support_DX: bool. Indicates if the x direction should be supported.
        :param support_DY: bool. Indicates if the y direction should be supported.
        :param support_DZ: bool. Indicates if the z direction should be supported.
        :param support_RX: bool. Indicates if the x rotation direction should be supported.
        :param support_RY: bool. Indicates if the y rotation direction should be supported.
        :param support_RZ: bool. Indicates if the z rotation direction should be supported.
        """

        self.nodes_support[node_index] = [support_DX, support_DY, support_DZ, support_RX, support_RY, support_RZ]

    def AddReleases(self, member_index: int, Dxi: bool = False, Dyi: bool = False, Dzi: bool = False,
                    Rxi: bool = False, Ryi: bool = False, Rzi: bool = False, Dxj: bool = False, Dyj: bool = False,
                    Dzj: bool = False, Rxj: bool = False, Ryj: bool = False, Rzj: bool = False) -> None:
        """
        Adds a member release to the frame.

        :param member_index: int. Index of the member to add releases to.
        :param Dxi: bool. Indicates if the i node in the x direction should be released.
        :param Dyi: bool. Indicates if the i node in the y direction should be released.
        :param Dzi: bool. Indicates if the i node in the z direction should be released.
        :param Rxi: bool. Indicates if the i node in the x rotation direction should be released.
        :param Ryi: bool. Indicates if the i node in the y rotation direction should be released.
        :param Rzi: bool. Indicates if the i node in the z rotation direction should be released.
        :param Dxj: bool. Indicates if the j node in the x direction should be released.
        :param Dyj: bool. Indicates if the j node in the y direction should be released.
        :param Dzj: bool. Indicates if the j node in the z rotation direction should be released.
        :param Rxj: bool. Indicates if the j node in the x rotation direction should be released.
        :param Ryj: bool. Indicates if the j node in the y rotation direction should be released.
        :param Rzj: bool. Indicates if the j node in the xz direction should be released.
        """

        self.members_releases[member_index] = [Dxi, Dyi, Dzi, Rxi, Ryi, Rzi, Dxj, Dyj, Dzj, Rxj, Ryj, Rzj]

    def AddNodeLoad(self, node_index: int, Px: float = 0, Py: float = 0, Pz: float = 0,
                    Mx: float = 0, My: float = 0, Mz: float = 0, case: int = 0) -> None:
        """
        Adds a node point load to the frame.

        :param node_index: int. Index of the node to add the load to.
        :param Px: float. Load in the x direction.
        :param Py: float. Load in the y direction.
        :param Pz: float. Load in the z direction.
        :param Mx: float. Moment in the x direction.
        :param My: float. Moment in the y direction.
        :param Mz: float. Moment in the z direction.
        :param case: int. Index of the case to add the load to.
        """

        case_found = False

        for load_case in self.nodes_loads[node_index]:

            caseINDEX, load = load_case[0], load_case[1]

            if caseINDEX == case:
                load[0] += Px
                load[1] += Py
                load[2] += Pz
                load[3] += Mx
                load[4] += My
                load[5] += Mz
                case_found = True

        if not case_found:

            self.nodes_loads[node_index].append([case, [Px, Py, Pz, Mx, My, Mz]])

            if case not in self.casses: self.casses.append(case)

    def AddMemberPointLoad(self, member_index: int, x: float, Px: float = 0, Py: float = 0, Pz: float = 0,
                           Mx: float = 0, My: float = 0, Mz: float = 0, case: int = 0) -> None:
        """
        Adds a member point load to the frame.

        :param member_index: int. Index of the member to add the load to.
        :param x: float. Location on the member for the load.
        :param Px: float. Load in the x direction.
        :param Py: float. Load in the y direction.
        :param Pz: float. Load in the z direction.
        :param Mx: float. Moment in the x direction.
        :param My: float. Moment in the y direction.
        :param Mz: float. Moment in the z direction.
        :param case: int. Index of the case to add the load to.
        """

        case_found = False

        for loads_case in self.members_point_loads[member_index]:

            case_index, loads = loads_case[0], loads_case[1]

            if case_index == case:

                case_found = True
                load_found = False
                loadsX, loadsMag = loads[0], loads[1]

                for i in range(len(loadsX)):

                    if loadsX == x:
                        loadsMag[i][0] += Px
                        loadsMag[i][1] += Py
                        loadsMag[i][2] += Pz
                        loadsMag[i][3] += Mx
                        loadsMag[i][4] += My
                        loadsMag[i][5] += Mz
                        load_found = True

                if not load_found:
                    loadsX.append(x)
                    loadsMag.append([Px, Py, Pz, Mx, My, Mz])

        if not case_found:

            self.members_point_loads[member_index].append([case, [[x], [[Px, Py, Pz, Mx, My, Mz]]]])

            if case not in self.casses: self.casses.append(case)

    def addMemberDistLoad(self, member_index: int, x1: float, x2: float, wx1: float = 0, wx2: float = 0,
                          wy1: float = 0, wy2: float = 0, wz1: float = 0, wz2: float = 0, case: int = 0) -> None:
        """
        Adds a member distribution load to the frame.

        :param member_index: int. Index of the member to add the load to.
        :param x1: float. Starting location of the distributed load.
        :param x2: float. Ending location of the distributed load.
        :param wx1: float. Load in the x direction at the start of the distributed load.
        :param wx2: float. Load in the x direction at the end of the distributed load.
        :param wy1: float. Load in the y direction at the start of the distributed load.
        :param wy2: float. Load in the y direction at the end of the distributed load.
        :param wz1: float. Load in the z direction at the start of the distributed load.
        :param wz2: float. Load in the z direction at the end of the distributed load.
        :param case: int. Index of the load case to add the load to.
        """

        case_found = False

        for loads_case in self.members_dist_loads[member_index]:

            case_index, loads = loads_case[0], loads_case[1]

            if case_index == case:
                case_found = True
                loads.append([[x1, x2], [wx1, wx2, wy1, wy2, wz1, wz2]])

        if not case_found:

            self.members_dist_loads[member_index].append([case, [[[x1, x2], [wx1, wx2, wy1, wy2, wz1, wz2]]]])

            if case not in self.casses:
                self.casses.append(case)

    def PreAnalysisLinear(self, log: bool = False) -> None:
        """
        Runs all analysises that can be run befor the optimizer.

        :param log: bool. If sub calculations should be printed to the consel. Used for debuging.
        """

        if log:
            print("--------------------------------------------------------------")
            print("---------------  Pre-analysis linear  ------------------------")
            print("--------------------------------------------------------------")

        self.materials = np.array(self.materials)
        self.nodes_support = np.array(self.nodes_support)
        self.members_releases = np.array(self.members_releases)
        self.members = np.array(self.members)
        self.nodes_cord = np.array(self.nodes_cord)

        self.nodes_dof_unknown, self.nodes_dof_known = hf.part_D(self.nodes_support)

        (self.members_dof, self.members_L, self.members_dof_unreleced, self.members_dof_releced, self.members_T,
         self.point_loads, self.dist_loads) = hf.prep_members(self.nodes_cord, self.members, self.members_releases,
                                                              self.members_point_loads, self.members_dist_loads,
                                                              len(self.casses))

        self.members_L = np.array(self.members_L)

        if log:
            print("nodes_dof_unknown: ", self.nodes_dof_unknown)
            print("nodes_dof_known: ", self.nodes_dof_known)
            print("members-DOF: ", self.members_dof)
            print("members_L: ", self.members_L)
            print("members_dof_unreleced: ", self.members_dof_unreleced)
            print("members_dof_releced: ", self.members_dof_releced)
            print("members_T: ", self.members_T)
            print("point_loads: ", self.point_loads)
            print("dist_loads: ", self.dist_loads)

    def AnalysisLinear(self, get_weight: bool = False, get_reactions: bool = False, get_internal_forces: bool = False,
                       log: bool = False) -> (np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray,
                                              np.ndarray, np.ndarray | None, np.ndarray | None, np.ndarray | None):
        """
        Runs a linear analysis of the frame. PreAnalysisLinear needs to be run befor this.

        :param get_weight: bool. If the weight of the members should be calculated
        :param get_reactions: bool. If the reaction should be calculated
        :param get_internal_forces: bool. If the internal forces should be calculated
        :param log: bool. If sub calculations should be printed to the consel. Used for debuging.
        :return:
            D: ndarray. Main node deflection array.
            DX: ndarray. Node deflection in the X direction.
            DY: ndarray. Node deflection in the Y direction.
            DZ: ndarray. Node deflection in the Z direction.
            RX: ndarray. Node rotation deflection in the X direction.
            RY: ndarray. Node rotation deflection in the Y direction.
            RZ: ndarray. Node rotation deflection in the Z direction.
            weight: ndarray or None. Weight of each member of the frame.
            reactions: ndarray or None. Reactions of the nodes to teh supports.
            internal_forces: ndarray or None. Maximum internal forces of each member.
        """

        if log:
            print("--------------------------------------------------------------")
            print("-------------------  Analysis Linear  ------------------------")
            print("--------------------------------------------------------------")

        if isinstance(self.members_L, np.ndarray) and isinstance(self.members_T, np.ndarray):
            num_n = len(self.nodes_cord)
            num_m = len(self.members)
            num_c = len(self.casses)

            self.members_cross_section_props = np.array(self.members_cross_section_props)

            k_local = hf.get_k_local_array(self.materials, self.members, self.members_cross_section_props,
                                           self.members_L, log)

            k11, k12, k21, k22 = hf.member_part_k_array(k_local, self.members_dof_unreleced, self.members_dof_releced,
                                                        num_m)
            if log:
                print("k11: ", k11)
                print("k12: ", k12)
                print("k21: ", k21)
                print("k22: ", k22)

            fer_unc_ARRAY = hf.get_member_fer_unc(self.members_L, self.point_loads, self.dist_loads, num_m, num_c)
            if log:
                print("fer_unc_ARRAY: ", fer_unc_ARRAY)

            fer1, fer2 = hf.member_part_fer(fer_unc_ARRAY, self.members_dof_unreleced, self.members_dof_releced,
                                            num_m, num_c)
            if log:
                print("fer1 dim: ", np.shape(fer1))
                print("fer1: ", fer1)
                print("fer2: ", fer2)

            fer_condensed = hf.get_fer(self.members_releases, k12, k22, fer1, fer2, num_m, num_c)
            if log:
                print("fer_condensed: ", fer_condensed)

            FER1, FER2 = hf.get_global_fixed_end_reaction_vector(self.nodes_cord, self.members_dof, self.members_T,
                                                                 self.nodes_dof_unknown, self.nodes_dof_known,
                                                                 fer_condensed, num_m, num_c)
            if log:
                print("FER1: ", FER1)
                print("FER2: ", FER2)

            P1, P2 = hf.get_parted_global_nodal_force_vector(self.nodes_loads, self.casses, self.nodes_dof_unknown,
                                                             self.nodes_dof_known, num_n)
            if log:
                print("P1: ", P1)
                print("P2: ", P2)

            k_global_members = hf.k_member_make_global(k_local, self.members_T)
            if log: print("k_global_members: ", k_global_members)

            K_global = hf.get_K_global(self.members_dof, k_global_members, num_n, num_m)
            if log: print("K_global: ", K_global)

            K11, K12, K21, K22 = hf.partition_K_gloabl(K_global, self.nodes_dof_unknown, self.nodes_dof_known)
            if log:
                print("K11: ", K11)
                print("K12: ", K12)
                print("K21: ", K21)
                print("K22: ", K22)

            D, DX, DY, DZ, RX, RY, RZ = hf.get_D(K11, K12, P1, FER1, self.nodes_dof_unknown, self.nodes_dof_known,
                                                 num_n, num_c, log)

            weight = None
            reactions = None
            internal_forces = None

            if get_weight:
                weight = hf.get_weight(self.materials, self.members, self.members_L, self.members_cross_section_props)

            if get_reactions or get_internal_forces:
                D_members = hf.get_member_direction_deflections(self, DX, DY, DZ, RX, RY, RZ, num_m, num_c)
                d = hf.get_d(self.members_T, D_members, num_m, num_c)
                f = hf.get_f(k_local, d, fer_condensed, num_m, num_c)

                if get_reactions:
                    F = hf.get_F(self.members_T, f, num_m, num_c)
                    reactions = hf.get_reactions(self.nodes_support, self.nodes_loads, self.members,
                                                 self.members_releases, F, num_c, num_m, num_n)

                if get_internal_forces:
                    abs_F, abs_M = hf.solve_internal_forces(self.members, self.members_L,
                                                            self.members_cross_section_props, self.materials,
                                                            self.casses, self.point_loads, self.dist_loads, f,
                                                            fer_unc_ARRAY, d, num_m, num_c)
                    internal_forces = [abs_F, abs_M]

            return D, DX, DY, DZ, RX, RY, RZ, weight, reactions, internal_forces
        else:
            raise Exception('Pre analysis has not been run. Aborting analysis.')

    def optimize(self, memberGroup: list, memberGroupType: list,  lowerBound: list | float, upperBound: list | float,
                 costFunction: str, getWeight: bool = False, getReactions: bool = False,
                 getInternalForces: bool = False, log: bool =False):
        """
        Finds the optimum cross-sections for members with their cross-sections not set.

        :param memberGroup: list. list of indices of member groups for non set members to be assigned to.
                            Must be length of non set members.
        :param memberGroupType: list. List of cross-section types for each member group to be assigned.
                                Must be length of number of member gorups.
        :param lowerBound: list or float. Lower bound on the optimization variables.
                           If list must be length of number of variables.
        :param upperBound: list or float. Upper bound on the optimization variables.
                           If list must be length of number of variables.
        :param costFunction: String. Cost function for optimization stored in a string.
        :param getWeight: bool. Weather the weighht of all the members is needed for the cost function.
        :param getReactions: bool. Weather the reactions are needed for the cost function.
        :param getInternalForces: bool. Weather the internal forces are needed for the cost function.
        :return: scipi optimization_results class: results of the optimization.
        """

        if log:
            print("--------------------------------------------------------------")
            print("-----------------  Global Optimization  ----------------------")
            print("--------------------------------------------------------------")

        hf.chack_inputs(self.members, memberGroup, memberGroupType)

        numVarables = hf.get_num_varables(memberGroupType)
        if log: print("num_varables: ", numVarables)

        self.PreAnalysisLinear(log=log)

        constants = [self, costFunction, memberGroup, memberGroupType, getWeight, getReactions, getInternalForces, log]

        bounds = hf.get_bounds(lowerBound, upperBound, numVarables)

        optimization_results = opt.shgo(hf.get_cost, bounds, args=[constants])

        if log: print("optimization results: ", optimization_results)

        return optimization_results


#     simple_beam = Frame3D()
#
#     simple_beam.AddNode(0, 0, 0)
#     simple_beam.AddNode(168, 0, 0)
#     simple_beam.AddNode(168, 5, 0)
#
#     simple_beam.AddMaterial(29000, 11200, 0.3, 2.836e-4)
#
#     simple_beam.AddMember(0, 1, 0, False, 20, 100, 150, 250)
#     simple_beam.AddMember(0, 2, 0, False, 20, 100, 150, 250)
#
#     simple_beam.AddSupport(0, True, True, True, True, False, False)
#     simple_beam.AddSupport(1, True, True, True, True, False, False)
#     # simple_beam.AddSupport(0, True, True, True, True, True, True)
#     # simple_beam.AddSupport(1, True, True, True, True, True, True)
#     # 'M1', 'Fy', -0.01, -0.01, 0, 168
#     # simple_beam.AddNodeLoad(0, Pz=1, case=0)
#     simple_beam.AddMemberPointLoad(0, 50, Pz=50, case=0)
#     # simple_beam.AddMemberPointLoad(0, 2, Pz=1, case=0)
#     # simple_beam.AddMemberPointLoad(0, 2, Pz=1, case=1)
#     # simple_beam.addMemberDistLoad(0,0,5,5,2,0,0)
#     # simple_beam.addMemberDistLoad(0, 0, 168, wy1 = -0.01, wy2 = -0.01, case=0)
#     # simple_beam.addMemberSelfWeight()
#     # simple_beam.addMemberSelfWeight(case=1)
#     simple_beam.PreAnalysisLinear(log=False)
#     #print(simple_beam.AnalysisLinear(get_weight=True,get_internal_forces=True))
#     simple_beam.AnalysisLinear(log=False)
#
#     member_group_type = ["Angle"]
#     member_group = [0, 0]
#     simple_beam.optimize(member_group, member_group_type, [1,1,0.1],[10,10,0.9], get_weight = True, get_reactions = False, get_internal_forces = False, log = False)


import helperFunctions as hf
import numpy as np

class Frame3D_T:
    def __init__(self) -> None:

        # general Lists
        self.casses = []

        # Node Lists
        self.nodes_cord = []
        self.nodes_loads = [] # [[load case 1, [Px, Py, Pz, Mx, My, Mz]], [load case 2 ...  ]]
        self.nodes_D = []
        self.nodes_Reactions = []
        self.nodes_Support = []
        self.nodes_D_Enforced = []

        # Material Lists
        self.materials = []

        # Member Lists
        self.members = []
        self.members_CrossSectionProps = []
        self.members_SegmentsZ = []
        self.members_SegmentsY = []
        self.members_SegmentsX = []
        self.members_Releases = []
        self.members_DOF = []
        self.members_L = []
        self.members_PartD_unreleced = []
        self.members_PartD_releced = []
        self.members_T = []
        self.members_K = []

        self.members_PointLoads = [] # [[load case 1, [ [x1, x2 ,...], [[Px1, Py1, Pz1, Mx1, My1, Mz1],[Px2, Py2, Pz2, Mx2, My2, Mz2],...]], [load case 2 ...  ]]
        self.members_DistLoads = [] # [[load case 1, [ [[x1_1,x2_1], [wz1_1,wz2_1,wy1_1,wy2_1]], [[x1_2,x2_2], [wz1_2,wz2_2,wy1_2,wy2_2]],...]], [load case 2 ... ]]
        self.members_SelfWeight = [] # [[load case 1, [factorX, faxtorY, factorZ]], [load case 2 ...]]

        # preped lists
        self.D_unknown = None
        self.D_known = None
        self.D_known_val = None

    def addNode(self, X: float, Y: float, Z: float):
        self.nodes_cord.append([X, Y, Z])
        self.nodes_loads.append([])
        self.nodes_D.append([])
        self.nodes_Reactions.append([])
        self.nodes_Support.append([False,False,False,False,False,False])
        self.nodes_D_Enforced.append([None,None,None,None,None,None])

    def addMaterial(self, E: float, G: float, nu: float, rho: float, fy: float | None = None):
        self.materials.append([E, G, nu, rho, fy])

    def addMember(self, i_node: int, j_node: int, material_index: int, setCrossSectionProps: bool, A: float, Iy: float, Iz: float, J: float):
        self.members.append([i_node, j_node, material_index, setCrossSectionProps])
        self.members_CrossSectionProps.append([A, Iy, Iz, J])
        self.members_SegmentsZ.append([])
        self.members_SegmentsY.append([])
        self.members_SegmentsX.append([])
        self.members_Releases.append([False, False, False, False, False, False, False, False, False, False, False, False])
        self.members_DOF.append([])
        self.members_L.append([])
        self.members_PartD_unreleced.append([])
        self.members_PartD_releced.append([])
        self.members_T.append([])
        self.members_K.append([])

        self.members_PointLoads.append([])
        self.members_DistLoads.append([])
        self.members_SelfWeight.append([])

    def defSupport(self, node_index: int, support_DX: bool = False, support_DY: bool = False, support_DZ: bool = False, support_RX: bool = False, support_RY: bool = False, support_RZ: bool = False):
        self.nodes_Support[node_index] = [support_DX, support_DY, support_DZ, support_RX, support_RY, support_RZ]

    def defReleases(self, member_index: int, Dxi: bool = False, Dyi: bool = False, Dzi: bool = False, Rxi: bool = False, Ryi: bool = False, Rzi: bool = False,
            Dxj: bool = False, Dyj: bool = False, Dzj: bool = False, Rxj: bool = False, Ryj: bool = False, Rzj: bool = False):
        self.members_Releases[member_index] = [Dxi, Dyi, Dzi, Rxi, Ryi, Rzi, Dxj, Dyj, Dzj, Rxj, Ryj, Rzj]

    def addNodeLoad(self, node_index: int, Px: float = 0, Py: float = 0, Pz: float = 0, Mx: float = 0, My: float = 0, Mz: float = 0, case: int = 0):
        caseFound = False
        for load_case in self.nodes_loads[node_index]:
            caseINDEX, load = load_case[0], load_case[1]
            if caseINDEX == case:
                load[0] += Px
                load[1] += Py
                load[2] += Pz
                load[3] += Mx
                load[4] += My
                load[5] += Mz
                caseFound = True
        if not caseFound:
            self.nodes_loads[node_index].append([case, [Px, Py, Pz, Mx, My, Mz]])
            if case not in self.casses:
                self.casses.append(case)

    def addMemberPointLoad(self, member_index: int, x: float, Px: float = 0, Py: float = 0, Pz: float = 0, Mx: float = 0, My: float = 0, Mz: float = 0, case: int = 0):
        caseFound = False
        for loads_case in self.members_PointLoads[member_index]:
            caseINDEX, loads = loads_case[0], loads_case[1]
            if caseINDEX == case:
                caseFound = True
                loadFound = False
                loadsX, loadsMag = loads[0], loads[1]
                for i in range(len(loadsX)):
                    if loadsX == x:
                        loadsMag[i][0] += Px
                        loadsMag[i][1] += Py
                        loadsMag[i][2] += Pz
                        loadsMag[i][3] += Mx
                        loadsMag[i][4] += My
                        loadsMag[i][5] += Mz
                        loadFound = True
                if not loadFound:
                    loadsX.append(x)
                    loadsMag.append([Px, Py, Pz, Mx, My, Mz])
        if not caseFound:
            self.members_PointLoads[member_index].append([case, [[x], [[Px, Py, Pz, Mx, My, Mz]]]])
            if case not in self.casses:
                self.casses.append(case)

    def addMemberDistLoad(self, member_index: int, x1: float, x2: float, wx1: float = 0, wx2: float = 0, wy1: float = 0, wy2: float = 0, wz1: float = 0, wz2: float = 0, case: int = 0):
        caseFound = False
        for loads_case in self.members_DistLoads[member_index]:
            caseINDEX, loads = loads_case[0], loads_case[1]
            if caseINDEX == case:
                caseFound = True
                loads.append([[x1, x2], [wx1, wx2, wy1, wy2, wz1, wz2]])
        if not caseFound:
            self.members_DistLoads[member_index].append([case, [[[x1, x2], [wx1, wx2, wy1, wy2, wz1, wz2]]]])
            if case not in self.casses:
                self.casses.append(case)

    # TODO needs to be implemented
    # def addMemberSelfWeight(self, case: int = 0, factorX: float = 0, factorY: float = 0, factorZ: float = -1):
    #     caseADDED = False
    #     for selfWeight in self.members_SelfWeight:
    #         caseFound = False
    #         for load_case in selfWeight:
    #             if load_case[0] == case:
    #                 load_case[1][0] += factorX
    #                 load_case[1][1] += factorY
    #                 load_case[1][2] += factorZ
    #                 caseFound = True
    #         if not caseFound:
    #             caseADDED = True
    #             selfWeight.append([case, [factorX, factorY, factorZ]])
    #
    #     if caseADDED and case not in self.casses:
    #         self.casses.append([case])

    def preAnalysis_linear(self):
        self.materials = np.array(self.materials)
        self.D_unknown, self.D_known = hf.partD(self)
        self.members_DOF, self.members_L, self.members_PartD_unreleced, self.members_PartD_releced, self.members_T = hf.prepMembers(self)

    def analysis_linear(self, getWeight = False, getReactions = False, getInternalForces = False, getStresses = False):
        numN = len(self.nodes_cord)
        numM = len(self.members)
        numC = len(self.casses)

        self.members_CrossSectionProps = np.array(self.members_CrossSectionProps)

        pointLoads, distLoads = hf.assembleLoads(self)

        k_local = hf.get_k_local_ARRAY(self, self.members_L)

        k11, k12, k21, k22 = hf.memberPart_k_ARRAY(k_local, self.members_PartD_unreleced, self.members_PartD_releced, numM)

        fer_unc_ARRAY = hf.get_fer_unc_ARRAY(self, pointLoads, distLoads)

        fer1, fer2 = hf.memberPart_fer_ARRAY(fer_unc_ARRAY, self.members_PartD_unreleced, self.members_PartD_releced)

        ferCondensed = hf.get_fer_ARRAY(self, k12, k22, fer1, fer2, numM, numC)

        FER1, FER2 = hf.getGlobalFixedEndReactionVector_ARRAY(self, self.members_T, self.D_unknown, self.D_known, ferCondensed)

        P1, P2 = hf.partitionedGlobalNodalForceVector(self, numN)

        k_global_members = hf.k_member_make_global(k_local, self.members_T)

        K_global = hf.get_K_Global(self, k_global_members, numN, numM)

        K11, K12, K21, K22 = hf.partition_K_gloabl(K_global, self.D_unknown, self.D_known)

        D, DX, DY, DZ, RX, RY, RZ = hf.get_D(K11, K12, P1, FER1, self.D_unknown, self.D_known, numN, numC)

        weight = None
        reactions = None
        internalForces = None
        stresses = None

        if getWeight:
            weight = hf.getWeight(self)

        if getReactions or getInternalForces:
            D_members = hf.get_member_direction_deflections(self, DX, DY, DZ, RX, RY, RZ, numM, numC)
            d = hf.getd(self.members_T, D_members, numM, numC)
            f = hf.getf(k_local, d, ferCondensed, numM, numC)

            if getReactions:
                F = hf.getF(self.members_T, f, numM, numC)
                reactions = hf.getReactions(self, F, numC, numM, numN)

            if getInternalForces:
                hf.solveInternalForces(self, pointLoads, distLoads, f, fer_unc_ARRAY, d, numM, numC)
                # TODO

        if getStresses:
            #TODO
            pass

        return D, DX, DY, DZ, RX, RY, RZ, weight, reactions, internalForces, stresses

if __name__ == '__main__':
    simple_beam = Frame3D_T()

    simple_beam.addNode(0, 0, 0)
    simple_beam.addNode(168, 0, 0)
    # simple_beam.addNode(168, 5, 0)

    simple_beam.addMaterial(29000, 11200, 0.3, 2.836e-4)

    simple_beam.addMember(0, 1, 0, True, 20, 100, 150, 250)
    # simple_beam.addMember(0, 2, 0, True, 20, 100, 150, 250)

    simple_beam.defSupport(0, True, True, True, False, False, False)
    simple_beam.defSupport(1, True, True, True, True, False, False)
    # simple_beam.defSupport(0, True, True, True, True, True, True)
    # simple_beam.defSupport(1, True, True, True, True, True, True)
#'M1', 'Fy', -0.01, -0.01, 0, 168
    #simple_beam.addNodeLoad(0, Pz=1, case=0)
    simple_beam.addMemberPointLoad(0, 50, Py=10, case=0)
    # simple_beam.addMemberPointLoad(0, 2, Pz=1, case=0)
    # simple_beam.addMemberPointLoad(0, 2, Pz=1, case=1)
    # simple_beam.addMemberDistLoad(0,0,5,5,2,0,0)
    #simple_beam.addMemberDistLoad(0, 0, 168, wy1 = -0.01, wy2 = -0.01, case=0)
    # simple_beam.addMemberSelfWeight()
    # simple_beam.addMemberSelfWeight(case=1)
    simple_beam.preAnalysis_linear()
    print(simple_beam.analysis_linear(getWeight=True,getInternalForces=True))


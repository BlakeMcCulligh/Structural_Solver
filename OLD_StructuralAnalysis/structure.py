import numpy as np

from OLD_StructuralAnalysis.Solvers.frame2DSolver import Frame2D
from OLD_StructuralAnalysis.Solvers.frame3DSolver import Frame3D
from OLD_StructuralAnalysis.Solvers.truss2DSolver import Truss2D
from OLD_StructuralAnalysis.Solvers.truss3DSolver import Truss3D
from OLD_StructuralAnalysis.member import Member
from OLD_StructuralAnalysis.crossSection import CrossSection
from OLD_StructuralAnalysis.node import Node

class Structure:
    def __init__(self):
        """
        Handels everything to do with the structure.
        """

        self.solveObject = None # the object that handels the solving of the structure (---Solver object)

        self.isTruss = False # is the structure a truss. (Boolean)
        self.is3D = False # is the structure 3D. (boolean)

        self.nodes = [] # list of Node objects for the structure (list)
        self.members = [] # list of member objects for the structure (list)
        self.crossSections = [] # list of cross-section objects for the structure (list)

        self.optimizationResults = None # the Results of a member cross-secction optimization (list)

    def setTruss(self):
        self.isTruss = True

    def set3D(self):
        self.is3D = True

    def setFrame(self):
        self.isTruss = False

    def set2D(self):
        self.is3D = False

    def addNode(self, cords: list):
        """
        Adds a Node to the structure.
        :param cords: list of cords [location, y] if 2D, [location, y, z] if 3D
        """
        self.nodes.append(Node(self.isTruss, self.is3D, cords))

    def addMember(self, nodeIndexs: list, crossSectionIndex: int):
        """
        Adds a member to the structure.
        :param nodeIndexs: list of the Node indexs the member is to conect [Node 1, Node 2].
        :param crossSectionIndex: i of the cross-section the member is to be.
        """
        self.members.append(Member(self.nodes, self.crossSections, self.is3D, self.isTruss, nodeIndexs, crossSectionIndex))

    def addCrossSection(self, Optimize: bool, E, G = None, A = None, I = None, I_weak = None, J = None, memberType = None, minBounds = None, maxBounds = None):
        """
        Adds a cross-section to the structure.
        :param Optimize: is the cross-seection an optimization cross-section. (boolean)
        :param E: Elastic modulus of the cross-section. (float)
        :param G: shear modulus of the cross-section. (float)
        :param A: Area of the cross-section. (float)
        :param I: moment of inertia of the cross-section around the strong axis if 3D. (float)
        :param I_weak: moment of inertia of the cross-section around the weak axis. Only needed for 3D. (float)
        :param J: moment of inertia of the cross-section for tortian. Only needed for 3D. (float)
        :param memberType: What type of member is the cross-section. Only needed for optimization (string)
        :param minBounds: Minimum bounds on the optimization properties of the cross-section. Only needed for optimization (list)
        :param maxBounds: Maximum bounds on the optimization properties of the cross-section. Only needed for optimization (list)
        """

        newCrossSection = CrossSection()
        if Optimize:
            if self.is3D:
                if self.isTruss:
                    newCrossSection.optAdd3DTruss(E, minBounds, maxBounds)
                else:
                    newCrossSection.optAdd3DFrame(E, G, memberType, minBounds, maxBounds)
            else:
                if self.isTruss:
                    newCrossSection.optAdd2DTruss(E, minBounds, maxBounds)
                else:
                    newCrossSection.optAdd2DFrame(E, memberType, minBounds, maxBounds)
        else:
            if self.is3D:
                if self.isTruss:
                    newCrossSection.add3DTruss(A, E)
                else:
                    newCrossSection.add3DFrame(A, E, G, I, I_weak, J)
            else:
                if self.isTruss:
                    newCrossSection.add2DTruss(A, E)
                else:
                    newCrossSection.add2DFrame(A, E, I)
        self.crossSections.append(newCrossSection)

    def addSupport(self, nodeIndex: int, support: list):
        """
        Adds a support to the structure.
        :param nodeIndex: i of the Node the support is to be applied to. (int)
        :param support: what degress of fredom the support is to be applied to. (list) of booleans
                        if 2D truss: length = 2
                        if 3D truss: length = 3
                        if 2D frame: length = 3
                        if 3D frame: length = 6
        """
        self.nodes[nodeIndex].addSupport(support)

    def addNodeLoad(self, nodeIndex: int, load: list):
        """
        Adds a point load at a Node to the structure.
        :param nodeIndex: i of the Node the load is to be applied to. (int)
        :param load: magnitude of the load in each degree of fredom. (list) of floats
                     if 2D truss: length = 2
                     if 3D truss: length = 3
                     if 2D frame: length = 3
                     if 3D frame: length = 6
        """
        self.nodes[nodeIndex].addLoad(load)

    def addPointLoad(self, memberIndex: int, load: list, location: float):
        self.members[memberIndex].addPointLoad(load, location)

    def addUniformlyDistributedLoad(self, memberIndex: int, load: list):
        self.members[memberIndex].addUniformlyDistributedLoad(load)

    def addDistributedLoad(self, memberIndex: int, load: list, location: list):
        self.members[memberIndex].addDistributedLoad(load, location)

    def addRelece(self, memberIndex: int, relece: list):
        """
        Adds a relece to the structure. Only used for frames.
        :param memberIndex: i of the member the relece is to be applied to. (int)
        :param relece: what degress of fredom are to be releced for both ends of the member. (list) of floats
                       if 2D frame: length = 6
                       if 3D frame: length = 12
        """
        self.members[memberIndex].addRelece(relece)

    def assembleGeoLists(self):
        """
        Assembles the lists of geometric properties to be sent to the solver.
        helper absFunction for the solve and optimize functions.
        :return: nodes, members, Loads, Supports
        """
        NODES = []
        for node in self.nodes:
            NODES.append(node.cords)

        MEMBERS = []
        for member in self.members:
            MEMBERS.append(member.nodeIndexs)

        LOADS = []
        for i, node in enumerate(self.nodes):
            if any(np.array(node.load) != 0):
                LOADS.append([i] + node.load)

        SUPPORTS = []
        for i, node in enumerate(self.nodes):
            if any(node.support):
                SUPPORTS.append([i] + node.support)

        return NODES, MEMBERS, LOADS, SUPPORTS

    def solve(self):
        """
        Linear solves the structure using stiffness matrices.
        """
        NODES, MEMBERS, LOADS, SUPPORTS = self.assembleGeoLists()

        CROSSSECTIONS = []
        for member in self.members:
            CROSSSECTIONS.append(member.crossSectionIndex)

        if self.isTruss:
            if self.is3D:
                A = []
                E = []
                for crossSection in self.crossSections:
                    A.append(crossSection.A)
                    E.append(crossSection.E)

                self.solveObject = Truss3D()
                self.solveObject.setGeom(NODES, MEMBERS, LOADS, SUPPORTS)
                self.solveObject.setA(A)
                self.solveObject.setE(E)
                self.solveObject.setMemberGroup(CROSSSECTIONS)
                self.solveObject.solveLinear()
            else:
                A = []
                E = []
                for crossSection in self.crossSections:
                    A.append(crossSection.A)
                    E.append(crossSection.E)

                self.solveObject = Truss2D()
                self.solveObject.setGeom(NODES, MEMBERS, LOADS, SUPPORTS)
                self.solveObject.setA(A)
                self.solveObject.setE(E)
                self.solveObject.setMemberGroup(CROSSSECTIONS)
                self.solveObject.solveLinear()

        else:

            RELEASES = []
            for i, member in enumerate(self.members):
                if any(member.relece):
                    RELEASES.append([i] + member.relece)

            if self.is3D:
                A = []
                E = []
                G = []
                J = []
                Iy = []
                Iz = []
                for crossSection in self.crossSections:
                    A.append(crossSection.A)
                    E.append(crossSection.E)
                    G.append(crossSection.G)
                    J.append(crossSection.J)
                    Iy.append(crossSection.I_weak)
                    Iz.append(crossSection.I_main)

                self.solveObject = Frame3D()
                self.solveObject.setGeom(NODES, MEMBERS, LOADS, SUPPORTS, RELEASES)
                self.solveObject.setA(A)
                self.solveObject.setMaterialProperties(E, G)
                self.solveObject.setStiffnesses(J, Iy, Iz)
                self.solveObject.setMemberGroups(CROSSSECTIONS)
                self.solveObject.solveLinear()

            else:
                A = []
                E = []
                I = []
                for crossSection in self.crossSections:
                    A.append(crossSection.A)
                    E.append(crossSection.E)
                    I.append(crossSection.I_main)

                self.solveObject = Frame2D()
                self.solveObject.setGeom(NODES, MEMBERS, LOADS, SUPPORTS, RELEASES)
                self.solveObject.setA(A)
                self.solveObject.setE(E)
                self.solveObject.setI(I)
                self.solveObject.setMemberGroups(CROSSSECTIONS)
                self.solveObject.solveLinear()

        self.solveInternalForces()

    def solveInternalForces(self):
        if self.isTruss:
            if self.is3D:
                self.solveObject.calcMemberDeflections()
                self.solveObject.calcInternalForces()
                #self.solveObject.calcNomalStresses()
            else:
                self.solveObject.calcMemberDeflections()
                self.solveObject.calcInternalForces()
                #self.solveObject.calcNomalStresses()
        else:
            if self.is3D:
                pass
                # TODO
            else:
                pass
                # TODO

    def optimize(self, initalGuess: list):
        """
        Optimizes the structures cross-sections.
        :param initalGuess: inital guess of the cross-section optimization values. (list of floats)
        """
        NODES, MEMBERS, LOADS, SUPPORTS = self.assembleGeoLists()

        CROSSSECTIONS = []
        for member in self.members:
            CROSSSECTIONS.append(member.crossSectionIndex)

        if self.isTruss:
            if self.is3D:
                E = []
                for crossSection in self.crossSections:
                    E.append(crossSection.E)

                self.solveObject = Truss3D()
                self.solveObject.setGeom(NODES, MEMBERS, LOADS, SUPPORTS)
                self.solveObject.setE(E)
                self.solveObject.setMemberGroup(CROSSSECTIONS)
                self.optimizationResults = self.solveObject.optimize(self.crossSections, initalGuess)

            else:
                E = []
                for crossSection in self.crossSections:
                    E.append(crossSection.E)

                self.solveObject = Truss2D()
                self.solveObject.setGeom(NODES, MEMBERS, LOADS, SUPPORTS)
                self.solveObject.setE(E)
                self.solveObject.setMemberGroup(CROSSSECTIONS)
                self.optimizationResults = self.solveObject.optimize(self.crossSections, initalGuess)

        else:
            RELEASES = []
            for i, member in enumerate(self.members):
                if any(member.relece):
                    RELEASES.append([i] + member.relece)

            if self.is3D:
                E = []
                G = []
                for crossSection in self.crossSections:
                    E.append(crossSection.E)
                    G.append(crossSection.G)

                self.solveObject = Frame3D()
                self.solveObject.setGeom(NODES, MEMBERS, LOADS, SUPPORTS, RELEASES)
                self.solveObject.setMaterialProperties(E, G)
                self.solveObject.setMemberGroups(CROSSSECTIONS)
                self.optimizationResults = self.solveObject.optimize(self.crossSections, initalGuess)

            else:
                E = []
                for crossSection in self.crossSections:
                    E.append(crossSection.E)

                self.solveObject = Frame2D()
                self.solveObject.setGeom(NODES, MEMBERS, LOADS, SUPPORTS, RELEASES)
                self.solveObject.setE(E)
                self.solveObject.setMemberGroups(CROSSSECTIONS)
                self.optimizationResults = self.solveObject.optimize(self.crossSections, initalGuess)

    def getDeflections(self):
        if self.solveObject is not None:
            #todo
            return self.solveObject.U
        else:
            print("The structure has not been solved yet.")
            return None

    def getForces(self):
        if self.solveObject is not None:
            # todo
            return self.solveObject.F_m_internal
        else:
            print("The structure has not been solved yet.")
            return None

    # def getStresses(self):
    #     if self.solveObject is not None:
    #         #todo
    #         return self.solveObject.sigma
    #     else:
    #         print("The structure has not been solved yet.")
    #         return None

    def printOptimizationResults(self):
        """
        prints the Results for the optimization variables returned from the optimizer.
        """
        if self.optimizationResults is not None:
            if self.isTruss:
                if self.is3D:
                    print("A: ", self.optimizationResults)
                else:
                    print("A: ", self.optimizationResults)
            else:
                if self.is3D:
                    print("Results: ", self.optimizationResults)
                else:
                    print("Results: ", self.optimizationResults)
        else:
            print("The structure has not been optimized yet.")






def testTruss2D():
    S = Structure()
    S.setTruss()

    S.addNode([0, 0])
    S.addNode([1, 0])
    S.addNode([1, 1])
    S.addNode([0, 1])

    S.addCrossSection(False, E=1, A = 2)

    S.addMember([0, 1], 0)
    S.addMember([1, 2], 0)
    S.addMember([2, 3], 0)
    S.addMember([3, 0], 0)
    S.addMember([0, 2], 0)

    S.addNodeLoad(2, [5, 0])
    S.addNodeLoad(3, [0, 2])

    S.addPointLoad(4, [5, 0], 0.5)
    S.addUniformlyDistributedLoad(2, [2,0])
    S.addDistributedLoad(3,[0,0,0,1], [0, 1])



    S.addSupport(0, [True, True])
    S.addSupport(1, [False, True])

    S.solve()
    print(S.getDeflections())
    print(S.getForces())
    #print(S.getStresses())

testTruss2D()

#
# S = Structure()
# #S.set3D()
# S.setTruss()
#
# S.AddPrintNode([0,0,0])
# S.AddPrintNode([1,0,0])
# S.AddPrintNode([1,1,0])
# S.AddPrintNode([0,1,0])
#
# S.addCrossSection(True, E=1, G=1, memberType="SquareHSS", minBounds=[0.1, 0.01], maxBounds=[10, 0.09])
#
# S.addMember([0,1], 0)
# S.addMember([1,2], 0)
# S.addMember([2,3], 0)
# S.addMember([3,0], 0)
# S.addMember([0,2], 0)
#
# S.addNodeLoad(2, [5, 0, 0, 0, 0, 0])
# S.addNodeLoad(3, [0, 2, 0, 0, 0, 0])
# #Supports = [[0,True,True,True,True,False,False],[1,False,True,False,True,False,False]]
# S.addSupport(0, [True,True, True,True, False, False])
# S.addSupport(1, [False,True, False, True, False, False])
#
# S.addRelece(1, [0,0,0,0,0,1,0,0,0,0,0,1])
# S.addRelece(3, [0,0,0,0,0,1,0,0,0,0,0,1])
# S.addRelece(4, [0,0,0,0,0,1,0,0,0,0,0,1])
# # S.addRelece(1, [0,0,1,0,0,1])
# # S.addRelece(3, [0,0,1,0,0,1])
# # S.addRelece(4, [0,0,1,0,0,1])
#
# # S.solve()
# # S.printDeflections()
#
# S.optimize([1,1])
# S.printOptimizationResults()




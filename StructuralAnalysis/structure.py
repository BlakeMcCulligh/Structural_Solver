import numpy as np

from StructuralAnalysis.LinearSolver.frame2DSolver import Frame2D
from StructuralAnalysis.LinearSolver.frame3DSolver import Frame3D
from StructuralAnalysis.LinearSolver.truss2DSolver import Truss2D
from StructuralAnalysis.LinearSolver.truss3DSolver import Truss3D
from StructuralAnalysis.member import Member
from StructuralAnalysis.crossSection import CrossSection
from StructuralAnalysis.node import Node

class Structure:
    def __init__(self):
        self.solveObject = None
        self.isTruss = False
        self.is3D = False

        self.nodes = []
        self.members = []
        self.crossSections = []

        self.optimizationResults = None

    def setTruss(self):
        self.isTruss = True

    def set3D(self):
        self.is3D = True

    def addNode(self, cords: list):
        self.nodes.append(Node(self.isTruss, self.is3D, cords))

    def addMember(self, nodeIndexs: list, memberGroupIndex: int):
        self.members.append(Member(self.nodes, self.crossSections, self.is3D, self.isTruss, nodeIndexs, memberGroupIndex))

    def addCrossSection(self, Optimize, E, A = None, I = None, I_weak = None, G = None, J = None, memberType = None, minBounds = None, maxBounds = None):
        newCrossSection = CrossSection()
        if Optimize:
            if self.is3D:
                if self.isTruss:
                    newCrossSection.optAdd3DTruss(E)
                else:
                    newCrossSection.optAdd3DFrame(E, G, memberType, minBounds, maxBounds)
            else:
                if self.isTruss:
                    newCrossSection.optAdd2DTruss(E)
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
        self.nodes[nodeIndex].addSupport(support)

    def addNodeLoad(self, nodeIndex: int, load: list):
        self.nodes[nodeIndex].addLoad(load)

    def addRelece(self, memberIndex: int, relece: list):
        self.members[memberIndex].addRelece(relece)

    def assembleGeoLists(self):
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

    def optimize(self, initalGuess: list, MinArea = None, MaxArea = None):
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
                self.optimizationResults = self.solveObject.optimize(MinArea, MaxArea, initalGuess)

            else:
                E = []
                for crossSection in self.crossSections:
                    E.append(crossSection.E)

                self.solveObject = Truss2D()
                self.solveObject.setGeom(NODES, MEMBERS, LOADS, SUPPORTS)
                self.solveObject.setE(E)
                self.solveObject.setMemberGroup(CROSSSECTIONS)
                self.optimizationResults = self.solveObject.optimize(MinArea, MaxArea, initalGuess)

        else:

            RELEASES = []
            for i, member in enumerate(self.members):
                if any(member.relece):
                    RELEASES.append([i] + member.relece)

            if self.is3D:
                pass # todo
            else:
                E = []
                for crossSection in self.crossSections:
                    E.append(crossSection.E)

                self.solveObject = Frame2D()
                self.solveObject.setGeom(NODES, MEMBERS, LOADS, SUPPORTS, RELEASES)
                self.solveObject.setE(E)
                self.solveObject.setMemberGroups(CROSSSECTIONS)
                self.optimizationResults = self.solveObject.optimize(self.crossSections, initalGuess)



    def printDeflections(self):
        if self.solveObject is not None:
            print("Defflections: ", self.solveObject.U)
        else:
            print("The structure has not been solved yet.")

    def printOptimizationResults(self):
        if self.optimizationResults is not None:
            if self.isTruss:
                if self.is3D:
                    print("A: ", self.optimizationResults)
                else:
                    print("A: ", self.optimizationResults)
            else:
                if self.is3D:
                    pass  # todo
                else:
                   # todo
                    print("Results: ", self.optimizationResults)
        else:
            print("The structure has not been optimized yet.")


S = Structure()
#S.set3D()
#S.setTruss()

S.addNode([0,0])
S.addNode([1,0])
S.addNode([1,1])
S.addNode([0,1])

S.addCrossSection(True, E=1, memberType="SquareHSS", minBounds=[0.1, 0.01], maxBounds=[10, 0.09])

S.addMember([0,1], 0)
S.addMember([1,2], 0)
S.addMember([2,3], 0)
S.addMember([3,0], 0)
S.addMember([0,2], 0)

S.addNodeLoad(2, [5, 0, 0])
S.addNodeLoad(3, [0, 2, 0])

S.addSupport(0, [True,True, False])
S.addSupport(1, [False,True, False])

# S.addRelece(1, [0,0,0,0,0,1,0,0,0,0,0,1])
# S.addRelece(3, [0,0,0,0,0,1,0,0,0,0,0,1])
# S.addRelece(4, [0,0,0,0,0,1,0,0,0,0,0,1])
S.addRelece(1, [0,0,1,0,0,1])
S.addRelece(3, [0,0,1,0,0,1])
S.addRelece(4, [0,0,1,0,0,1])

# S.solve()
#
# S.printDeflections()

S.optimize([1,1])
S.printOptimizationResults()




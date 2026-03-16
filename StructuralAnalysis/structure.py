import numpy as np

from StructuralAnalysis.LinearSolver.frame2DSolver import Frame2D
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



    def setTruss(self):
        self.isTruss = True

    def set3D(self):
        self.is3D = True

    def addNode(self, cords: list):
        self.nodes.append(Node(self.isTruss, self.is3D, cords))

    def addMember(self, nodeIndexs: list, memberGroupIndex: int):
        self.members.append(Member(self.nodes, self.crossSections, self.is3D, self.isTruss, nodeIndexs, memberGroupIndex))

    def addCrossSection(self, Optimize, E, A = None, I = None, I_weak = None):
        newCrossSection = CrossSection()

        if Optimize:
            if self.is3D:
                if self.isTruss:
                    newCrossSection.optAdd3DTruss(E)
                else:
                    newCrossSection.optAdd3DFrame(E)
            else:
                if self.isTruss:
                    newCrossSection.optAdd2DTruss(E)
                else:
                    newCrossSection.optAdd2DFrame(E)

        else:
            if self.is3D:
                if self.isTruss:
                    newCrossSection.add3DTruss(A, E)
                else:
                    newCrossSection.add3DFrame(A, E, I, I_weak)
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

    def solve(self):
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
                print("Yet To Be Implemented")
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

    def printDeflections(self):
        if self.solveObject is not None:
            print("Defflections: ", self.solveObject.U)
        else:
            print("The structure has not been solved yet.")


# Nodes = [[0,0],[1,0],[1,1],[0,1]]
# Members = [[0,1],[1,2],[2,3],[3,0],[0,2]]
# Loads = [[2,5,0],[3,0,2]]
# Supports = [[0,True,True],[1,False,True]]
# MemberGroups_set = [0,0,0,0,0]

# Nodes = [[0,0],[1,0],[1,1],[0,1]]
# Members = [[0,1],[1,2],[2,3],[3,0],[0,2]]
# Loads = [[2,5,0,0],[3,0,2,0]]
# Supports = [[0,True,True,False],[1,False,True,False]]
# Releases = [[0,0,0,0,0,0,0],[1,0,0,1,0,0,1],[2,0,0,0,0,0,0],[3,0,0,1,0,0,1],[4,0,0,1,0,0,1]]
# MemberGroups_set = [0,0,0,0,0]
# E_set = [1]
# A_set = [1]
# I_set = [1]

S = Structure()
S.setTruss()
S.set3D()

S.addNode([0,0,0])
S.addNode([1,0,0])
S.addNode([1,1,0])
S.addNode([0,1,0])

S.addCrossSection(False, 1, 1)

S.addMember([0,1], 0)
S.addMember([1,2], 0)
S.addMember([2,3], 0)
S.addMember([3,0], 0)
S.addMember([0,2], 0)

# S.addNodeLoad(2, [5,0])
# S.addNodeLoad(3, [0,2])
S.addNodeLoad(2, [5, 0, 0])
S.addNodeLoad(3, [0, 2, 0])

# S.addSupport(0, [True,True])
# S.addSupport(1, [False,True])
# S.addSupport(0, [True,True, False])
# S.addSupport(1, [False,True,False])
S.addSupport(0, [True,True, True])
S.addSupport(1, [False,True,True])
S.addSupport(2, [False,False, True])
S.addSupport(3, [False,False,True])

# S.addRelece(1, [False, False, True, False, False, True])
# S.addRelece(3, [False, False, True, False, False, True])
# S.addRelece(4, [False, False, True, False, False, True])

S.solve()

S.printDeflections()




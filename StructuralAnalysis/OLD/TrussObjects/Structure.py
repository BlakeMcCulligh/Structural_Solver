import numpy as np

from StructuralAnalysis.OLD.GeneralObjects.CrossSection import CrossSection
from StructuralAnalysis.OLD.GeneralObjects.Material import Material
from StructuralAnalysis.OLD.TrussObjects.LoadCase import LoadCase
from StructuralAnalysis.OLD.TrussObjects.Member import Member
from StructuralAnalysis.OLD.TrussObjects.Node import Node


class Truss:
    def __init__(self):

        self.memberGroupings = []

        self.members = []
        self.nodes = []
        self.suports = [] # [index of node that is supported, direction of the support (x = 1, y = 2, z = 3)]

        self.StiffnessMatrix = None

        self.nodesFree = None
        self.nodesSupported = None

        self.weight = None

        self.loadCases = []

        self.unitLoadNodeIndex = None

    def addNode(self, cords):

        index = None
        nodeAlreadyExists = False
        if len(cords) == 3:
            for i in range(len(self.nodes)):
                if self.nodes[i].cords[0] == cords[0] and self.nodes[i].cords[1] == cords[1] and self.nodes[i].cords[2] == cords[2]:
                    index = i
                    nodeAlreadyExists = True
        else:
            for i in range(len(self.nodes)):
                if self.nodes[i].cords[0] == cords[0] and self.nodes[i].cords[1] == cords[1] and self.nodes[i].cords[2] == cords[2]:
                    index = i
                    nodeAlreadyExists = True

        if not nodeAlreadyExists:
            index = len(self.nodes)
            newNode = Node(cords, index)
            self.nodes.append(newNode)

        return index

    def addMember(self, nodesIndex, crossSection, material):
        newMember = Member([self.nodes[nodesIndex[0]], self.nodes[nodesIndex[1]]], crossSection, material)
        self.members.append(newMember)

    def addSupport(self, nodeIndex, direction):
        """
        :param nodeIndex: The index of the node being supported
        :param direction: x = 1, y = 2, z = 3
        """

        self.suports.append([nodeIndex, direction])

    # def addNodePointLoad(self, nodeIndex, direction, magnitued):
    #     """
    #
    #     :param nodeIndex: The index of the node the load is applyed to
    #     :param direction: (x = 1, y = 2, z = 3)
    #     :param magnitued: the magnitued of the load being applyed. A negative value indecated in the negative direcction of the direction indicated.
    #     :return:
    #     """
    #
    #     self.nodePointLoads.append([nodeIndex, direction, magnitued])

    def AssembleGlobalStiffnessMatrix(self):
        """
        Creates the global stiffness matrix
        """

        self.arangeNodeFixarity()
        for loadCase in self.loadCases:
            loadCase.arangeForces(self.nodes[0].numDimensions, len(self.nodes))

        if self.nodes[0].numDimensions == 3:
            K = np.zeros((3 * len(self.nodes),3 * len(self.nodes)))
            for i in range(len(self.members)):

                eldofs = np.concatenate([np.arange(3 * self.members[i].nodes[0].index, 3 * (self.members[i].nodes[0].index + 1)),
                                         np.arange(3 * self.members[i].nodes[1].index, 3 * (self.members[i].nodes[1].index + 1))])

                K[np.ix_(eldofs, eldofs)] += self.members[i].GlobalStiffnessMatrix

            self.StiffnessMatrix = K

        else:
            K = np.zeros((2 * len(self.nodes),2 * len(self.nodes)))
            for i in range(len(self.members)):

                aux = [2 * self.members[i].nodes[0].index, 2 * self.members[i].nodes[1].index]
                index = np.r_[aux[0]:aux[0] + 2, aux[1]:aux[1] + 2]
                K[np.ix_(index, index)] += self.members[i].GlobalStiffnessMatrix

            self.StiffnessMatrix = K

    def solveStiffnessMatix(self):

        bars = []
        L = []
        A = []
        E = []
        for m in self.members:
            bars.append([m.nodes[0].index, m.nodes[1].index])
            L.append(m.length)
            A.append(m.crossSection.A)
            E.append(m.material.E)
        bars = np.array(bars)
        L = np.array(L)
        A = np.array(A)
        E = np.array(E)

        nodes = []
        for n in self.nodes:
            nodes.append(n.cords)
        nodes = np.array(nodes)

        d = nodes[bars[:, 1], :] - nodes[bars[:, 0], :]
        angle = d.T / L
        a = np.concatenate((-angle.T, angle.T), axis=1)

        for j in range(len(self.members)):
            self.members[j].normalForces = [0] * len(self.loadCases)

        for j in range(len(self.nodes)):
            self.nodes[j].deformation = [0] * len(self.loadCases)

        for i in range(len(self.loadCases)):

            self.loadCases[i].deformation[self.nodesFree] = np.linalg.solve(self.StiffnessMatrix[np.ix_(self.nodesFree, self.nodesFree)],
                                           self.loadCases[i].applyedForce[self.nodesFree])
            print(self.loadCases[i].deformation)
            U =  self.loadCases[i].deformation.reshape(len(self.nodes), self.nodes[0].numDimensions)
            print(U)
            for j in range(len(self.nodes)):
                self.nodes[j].deformation[i] = U[j]
                print(U[j])

            u = np.concatenate((U[bars[:, 0]], U[bars[:, 1]]), axis=1)
            N = E[:] * A[:] / L[:] * (a[:] * u[:]).sum(axis=1)

            self.loadCases[i].internalForces = np.array(N)[np.newaxis].T

            for j in range(len(self.members)):
                self.members[j].normalForces[i] = self.loadCases[i].internalForces[j]

            self.loadCases[i].reactions = np.zeros((len(self.nodes)* self.nodes[0].numDimensions,1))
            self.loadCases[i].reactions[self.nodesSupported] = self.StiffnessMatrix[self.nodesSupported, :].dot(self.loadCases[i].deformation)

    def arangeNodeFixarity(self):
        nDim = self.nodes[0].numDimensions

        for loadCase in self.loadCases:
            loadCase.deformation = np.zeros(shape=(nDim * len(self.nodes), 1))

        alldofs = []
        for i in range(nDim * len(self.nodes)):
            alldofs.append(i)

        self.nodesSupported = []
        for i in range(len(self.suports)):
            thisSupport = nDim * (self.suports[i][0]) + self.suports[i][1]-1

            self.nodesSupported.append(thisSupport)
            for loadCase in self.loadCases:
                loadCase.deformation[thisSupport] = 0
        self.nodesFree = alldofs.copy()
        self.nodesSupported.sort(reverse=True)
        for index in self.nodesSupported:
            del self.nodesFree[index]

    def calcUnitLoadDisplacements(self):

        ULCasses = []
        for i in range(self.nodes[0].numDimensions):
            ULCasses.append(LoadCase())
            ULCasses[i].addNodePointLoad(self.unitLoadNodeIndex, i + 1, 1)

        for j in range(len(self.nodes)):
            self.nodes[j].unitDeformation = [0] * len(ULCasses)

        for i in range(len(ULCasses)):

            ULCasses[i].deformation[self.nodesFree] = np.linalg.solve(
                self.StiffnessMatrix[np.ix_(self.nodesFree, self.nodesFree)],
                ULCasses[i].applyedForce[self.nodesFree])

            U = ULCasses[i].deformation.reshape(len(self.nodes), self.nodes[0].numDimensions)

            for j in range(len(self.nodes)):
                self.nodes[j].unitDeformation[i] = U[j]

    def calcStresses(self):
        for member in self.members:
            member.calcNormalStress()

    def calcWeigth(self):
        W = 0
        for m in self.members:
            W += m.material.desity * m.crossSection.A * m.length

        self.weight = W

    def checkActiveDisplacementConstraints(self):
        for n in self.nodes:
            n.checkActiveDisplacementConstraint()

truss = Truss()
typMaterial = Material()
typMaterial.E = 1e4
typCorssSection = CrossSection()
typCorssSection.A = 0.111

truss.addNode([0, 120])
truss.addNode([120, 120])
truss.addNode([240, 120])
truss.addNode([360, 120])
truss.addNode([0, 0])
truss.addNode([120, 0])
truss.addNode([240, 0])
truss.addNode([360, 0])

truss.addMember([0,1], typCorssSection, typMaterial)
truss.addMember([1,2], typCorssSection, typMaterial)
truss.addMember([2,3], typCorssSection, typMaterial)
truss.addMember([4,5], typCorssSection, typMaterial)
truss.addMember([5,6], typCorssSection, typMaterial)
truss.addMember([6,7], typCorssSection, typMaterial)
truss.addMember([5,1], typCorssSection, typMaterial)
truss.addMember([6,2], typCorssSection, typMaterial)
truss.addMember([7,3], typCorssSection, typMaterial)
truss.addMember([0,5], typCorssSection, typMaterial)
truss.addMember([4,1], typCorssSection, typMaterial)
truss.addMember([1,6], typCorssSection, typMaterial)
truss.addMember([5,2], typCorssSection, typMaterial)
truss.addMember([2,7], typCorssSection, typMaterial)
truss.addMember([7,3], typCorssSection, typMaterial)

truss.addSupport(0,1)
truss.addSupport(0,2)
truss.addSupport(4,1)
truss.addSupport(4,2)

LC1 = LoadCase()
LC1.addNodePointLoad(7,2,-10)
truss.loadCases.append(LC1)

truss.AssembleGlobalStiffnessMatrix()
truss.solveStiffnessMatix()

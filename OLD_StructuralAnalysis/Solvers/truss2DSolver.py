import numpy as np

from OLD_StructuralAnalysis.CrossSectionOptimization import truss2DOptimizer


class Truss2D:
    def __init__(self):
        self.nodes = None # [location,y]
        self.members = None # [Node 1, Node 2]
        self.loads = None # [Node, location, y]
        self.supports = None # [Node, location, y]
        self.memberGroup = None
        self.A = None # [A] * length members
        self.E = None # [E] * length members

        self.L = None
        self.K_m_local = None
        self.K_m_global = None
        self.T = None
        self.T_T = None
        self.K = None
        self.freeDOF = None
        self.fixedDOF = None
        self.F = None
        self.U = None
        self.U_m_global = None
        self.U_m_local = None
        self.Reactions = None
        self.F_m_internal = None
        self.sigma = None

    def setGeom(self, nodes, members, loads, supports):
        self.nodes = np.array(nodes)
        self.members = np.array(members)
        self.loads = np.array(loads)
        self.supports = np.array(supports)

    def setA(self, A):
        self.A = np.array(A)

    def setE(self, E):
        self.E = np.array(E)

    def setMemberGroup(self, memberGroup):
        self.memberGroup = memberGroup

    def calcLengths(self):
        dx = self.nodes[self.members[:,1],0] - self.nodes[self.members[:,0],0]
        dy = self.nodes[self.members[:,1],1] - self.nodes[self.members[:,0],1]
        self.L = np.sqrt(dx**2 + dy**2)

    def calcLocalStiffnessMatrices(self):
        A = self.A[self.memberGroup]
        E = self.E[self.memberGroup]
        Mult = A * E / self.L
        MultMatrix = np.array([[1,0,-1,0],
         [0,0,0,0],
         [-1,0,1,0],
         [0,0,0,0]])
        self.K_m_local = Mult[:, np.newaxis, np.newaxis] * MultMatrix

    def calcTransformationMatrices(self):
        C = (self.nodes[self.members[:, 1], 0] - self.nodes[self.members[:, 0], 0]) /self.L
        S = (self.nodes[self.members[:, 1], 1] - self.nodes[self.members[:, 0], 1]) /self.L
        T = []
        for c, s in zip(C, S):
            T.append(np.array([
                [c,s,0,0],
                [-s,c,0,0],
                [0,0,c,s],
                [0,0,-s,c]
            ]))
        self.T = np.array(T)
        self.T_T = np.array([t.T for t in self.T])

    def calcTransfromLocalToGlobalStiffnessMatrix(self):
        self.K_m_global = self.T_T @ self.K_m_local @ self.T

    def calcGlobalStiffnessMatrix(self):
        self.calcLocalStiffnessMatrices()
        self.calcTransformationMatrices()
        self.calcTransfromLocalToGlobalStiffnessMatrix()

    def assembleGlobalStiffnessMatrix(self):
        K_global = np.zeros((len(self.nodes) * 2, len(self.nodes) * 2)).tolist()
        for index, m in enumerate(self.members):
            global_dofs = [2 * m[0], 2 * m[0] + 1, 2 * m[1], 2 * m[1] + 1]
            for i in range(4):
                for j in range(4):
                    K_global[int(global_dofs[i])][int(global_dofs[j])] += self.K_m_global[index][i, j]
        self.K = np.array(K_global)

    def calcDegressOfFreedom(self):
        self.fixedDOF = []
        for support in self.supports:
            if support[1]:
                self.fixedDOF.append(support[0] * 2)
            if support[2]:
                self.fixedDOF.append(support[0] * 2 + 1)
        self.freeDOF = list(set(range(len(self.nodes) * 2)) - set(self.fixedDOF))

    def calcForces(self):
        F = np.zeros(len(self.nodes) * 2)
        for load in self.loads:
            F[load[0] * 2] = load[1]
            F[load[0] * 2 + 1] = load[2]
        self.F = F

    def calcDeflections(self):
        Kff = self.K[np.ix_(self.freeDOF, self.freeDOF)]
        Uf = np.linalg.solve(Kff, self.F[self.freeDOF])
        self.U = np.zeros(len(self.nodes) * 2)
        self.U[self.freeDOF] = Uf

    def calcReactions(self):
        self.Reactions = self.F
        self.Reactions[self.fixedDOF] = self.K[self.fixedDOF, :] @ self.U

    def solveLinear(self):
        self.calcLengths()
        self.calcGlobalStiffnessMatrix()
        self.assembleGlobalStiffnessMatrix()
        self.calcDegressOfFreedom()
        self.calcForces()
        self.calcDeflections()
        self.calcReactions()

    def calcMemberDeflections(self):
        self.U_m_global = []
        for i, m in enumerate(self.members):
            self.U_m_global.append([self.U[m[0] * 2], self.U[m[0] * 2 + 1], self.U[m[1] * 2],
                                    self.U[m[1] * 2 + 1]])
        self.U_m_local = np.array([self.T[i] @ self.U_m_global[i] for i in range(len(self.members))])

    def calcInternalForces(self):
        self.F_m_internal = np.array([k @ u for k, u in zip(self.K_m_local, self.U_m_local)])
        print(self.F_m_internal)
        self.F_m_internal = self.F_m_internal[:, 0]

    # def calcNomalStresses(self):
    #     self.sigma = self.F_m_internal / self.A

    def optimizeSolve(self, A):
        self.A = A
        self.calcLocalStiffnessMatrices()
        self.calcTransfromLocalToGlobalStiffnessMatrix()
        self.assembleGlobalStiffnessMatrix()
        self.calcDeflections()

    def optimize(self, crossSections, initalGuess: list):

        minBounds = []
        maxBounds = []

        for crossSection in crossSections:
            minBounds = minBounds + crossSection.minBounds
            maxBounds = maxBounds + crossSection.maxBounds

        self.calcLengths()
        self.calcTransformationMatrices()
        self.calcDegressOfFreedom()
        self.calcForces()

        initalGuess = np.array(initalGuess)

        areas = truss2DOptimizer.optimize(self, [minBounds, maxBounds], initalGuess)

        return areas



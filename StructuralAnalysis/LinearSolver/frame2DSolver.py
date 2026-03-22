import numpy as np

from CrossSectionOptimization import frame2DOptimizer
from StructuralAnalysis.CrossSectionCalculaters import SquareHSS, RectHSS, TubeHSS

class Frame2D:
    def __init__(self):
        self.nodes = None # [x,y]
        self.members = None # [node 1, node 2]
        self.loads = None # [node, x, y, m]
        self.supports = None # [node, x, y, m]
        self.releases = None # [node, x_i, y_i, m_i, x_j, y_j, m_j]
        self.memberGroups = None
        self.A = None # [A] * length members
        self.E = None # [E] * length members
        self.I = None # [I] * length members

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
        self.F_m_internal = None
        self.sigma = None

        self.AFormulas = None
        self.IFormulas = None

        self.numVariables = None

    def setGeom(self, nodes, members, loads, supports, releases):
        self.nodes = np.array(nodes)
        self.members = np.array(members)
        self.loads = np.array(loads)
        self.supports = np.array(supports)
        self.releases = np.array(releases)

    def setMemberGroups(self, MemberGroups):
        self.memberGroups = MemberGroups

    def setA(self, A):
        self.A = np.array(A)

    def setE(self, E):
        self.E = np.array(E)

    def setI(self, I):
        self.I = np.array(I)

    def calcLengths(self):
        dx = self.nodes[self.members[:,1],0] - self.nodes[self.members[:,0],0]
        dy = self.nodes[self.members[:,1],1] - self.nodes[self.members[:,0],1]
        self.L = np.sqrt(dx**2 + dy**2)

    def calcLocalStiffnessMatrices(self):
        self.K_m_local = []
        for i in range(len(self.members)):
            A = self.A[self.memberGroups[i]]
            E = self.E[self.memberGroups[i]]
            I = self.I[self.memberGroups[i]]
            L = self.L[i]
            self.K_m_local.append(np.array([
                [A * E / L, 0, 0, -A * E / L, 0, 0],
                [0, 12 * E * I / L ** 3, 6 * E * I / L ** 2, 0, -12 * E * I / L ** 3, 6 * E * I / L ** 2],
                [0, 6 * E * I / L ** 2, 4 * E * I / L, 0, -6 * E * I / L ** 2, 2 * E * I / L],
                [-A * E / L, 0, 0, A * E / L, 0, 0],
                [0, -12 * E * I / L ** 3, -6 * E * I / L ** 2, 0, 12 * E * I / L ** 3, -6 * E * I / L ** 2],
                [0, 6 * E * I / L ** 2, 2 * E * I / L, 0, -6 * E * I / L ** 2, 4 * E * I / L]
            ]))

    def applyRelecesTOLocalStiffnessMatrices(self):
        for MemberReleases in self.releases:
            index = MemberReleases[0]
            k = self.K_m_local[index]
            MemberReleases = MemberReleases[1:]
            R1_indices = []
            R2_indices = []
            for i in range(6):
                if not MemberReleases[i]:
                    R1_indices.append(i)
                else:
                    R2_indices.append(i)
            k11 = k[R1_indices, :][:, R1_indices]
            k12 = k[R1_indices, :][:, R2_indices]
            k21 = k[R2_indices, :][:, R1_indices]
            k22 = k[R2_indices, :][:, R2_indices]
            k = np.subtract(k11, np.matmul(np.matmul(k12, np.linalg.inv(k22)), k21))
            k[abs(k) < 1e-5] = 0
            i = 0
            for DOF in MemberReleases:
                if DOF:
                    k = np.insert(k, i, 0, axis=0)
                    k = np.insert(k, i, 0, axis=1)
                i += 1
            self.K_m_local[index] = k

        self.K_m_local = np.array(self.K_m_local)

    def calcTransformationMatrices(self):
        C = (self.nodes[self.members[:, 1], 0] - self.nodes[self.members[:, 0], 0]) /self.L
        S = (self.nodes[self.members[:, 1], 1] - self.nodes[self.members[:, 0], 1]) /self.L
        T = []
        for c, s in zip(C, S):
            T.append(np.array([
                [ c, s, 0,  0, 0, 0],
                [-s, c, 0,  0, 0, 0],
                [ 0, 0, 1,  0, 0, 0],
                [ 0, 0, 0,  c, s, 0],
                [ 0, 0, 0, -s, c, 0],
                [ 0, 0, 0,  0, 0, 1]
                ]))
        self.T = np.array(T)
        self.T_T = np.array([t.T for t in self.T])

    def calcTransfromLocalToGlobalStiffnessMatrix(self):
        self.K_m_global = self.T_T @ self.K_m_local @ self.T

    def calcGlobalStiffnessMatrix(self):
        self.calcLocalStiffnessMatrices()
        self.applyRelecesTOLocalStiffnessMatrices()
        self.calcTransformationMatrices()
        self.calcTransfromLocalToGlobalStiffnessMatrix()

    def assembleGlobalStiffnessMatrix(self):
        K_global = np.zeros((len(self.nodes) * 3, len(self.nodes) * 3)).tolist()
        for index, m in enumerate(self.members):
            global_dofs = [3 * m[0], 3 * m[0] + 1, 3 * m[0] + 2, 3 * m[1], 3 * m[1] + 1, 3 * m[1] + 2]
            for i in range(6):
                for j in range(6):
                    K_global[int(global_dofs[i])][int(global_dofs[j])] += self.K_m_global[index][i, j]
        self.K = np.array(K_global)

    def calcDegressOfFreedom(self):
        self.fixedDOF = []
        for support in self.supports:
            if support[1]:
                self.fixedDOF.append(support[0] * 3)
            if support[2]:
                self.fixedDOF.append(support[0] * 3 + 1)
            if support[3]:
                self.fixedDOF.append(support[0] * 3 + 2)
        self.freeDOF = list(set(range(len(self.nodes) * 3)) - set(self.fixedDOF))

    def calcForces(self):
        F = np.zeros(len(self.nodes) * 3)
        for load in self.loads:
            F[load[0] * 3] = load[1]
            F[load[0] * 3 + 1] = load[2]
            F[load[0] * 3 + 2] = load[3]
        self.F = F

    def calcDeflections(self):
        Kff = self.K[np.ix_(self.freeDOF, self.freeDOF)]
        Uf = np.linalg.solve(Kff, self.F[self.freeDOF])
        self.U = np.zeros(len(self.nodes) * 3)
        self.U[self.freeDOF] = Uf

    def calcReactions(self):
        self.F[self.fixedDOF] = self.K[self.fixedDOF, :] @ self.U

    def solveLinear(self):
        self.calcLengths()
        self.calcGlobalStiffnessMatrix()
        self.assembleGlobalStiffnessMatrix()
        self.calcDegressOfFreedom()
        self.calcForces()
        self.calcDeflections()
        self.calcReactions()

    def optimizeSolve(self, X):
        self.A = []
        self.I = []

        j = 0
        for i in range(max(self.memberGroups)+1):
            if self.numVariables[i] == 2:
                self.A.append(self.AFormulas[i](X[j], X[j+1]))
                self.I.append(self.IFormulas[i](X[j], X[j+1]))
                j += 2
            elif self.numVariables[i] == 3:
                self.A.append(self.AFormulas[i](X[j], X[j+1], X[j+2]))
                self.I.append(self.IFormulas[i](X[j], X[j+1], X[j+2]))
                j += 3
            else:
                print("Error: incorect number of variables")

        self.A = np.array(self.A)
        self.I = np.array(self.I)

        self.calcLocalStiffnessMatrices()
        self.calcTransfromLocalToGlobalStiffnessMatrix()
        self.assembleGlobalStiffnessMatrix()
        self.calcDeflections()

    def optimize(self, crossSections: list, initalGuess: list):

        self.AFormulas = []
        self.IFormulas = []
        self.numVariables = []

        minBounds = []
        maxBounds = []

        for crossSection in crossSections:
            minBounds = minBounds + crossSection.minBounds
            maxBounds = maxBounds + crossSection.maxBounds

            if crossSection.Type == "SquareHSS":
                self.AFormulas.append(SquareHSS.getA)
                self.IFormulas.append(SquareHSS.getI)
                self.numVariables.append(2)
            elif crossSection.Type == "RectHSS":
                self.AFormulas.append(RectHSS.getA)
                self.IFormulas.append(RectHSS.getIx)
                self.numVariables.append(3)
            elif crossSection.Type == "TubeHSS":
                self.AFormulas.append(TubeHSS.getA)
                self.IFormulas.append(TubeHSS.getI)
                self.numVariables.append(2)
            elif crossSection.Type == "Angle":
                self.AFormulas.append(RectHSS.getA)
                self.IFormulas.append(RectHSS.getIx)
                self.numVariables.append(3)

        self.calcLengths()
        self.calcTransformationMatrices()
        self.calcDegressOfFreedom()
        self.calcForces()

        initalGuess = np.array(initalGuess)

        areas = frame2DOptimizer.optimize(self, [minBounds, maxBounds], initalGuess)

        return areas

    # def calcMemberDeflections(self):
    #     self.U_m_global = []
    #     for i, m in enumerate(self.members):
    #         self.U_m_global.append([self.U[m[0] * 2], self.U[m[0] * 2 + 1], self.U[m[0] * 2 + 2], self.U[m[1] * 2],
    #                            self.U[m[1] * 2 + 1], self.U[m[1] * 2 + 2]])
    #     self.U_m_local = np.array([self.T[i] @ self.U_m_global[i] for i in range(len(self.members))])
    #
    # def calcInternalForces(self):
    #     self.F_m_internal = np.array([k @ u for k, u in zip(self.K_m_local, self.U_m_local)])

# Nodes = [[0,0],[1,0],[1,1],[0,1]]
# Members = [[0,1],[1,2],[2,3],[3,0],[0,2]]
# Loads = [[2,5,0,0],[3,0,2,0]]
# Supports = [[0,True,True,False],[1,False,True,False]]
# Releases = [[0,0,0,0,0,0,0],[1,0,0,1,0,0,1],[2,0,0,0,0,0,0],[3,0,0,1,0,0,1],[4,0,0,1,0,0,1]]
# MemberGroups_set = [0,0,0,0,0]
# E_set = [1]
# A_set = [1]
# I_set = [1]
#
# Frame = Frame2D()
# Frame.setGeom(Nodes, Members, Loads, Supports, Releases)
# Frame.setA(A_set)
# Frame.setE(E_set)
# Frame.setI(I_set)
# Frame.setMemberGroups(MemberGroups_set)
# Frame.solveLinear()
# Frame.calcMemberDeflections()
# Frame.calcInternalForces()
# print(Frame.U)
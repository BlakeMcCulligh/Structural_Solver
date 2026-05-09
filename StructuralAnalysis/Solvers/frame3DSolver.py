import math

import numpy as np

from CrossSectionOptimization import frame3DOptimizer
from StructuralAnalysis.CrossSectionCalculaters import SquareHSS, RectHSS, TubeHSS, Angle


class Frame3D:
    def __init__(self):
        self.nodes = None # [location,y, z]
        self.members = None # [Node 1, Node 2]
        self.loads = None # [Node, location, y, z, mx, my, mz]
        self.supports = None # [Node, location, y, location, mx, my, mz]
        self.releases = None # [Node, x_i, y_i, z_i, mx_i, my_i, mz_i, x_j, y_j, z_j, mx_j, my_j, mz_j]
        self.memberGroups = None
        self.A = None # [A] * length Members
        self.E = None # [E] * length Members
        self.G = None
        self.J = None
        self.Iy = None
        self.Iz = None

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
        self.Reactions = None

        self.AFormulas = None
        self.IxFormulas = None
        self.IyFormulas = None
        self.JFormulas = None
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

    def setMaterialProperties(self, E, G):
        self.E = np.array(E)
        self.G = np.array(G)

    def setStiffnesses(self, J, Iy, Iz):
        self.J = np.array(J)
        self.Iy = np.array(Iy)
        self.Iz = np.array(Iz)

    def calcLengths(self):
        dx = self.nodes[self.members[:,1],0] - self.nodes[self.members[:,0],0]
        dy = self.nodes[self.members[:,1],1] - self.nodes[self.members[:,0],1]
        dz = self.nodes[self.members[:,1],2] - self.nodes[self.members[:,0],2]
        self.L = np.sqrt(dx**2 + dy**2 + dz**2)

    def calcLocalStiffnessMatrices(self):
        self.K_m_local = []
        for i in range(len(self.members)):
            A = self.A[self.memberGroups[i]]
            E = self.E[self.memberGroups[i]]
            G = self.G[self.memberGroups[i]]
            J = self.J[self.memberGroups[i]]
            Iy = self.Iy[self.memberGroups[i]]
            Iz = self.Iz[self.memberGroups[i]]
            L = self.L[i]
            self.K_m_local.append(np.array([
                [A * E / L, 0, 0, 0, 0, 0, -A * E / L, 0, 0, 0, 0, 0],
                [0, 12 * E * Iz / L ** 3, 0, 0, 0, 6 * E * Iz / L ** 2, 0, -12 * E * Iz / L ** 3, 0, 0, 0, 6 * E * Iz / L ** 2],
                [0, 0, 12 * E * Iy / L ** 3, 0, -6 * E * Iy / L ** 2, 0, 0, 0, -12 * E * Iy / L ** 3, 0, -6 * E * Iy / L ** 2, 0],
                [0, 0, 0, G * J / L, 0, 0, 0, 0, 0, -G * J / L, 0, 0],
                [0, 0, -6 * E * Iy / L ** 2, 0, 4 * E * Iy / L, 0, 0, 0, 6 * E * Iy / L ** 2, 0, 2 * E * Iy / L, 0],
                [0, 6 * E * Iz / L ** 2, 0, 0, 0, 4 * E * Iz / L, 0, -6 * E * Iz / L ** 2, 0, 0, 0, 2 * E * Iz / L],
                [-A * E / L, 0, 0, 0, 0, 0, A * E / L, 0, 0, 0, 0, 0],
                [0, -12 * E * Iz / L ** 3, 0, 0, 0, -6 * E * Iz / L ** 2, 0, 12 * E * Iz / L ** 3, 0, 0, 0, -6 * E * Iz / L ** 2],
                [0, 0, -12 * E * Iy / L ** 3, 0, 6 * E * Iy / L ** 2, 0, 0, 0, 12 * E * Iy / L ** 3, 0, 6 * E * Iy / L ** 2, 0],
                [0, 0, 0, -G * J / L, 0, 0, 0, 0, 0, G * J / L, 0, 0],
                [0, 0, -6 * E * Iy / L ** 2, 0, 2 * E * Iy / L, 0, 0, 0, 6 * E * Iy / L ** 2, 0, 4 * E * Iy / L, 0],
                [0, 6 * E * Iz / L ** 2, 0, 0, 0, 2 * E * Iz / L, 0, -6 * E * Iz / L ** 2, 0, 0, 0, 4 * E * Iz / L]]))

    def applyRelecesTOLocalStiffnessMatrices(self):
        for MemberReleases in self.releases:
            index = MemberReleases[0]
            k = self.K_m_local[index]
            MemberReleases = MemberReleases[1:]
            R1_indices = []
            R2_indices = []
            for i in range(12):
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
        T_list = []
        for i, m in enumerate(self.members):
            Xi = self.nodes[m[0], 0]
            Xj = self.nodes[m[1], 0]
            Yi = self.nodes[m[0], 1]
            Yj = self.nodes[m[1], 1]
            Zi = self.nodes[m[0], 2]
            Zj = self.nodes[m[1], 2]
            x = [(Xj - Xi) / self.L[i], (Yj - Yi) / self.L[i], (Zj - Zi) / self.L[i]]
            # Vertical Members
            if math.isclose(Xi, Xj) and math.isclose(Zi, Zj):
                if Yj > Yi:
                    y = [-1, 0, 0]
                    z = [0, 0, 1]
                else:
                    y = [1, 0, 0]
                    z = [0, 0, 1]
            # Horizontal Members
            elif math.isclose(Yi, Yj):
                y = [0, 1, 0]
                z = np.cross(x, y)
                z = np.divide(z, (z[0] ** 2 + z[1] ** 2 + z[2] ** 2) ** 0.5)
            # Members neither vertical or horizontal
            else:
                proj = [Xj - Xi, 0, Zj - Zi]
                if Yj > Yi:
                    z = np.cross(proj, x)
                else:
                    z = np.cross(x, proj)
                z = np.divide(z, (z[0] ** 2 + z[1] ** 2 + z[2] ** 2) ** 0.5)
                y = np.cross(z, x)
                y = np.divide(y, (y[0] ** 2 + y[1] ** 2 + y[2] ** 2) ** 0.5)
            dirCos = np.array([x, y, z])
            T = np.zeros((12, 12))
            T[0:3, 0:3] = dirCos
            T[3:6, 3:6] = dirCos
            T[6:9, 6:9] = dirCos
            T[9:12, 9:12] = dirCos
            T_list.append(T)
        self.T = np.array(T_list)
        self.T_T = np.array([t.T for t in self.T])

    def calcTransfromLocalToGlobalStiffnessMatrix(self):
        self.K_m_global = self.T_T @ self.K_m_local @ self.T

    def calcGlobalStiffnessMatrix(self):
        self.calcLocalStiffnessMatrices()
        self.applyRelecesTOLocalStiffnessMatrices()
        self.calcTransformationMatrices()
        self.calcTransfromLocalToGlobalStiffnessMatrix()

    def assembleGlobalStiffnessMatrix(self):
        K_global = np.zeros((len(self.nodes) * 6, len(self.nodes) * 6)).tolist()
        for index, m in enumerate(self.members):
            global_dofs = [6 * m[0], 6 * m[0] + 1, 6 * m[0] + 2, 6 * m[0] + 3, 6 * m[0] + 4, 6 * m[0] + 5, 6 * m[1], 6 * m[1] + 1, 6 * m[1] + 2, 6 * m[1] + 3, 6 * m[1] + 4, 6 * m[1] + 5]
            for i in range(12):
                for j in range(12):
                    K_global[int(global_dofs[i])][int(global_dofs[j])] += self.K_m_global[index][i, j]
        self.K = np.array(K_global)

    def calcDegressOfFreedom(self):
        self.fixedDOF = []
        for support in self.supports:
            if support[1]:
                self.fixedDOF.append(support[0] * 6)
            if support[2]:
                self.fixedDOF.append(support[0] * 6 + 1)
            if support[3]:
                self.fixedDOF.append(support[0] * 6 + 2)
            if support[4]:
                self.fixedDOF.append(support[0] * 6 + 3)
            if support[5]:
                self.fixedDOF.append(support[0] * 6 + 4)
            if support[6]:
                self.fixedDOF.append(support[0] * 6 + 5)
        self.freeDOF = list(set(range(len(self.nodes) * 6)) - set(self.fixedDOF))

    def calcForces(self):
        F = np.zeros(len(self.nodes) * 6)
        for load in self.loads:
            F[load[0] * 6] = load[1]
            F[load[0] * 6 + 1] = load[2]
            F[load[0] * 6 + 2] = load[3]
            F[load[0] * 6 + 3] = load[4]
            F[load[0] * 6 + 4] = load[5]
            F[load[0] * 6 + 5] = load[6]
        self.F = F

    def calcDeflections(self):
        Kff = self.K[np.ix_(self.freeDOF, self.freeDOF)]
        Uf = np.linalg.solve(Kff, self.F[self.freeDOF])
        self.U = np.zeros(len(self.nodes) * 6)
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

    def optimizeSolve(self, X):
        self.A = []
        self.Iz = []
        self.Iy = []
        self.J = []

        j = 0
        for i in range(max(self.memberGroups)+1):
            if self.numVariables[i] == 2:
                self.A.append(self.AFormulas[i](X[j], X[j+1]))
                self.Iz.append(self.IxFormulas[i](X[j], X[j+1]))
                self.Iy.append(self.IyFormulas[i](X[j], X[j+1]))
                self.J.append(self.JFormulas[i](X[j], X[j+1]))
                j += 2
            elif self.numVariables[i] == 3:
                self.A.append(self.AFormulas[i](X[j], X[j+1], X[j+2]))
                self.Iz.append(self.IxFormulas[i](X[j], X[j+1], X[j+2]))
                self.Iy.append(self.IyFormulas[i](X[j], X[j+1], X[j+2]))
                self.J.append(self.JFormulas[i](X[j], X[j+1], X[j+2]))
                j += 3
            else:
                print("Error: incorect number of variables")

        self.A = np.array(self.A)
        self.Iz = np.array(self.Iz)
        self.Iy = np.array(self.Iy)
        self.J = np.array(self.J)

        self.calcLocalStiffnessMatrices()
        self.calcTransfromLocalToGlobalStiffnessMatrix()
        self.assembleGlobalStiffnessMatrix()
        self.calcDeflections()

    def optimize(self, crossSections: list, initalGuess: list):

        self.AFormulas = []
        self.IxFormulas = []
        self.IyFormulas = []
        self.JFormulas = []
        self.numVariables = []

        minBounds = []
        maxBounds = []

        for crossSection in crossSections:
            minBounds = minBounds + crossSection.minBounds
            maxBounds = maxBounds + crossSection.maxBounds

            if crossSection._type == "SquareHSS":
                self.AFormulas.append(SquareHSS.getA)
                self.IxFormulas.append(SquareHSS.getI)
                self.IyFormulas.append(SquareHSS.getI)
                self.JFormulas.append(SquareHSS.getJ)
                self.numVariables.append(2)
            elif crossSection._type == "RectHSS":
                self.AFormulas.append(RectHSS.getA)
                self.IxFormulas.append(RectHSS.getIx)
                self.IyFormulas.append(RectHSS.getIy)
                self.JFormulas.append(RectHSS.getJ)
                self.numVariables.append(3)
            elif crossSection._type == "TubeHSS":
                self.AFormulas.append(TubeHSS.getA)
                self.IxFormulas.append(TubeHSS.getI)
                self.IyFormulas.append(TubeHSS.getI)
                self.JFormulas.append(TubeHSS.getJ)
                self.numVariables.append(2)
            elif crossSection._type == "Angle":
                self.AFormulas.append(Angle.getA)
                self.IxFormulas.append(Angle.getIx)
                self.IyFormulas.append(Angle.getIy)
                self.JFormulas.append(Angle.getJ)
                self.numVariables.append(3)

        self.calcLengths()
        self.calcTransformationMatrices()
        self.calcDegressOfFreedom()
        self.calcForces()

        initalGuess = np.array(initalGuess)

        areas = frame3DOptimizer.optimize(self, [minBounds, maxBounds], initalGuess)

        return areas
import numpy as np

from CrossSectionOptimization import frame2DOptimizer
from frame3DSolver.cross_section_calculaters import square_hss, tube_hss, angle, rect_hss


class Frame2D:
    def __init__(self):
        """
        Handels the large computations for 2D frames.
        """
        self.nodes = None # [location,y]
        self.members = None # [Node 1, Node 2]
        self.loads = None # [Node, location, y, m]
        self.supports = None # [Node, location, y, m]
        self.releases = None # [Node, x_i, y_i, m_i, x_j, y_j, m_j]
        self.memberGroups = None
        self.A = None
        self.E = None
        self.I = None

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
        """
        Calculates the lengths of the members.
        """
        dx = self.nodes[self.members[:,1],0] - self.nodes[self.members[:,0],0]
        dy = self.nodes[self.members[:,1],1] - self.nodes[self.members[:,0],1]
        self.L = np.sqrt(dx**2 + dy**2)

    def calcLocalStiffnessMatrices(self):
        """
        Calculates the local stiffness matrices of the members.
        """
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
        """
        Applies the releases to the local stiffness matrices.
        """
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
        """
        Calculates the transformation matrices for the members.
        """
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
        """
        Converts the local stiffness matrices to global stiffness matrices for the members
        :return:
        """
        self.K_m_global = self.T_T @ self.K_m_local @ self.T

    def calcGlobalStiffnessMatrix(self):
        """
        Calcs the global stiffness matrices for the members from strart to finish.
        """
        self.calcLocalStiffnessMatrices()
        self.applyRelecesTOLocalStiffnessMatrices()
        self.calcTransformationMatrices()
        self.calcTransfromLocalToGlobalStiffnessMatrix()

    def assembleGlobalStiffnessMatrix(self):
        """
        Assmebles the global stiffness matrix for the structure
        """
        K_global = np.zeros((len(self.nodes) * 3, len(self.nodes) * 3)).tolist()
        for index, m in enumerate(self.members):
            global_dofs = [3 * m[0], 3 * m[0] + 1, 3 * m[0] + 2, 3 * m[1], 3 * m[1] + 1, 3 * m[1] + 2]
            for i in range(6):
                for j in range(6):
                    K_global[int(global_dofs[i])][int(global_dofs[j])] += self.K_m_global[index][i, j]
        self.K = np.array(K_global)

    def calcDegressOfFreedom(self):
        """
        Assembles the list of the free degrees of freedom for the structure
        """
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
        """
        Assembles the vector of forces aplyed to the structure
        """
        F = np.zeros(len(self.nodes) * 3)
        for load in self.loads:
            F[load[0] * 3] = load[1]
            F[load[0] * 3 + 1] = load[2]
            F[load[0] * 3 + 2] = load[3]
        self.F = F

    def calcDeflections(self):
        """
        Calculates the deflections of each of the degrees of freedom for the structure
        """
        Kff = self.K[np.ix_(self.freeDOF, self.freeDOF)]
        Uf = np.linalg.solve(Kff, self.F[self.freeDOF])
        self.U = np.zeros(len(self.nodes) * 3)
        self.U[self.freeDOF] = Uf

    def calcReactions(self):
        """
        Calculates the Reactions for the structure.
        """
        self.Reactions = self.F
        self.Reactions[self.fixedDOF] = self.K[self.fixedDOF, :] @ self.U

    def solveLinear(self):
        """
        solves the structure linearly using the stiffness methodd.
        """
        self.calcLengths()
        self.calcGlobalStiffnessMatrix()
        self.assembleGlobalStiffnessMatrix()
        self.calcDegressOfFreedom()
        self.calcForces()
        self.calcDeflections()
        self.calcReactions()

    def optimizeSolve(self, X):
        """
        called from the optimizer to solve the frame
        :param X: The optimization variables
        """
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
        """
        Runs the optimization of the frame.
        :param crossSections: list of cross-section objects
        :param initalGuess: the initial guess for the optimization
        :return Results
        """

        self.AFormulas = []
        self.IFormulas = []
        self.numVariables = []

        minBounds = []
        maxBounds = []

        for crossSection in crossSections:
            minBounds = minBounds + crossSection.minBounds
            maxBounds = maxBounds + crossSection.maxBounds

            if crossSection._type == "SquareHSS":
                self.AFormulas.append(SquareHSS.get_A)
                self.IFormulas.append(SquareHSS.get_I)
                self.numVariables.append(2)
            elif crossSection._type == "RectHSS":
                self.AFormulas.append(RectHSS.get_A)
                self.IFormulas.append(RectHSS.get_Ix)
                self.numVariables.append(3)
            elif crossSection._type == "TubeHSS":
                self.AFormulas.append(TubeHSS.get_A)
                self.IFormulas.append(TubeHSS.get_I)
                self.numVariables.append(2)
            elif crossSection._type == "Angle":
                self.AFormulas.append(Angle.get_A)
                self.IFormulas.append(Angle.get_Ix)
                self.numVariables.append(3)

        self.calcLengths()
        self.calcTransformationMatrices()
        self.calcDegressOfFreedom()
        self.calcForces()

        initalGuess = np.array(initalGuess)

        results = frame2DOptimizer.optimize(self, [minBounds, maxBounds], initalGuess)

        return results
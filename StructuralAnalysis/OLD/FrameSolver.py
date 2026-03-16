import numpy as np
from matplotlib import pyplot as plt

class Members:
    def __init__(self, Nodes, E, A, I):
        self.Nodes = Nodes
        self.E = E
        self.A = A
        self.I = I


def solveUF(K, f, u, JointFree, JointSupported):

    # solves the global stiffness equation to get the deformations, and the forces acting on the members
    # K : global stiffness matrix
    # f : forces acting on the members
    # u : deformations
    # JointFree : joints that are not restrained
    # JointSupported : joints that are restrained
    # noinspection SpellCheckingInspection
    print(K[np.ix_(JointFree, JointFree)])
    print(f[JointFree] - K[np.ix_(JointFree, JointSupported)] @ u[JointSupported])
    u[JointFree] = np.linalg.solve(K[np.ix_(JointFree, JointFree)], f[JointFree] - K[np.ix_(JointFree, JointSupported)] @ u[JointSupported])
    f[JointSupported] = K[JointSupported, :].dot(u)
    return u, f


def FrameElement2D(nodexy, E, A, I):

    # Creates a 4x4 element stiffness matrix, in the global coordinate system
    # input:
    # nodexy : [[x1,y1],[x2,y2]]
    # E : Young's modulus
    # A : Cross-section Area
    # I : Second moment of Area

    E1 = np.array([(nodexy[1][0] - nodexy[0][0]), (nodexy[1][1] - nodexy[0][1])])
    L = float(np.linalg.norm(E1))
    E1 = E1 / L
    E2 = np.array([-E1[1], E1[0]])
    Kel_bend = np.array(
        [[12 * E * I / (L ** 3), 6 * E * I / (L ** 2), -12 * E * I / (L ** 3), 6 * E * I / (L ** 2)],
         [6 * E * I / (L ** 2), 4 * E * I / L, -6 * E * I / (L ** 2), 2 * E * I / L],
         [-12 * E * I / (L ** 3), -6 * E * I / (L ** 2), 12 * E * I / (L ** 3), -6 * E * I / (L ** 2)],
         [6 * E * I / (L ** 2), 2 * E * I / L, -6 * E * I / (L ** 2), 4 * E * I / L]])
    Kel_axial = E * A / L * np.array([[1, -1], [-1, 1]])

    Kel_LOC = np.zeros(shape=(6, 6))
    Kel_LOC[0][0] = Kel_axial[0][0]
    Kel_LOC[0][3] = Kel_axial[0][1]
    Kel_LOC[3][0] = Kel_axial[1][0]
    Kel_LOC[3][3] = Kel_axial[1][1]
    for i in range(4):
        for j in range(4):
            if i == 0:
                ii = 1
            elif i == 1:
                ii = 2
            elif i == 2:
                ii = 4
            else:
                ii = 5
            if j == 0:
                jj = 1
            elif j == 1:
                jj = 2
            elif j == 2:
                jj = 4
            else:
                jj = 5
            Kel_LOC[ii][jj] = Kel_bend[i][j]
    Qrot = [[E1[0], E1[1], 0], [E2[0], E2[1], 0], [0, 0, 1]]
    TmUP = np.append(Qrot, np.zeros(shape=(3, 3)), axis=1)
    TmDown = np.append(np.zeros(shape=(3, 3)), Qrot, axis=1)
    Tmatrix = np.append(TmUP, TmDown, axis=0)
    Kel = Tmatrix.T @ Kel_LOC @ Tmatrix
    return Kel

def MakeGlobalStiffnessMatrix(K, elems, Nel, nodes, E, A, I):

    # Creates the global stiffness matrix
    # K : Global stiffness Matrix
    # elems : elements
    # Nel : number of elements
    # nodes : nodes
    # Nnodes : number of nodes
    # E : Young's Modulus
    # A : cross-section area
    # I : second moment of Area
    # noinspection SpellCheckingInspection

    for i in range(Nel):
        elnodes = [elems[i][0] - 1, elems[i][1] - 1]
        nodexy = [[nodes[elnodes[0]][0], nodes[elnodes[0]][1]], [nodes[elnodes[1]][0], nodes[elnodes[1]][1]]]

        Kel = FrameElement2D(nodexy, E[i], A[i], I[i])

        eldofs = np.concatenate([np.arange(3 * (elnodes[0]), 3 * (elnodes[0] + 1)),
                                 np.arange(3 * (elnodes[1]), 3 * (elnodes[1] + 1))])
        K[np.ix_(eldofs, eldofs)] += Kel
    return K

class Frame:
    def __init__(self):
        self.nodes = []
        self.members = []
        self.Suports = []
        self.Loads = []

        self.Deformation = None
        self.Force = None

        self.mag = 20
        self.numDivs = 20

    def addNode(self, x, y):
        node = [x,y]
        self.nodes.append(node)

    def addMember(self, Node1, Node2, E, A, I):
        member = Members([Node1, Node2], E, A, I)
        self.members.append(member)

    def addSupports(self, Node, Direction):
        # Direction: 1 = x direction, 2 = y direction, 3 = moment
        Suport = [Node,Direction]
        self.Suports.append(Suport)

    def addDistributedLoad(self, Node1, Node2, Magnatude):
        dLoad = [Node1, Node2, Magnatude]
        self.Loads.append(dLoad)

    def solve(self):
        Nel = len(self.members)
        Nnodes = len(self.nodes)

        alldofs = []
        for i in range(3 * Nnodes):
            alldofs.append(i)
        K = np.zeros(shape=(3 * Nnodes, 3 * Nnodes))
        self.Deformation = np.zeros(shape=(3 * Nnodes, 1))
        self.Force = np.zeros(shape=(3 * Nnodes, 1))

        JointsSupported = []
        for i in range(len(self.Suports)):
            thisSupport = 3 * (self.Suports[i][0] - 1) + self.Suports[i][1]
            JointsSupported.append(thisSupport)
            self.Deformation[thisSupport] = 0
        JointsFree = alldofs.copy()
        JointsSupported.sort(reverse=True)
        for index in JointsSupported:
            del JointsFree[index - 1]

        for i in range(len(self.Loads)):
            self.Force[3 * (self.Loads[i][0] - 1) + self.Loads[i][1] - 1] = self.Loads[i][2]

        MemberNodes = []
        E = []
        A = []
        I = []
        for i in range(len(self.members)):
            MemberNodes.append(self.members[i].Nodes)
            E.append(self.members[i].E)
            A.append(self.members[i].A)
            I.append(self.members[i].I)

        K = MakeGlobalStiffnessMatrix(K, MemberNodes, Nel, self.nodes, E, A, I)
        self.Deformation, self.Force = solveUF(K, self.Force, self.Deformation, JointsFree, JointsSupported)

        self.Deformation = self.Deformation
        self.Force = self.Force

    def printResults(self):

        print("Deformation: ", self.Deformation)
        print("Forces in each member: ", self.Force)

        MemberNodes = []
        for i in range(len(self.members)):
            MemberNodes.append(self.members[i].Nodes)

        Nel = len(MemberNodes)
        Nnodes = len(self.nodes)
        # ploting
        # undeformed
        for i in range(Nnodes):
            plt.plot(self.nodes[i][0], self.nodes[i][1], 'o', color='black')
        for i in range(Nel):
            plt.plot([self.nodes[MemberNodes[i][0] - 1][0], self.nodes[MemberNodes[i][1] - 1][0]],
                     [self.nodes[MemberNodes[i][0] - 1][1], self.nodes[MemberNodes[i][1] - 1][1]], color='black')

        # deformed
        for i in range(Nnodes):
            plt.plot(self.nodes[i][0] + self.mag * self.Deformation[i * 3], self.nodes[i][1] + self.mag * self.Deformation[i * 3 + 1], 'o', color='red')

        for i in range(Nel):
            elnodes = [MemberNodes[i][0] - 1, MemberNodes[i][1] - 1]
            E1 = np.array(
                [(self.nodes[elnodes[1]][0] - self.nodes[elnodes[0]][0]), (self.nodes[elnodes[1]][1] - self.nodes[elnodes[0]][1])])
            le = float(np.linalg.norm(E1))
            E1 = E1 / le
            E1 = E1.tolist()
            E2 = [-E1[1], E1[0]]

            eldofs = np.concatenate([np.arange(3 * (elnodes[0]), 3 * (elnodes[0] + 1)),
                                     np.arange(3 * (elnodes[1]), 3 * (elnodes[1] + 1))])
            ut = self.Deformation.tolist()
            eldisp = []
            for j in range(len(eldofs)):
                eldisp.append(ut[eldofs[j]])
            Qrot = [[E1[0], E1[1], 0], [E2[0], E2[1], 0], [0, 0, 1]]
            TmUP = np.append(Qrot, np.zeros(shape=(3, 3)), axis=1)
            TmDown = np.append(np.zeros(shape=(3, 3)), Qrot, axis=1)
            Tmatrix = np.append(TmUP, TmDown, axis=0)
            eldispLOC = Tmatrix @ eldisp
            plotpts = []
            for j in range(self.numDivs + 1):
                xi = j / self.numDivs
                xdispLOC = eldispLOC[0] * (1 - xi) + eldispLOC[3] * xi
                ydispLOC = eldispLOC[1] * (1 - 3 * xi ** 2 + 2 * xi ** 3) + eldispLOC[4] * (3 * xi ** 2 - 2 * xi ** 3) + \
                           eldispLOC[2] * le * (xi - 2 * xi ** 2 + xi ** 3) + eldispLOC[5] * le * (-xi ** 2 + xi ** 3)

                Q = np.array([[Qrot[0][0], Qrot[0][1]], [Qrot[1][0], Qrot[1][1]]])
                xydisp = (Q.T @ np.array([xdispLOC, ydispLOC]))
                x = self.nodes[MemberNodes[i][0] - 1][0] + xi * le * E1[0] + self.mag * xydisp[0]
                y = self.nodes[MemberNodes[i][0] - 1][1] + xi * le * E1[1] + self.mag * xydisp[1]
                plotpts.append([x.tolist()[0], y.tolist()[0]])
            for j in range(len(plotpts) - 1):
                plt.plot([plotpts[j][0], plotpts[j + 1][0]], [plotpts[j][1], plotpts[j + 1][1]], color='red')
        plt.show()

# inputs

F = Frame()
F.addNode(0.0,0.0)
F.addNode(1.0,0.0)
F.addNode(1,1)
F.addNode(0,1.0)

allE = 200000
allA = 0.01
allI = 0.01

F.addMember(1,4, allE, allA, allI)
F.addMember(2,3, allE, allA, allI)
F.addMember(3,4, allE, allA, allI)

F.addSupports(1,1)
F.addSupports(1,2)
F.addSupports(2,1)
F.addSupports(2,2)

F.addDistributedLoad(3,4,-2)

F.solve()
F.printResults()
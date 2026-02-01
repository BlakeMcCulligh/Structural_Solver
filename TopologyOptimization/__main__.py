from math import gcd, ceil
import itertools

from scipy import sparse
from scipy.spatial import Delaunay
import numpy as np
import cvxpy as cvx
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString, Polygon
from time import time

def createNodeGrid(zone, nodeSpasing):
    """
    Crates the grid of nodes to be used in the optimization

    :param zone: A shapely polygon that the truss must be within
    :param nodeSpasing: How far apart the nodes should be in the x and y directions
    :return: array of nodes, number of columbes of nodes, number of rows of nodes
    """

    minX, minY, maxX, maxY = zone.bounds
    dx, dy = maxX - minX, maxY - minY
    numColumbs, numRows = dx / nodeSpasing[0], dy / nodeSpasing[1]

    xv, yv = np.meshgrid(range(int(numColumbs) + 1), range(int(numRows) + 1))

    pts = [Point(xv.flat[i], yv.flat[i]) for i in range(xv.size)]

    # only add nodes that are inside the zone
    Nodes = np.array([[pt.x, pt.y] for pt in pts if zone.intersects(pt)])

    return Nodes, numColumbs, numRows

def assignLoadsAndSupports(Nodes, loadCasses, supports):
    """
    Makes the vectors for the forces and deflections to be used in sitffness matrices

    :param Nodes: Array of nodes
    :param loadCasses: list of load casses: a list of loads: [x cordanit, y cordanit, x loade, y load]
    :param supports: list of supports: [x cordanit, y cordanit, x support, y support]
    :return: vector of forces, vector of deflections
    """

    f = [np.zeros(len(Nodes) * 2)] * len(loadCasses)
    for j, loadCase in enumerate(loadCasses):
        for load in loadCase:
            for i, nd in enumerate(Nodes):
                if nd[0] == load[0] and nd[1] == load[1]:
                    f[j][2*i] = load[2]
                    f[j][2*i+1] = load[3]

    dof = np.ones((len(Nodes), 2))
    for support in supports:
        for i, nd in enumerate(Nodes):
            if nd[0] == support[0] and nd[1] == support[1]:
                if support[2] and support[3]: dof[i, :] = [0, 0]
                elif support[2]:              dof[i, :] = [1, 0]
                elif support[3]:              dof[i, :] = [0, 1]
                else:                         dof[i, :] = [1, 1]

    dof = np.array(dof).flatten()

    return f, dof

def createInitialStructure(Nodes, zone):
    """
    Creates the list of members and findes what ones should be initaly active

    :param Nodes: Array of nodes
    :param zone: A shapely polygon that the truss must be within
    :return: Array of members: [node 1 index, node 2 index, length, Active Member Boolean]
    """

    convex = True if zone.convex_hull.area == zone.area else False

    Members = []

    tri = Delaunay(Nodes)
    for simplex in tri.simplices:
        for i in range(3):
            line = tuple(sorted((simplex[i], simplex[(i + 1) % 3])))
            dx, dy = abs(Nodes[line[0]][0] - Nodes[line[1]][0]), abs(Nodes[line[0]][1] - Nodes[line[1]][1])
            seg = [] if convex else LineString([Nodes[line[0]], Nodes[line[1]]])

            if convex or zone.contains(seg) or zone.boundary.contains(seg):
                Members.append([line[0], line[1], np.sqrt(dx ** 2 + dy ** 2), True])

    Members = np.array(Members)
    A1 = np.zeros((len(Nodes), len(Nodes)))

    A1[Members[:, 0].astype(int) , Members[:, 1].astype(int) ] = 1
    for k in range(5):
        A1 = A1 @ A1
        np.fill_diagonal(A1, 0)
        A1[A1 > 0] = 1
        A1[A1 != 1] = 0

        A1 = A1 + A1.T
        A1[A1 > 0] = 1

    Members = []
    for i in range(len(Nodes)):
        for j in range(len(Nodes)):
            if i != j and A1[i,j] == 1:
                dx, dy = abs(Nodes[j][0] - Nodes[i][0]), abs(Nodes[j][1] - Nodes[i][1])
                Members.append([i, j, np.sqrt(dx ** 2 + dy ** 2), True])

    Members = np.array(Members)

    # setting active membes
    # for pm in [p for p in Members if p[2] <= 1.500]:
    #      pm[3] = True

    return Members

def simplify(Nodes, Members, f, dof):

    ParTolerance1 = 0.01

    ChangeMaster = True
    count = 0
    while ChangeMaster:
        if count > 10:
             break
        count += 1
        ChangeMaster = False
        print("looping")

        ActiveMembers = Members[Members[:, 3] == True]
        vol, a, q, u = solveOptimumProblem(Nodes, ActiveMembers, f, dof)
        #plotTruss(Nodes, ActiveMembers, a, max(a) * 1e-3)
        minA = max(a) * 1e-3

        V = [Nodes[Members[:, 1].astype(int), 0] - Nodes[Members[:, 0].astype(int), 0], Nodes[Members[:, 1].astype(int), 1] - Nodes[Members[:, 0].astype(int), 1]]

        Area = []
        Needed = []
        j = 0

        for i, m1 in enumerate(Members):
            if m1[3]:
                Area.append(a[j])
                if Area[-1] >= minA:
                    Needed.append(True)
                else:
                    Needed.append(False)
                j += 1
            else:
                Area.append(0)
                Needed.append(False)

        V = np.array(V)
        Area = np.array(Area)
        Needed = np.array(Needed)

        M = np.column_stack((Members[:, [0,1]].astype(int), Members[:, [2,3]], V.T, Area, Needed))

        for i, M1 in enumerate(M):
            for j, M2 in enumerate(M):

                if M1[3] and M2[3] and i != j and M1[6] > 1 and M2[6] > 1:

                    C = M1[4] * M2[5] - M1[5] * M2[4]
                    if abs(C) < ParTolerance1:

                        if M1[0] == M2[0]:
                            Change = parSharedPoint(Nodes, ActiveMembers, M, M1, i, M2, j, 0, 0)
                            if Change: ChangeMaster = True
                        elif M1[1] == M2[1]:
                            Change = parSharedPoint(Nodes, ActiveMembers, M, M1, i, M2, j, 1, 1)
                            if Change: ChangeMaster = True
                        elif M1[1] == M2[0]:
                            Change = parSharedPoint(Nodes, ActiveMembers, M, M1, i, M2, j, 1, 0)
                            if Change: ChangeMaster = True
                        elif M1[0] == M2[1]:
                            Change = parSharedPoint(Nodes, ActiveMembers, M, M1, i, M2, j, 0, 1)
                            if Change: ChangeMaster = True

        Members = M[:, [0,1,2,3]]

    return Members


def parSharedPoint(Nodes, ActiveMembers, M, M1, i, M2, j, sp1i, sp2i):

    minA = max(M[:,6]) * 1e-3

    ParTolerance2 = 0.01
    PerpTolerance = 0.01
    maxLength = 20

    nsp1i, nsp2i = abs(sp1i-1), abs(sp2i-1)
    v3 = [Nodes[int(M1[nsp1i])][0] - Nodes[int(M2[nsp2i])][0], Nodes[int(M1[nsp1i])][1] - Nodes[int(M2[nsp2i])][1]]
    C1 = M1[4] * v3[1] - M1[5] * v3[0]
    C2 = M2[4] * v3[1] - M2[5] * v3[0]
    D1 = M1[4] * v3[0] - M1[5] * v3[1]
    D2 = M2[4] * v3[0] - M2[5] * v3[1]

    # if a member from the ends that are not shared is also paralel (with tolerance)
    if abs(C1) < ParTolerance2 and abs(C2) < ParTolerance2:

        sharedPointI = M1[sp1i]
        sharedWithOther = False

        for k, M3 in enumerate(M):
            if (M3[0] == sharedPointI or M3[1] == sharedPointI) and (k != i or k != j) and M3[3] and M3[6] > minA:
                C3 = M1[4] * M3[5] - M1[5] * M3[4]
                if abs(C3) > ParTolerance2:
                    sharedWithOther = True

        if not sharedWithOther:
            dx = Nodes[int(M2[nsp2i]), 0] - Nodes[int(M1[nsp1i]), 0]
            dy = Nodes[int(M2[nsp2i]), 1] - Nodes[int(M1[nsp1i]), 1]
            newMember = [M1[nsp1i], M2[nsp2i], np.sqrt( dx**2 + dy**2 ), True, dx, dy, 100, True]

            if newMember[2] <= maxLength:
                print("simplifying 1")
                M[i] = newMember
                M[j][3] = False
                M[j][7] = False
                return True
        #else:

            # if M1[2] <= M2[2]:
            #     m = M[j]
            #     m[nsp2i] = M1[nsp1i]
            #
            #     lnew = ((Nodes[int(m[1])][0] - Nodes[int(m[0])][0]) ** 2 +
            #             (Nodes[int(m[1])][1] - Nodes[int(m[0])][1]) ** 2) **0.5
            #     if lnew < M2[2]:
            #         M[j][nsp2i] = M1[nsp1i]
            #         print("simplifying 4")
            #         return True
            #     else:
            #         return False
            #
            # else:
            #     m = M[i]
            #     m[nsp1i] = M2[nsp2i]
            #
            #     lnew = ((Nodes[int(m[1])][0] - Nodes[int(m[0])][0]) ** 2 +
            #             (Nodes[int(m[1])][1] - Nodes[int(m[0])][1]) ** 2) ** 0.5
            #     if lnew < M1[2]:
            #         M[i][nsp1i] = M1[nsp2i]
            #         print("simplifying 4")
            #         return True
            #     else:
            #         return False

    # if a member form the ends that are not shared is perpindicular (with tolerance)
    elif abs(D1) < PerpTolerance and abs(D2) < PerpTolerance:
        #print("Canadit 2.1")
        if M1[2] <= M2[2]:
            M[i][3] = False # TODO not sure if this will work
            M[i][7] = False
        else:
            M[j][3] = False # TODO not sure if this will work
            M[j][7] = False
        print("simplifying 2")
        return True

    # if a member from the ends is nether parallel or perpindicular (with tolerance)
    else:
        #print("Canadit 3.1")
        if M1[2] <= M2[2]:
            M[j][nsp2i] = M1[nsp1i]
        else:
            M[i][nsp1i] = M2[nsp2i]
        print("simplifying 3")
        return True
    return False

def solveOptimumProblem(Nodes, ActiveMembers, f, dof):
    """
    Solves for the optimum areas for all the acitve members

    :param Nodes: Array of nodes
    :param ActiveMembers: Array of active members: [node 1 index, node 2 index, length, Active Member Boolean]
    :param f: list of forces actcting on each node for each load case
    :param dof: list of degress of freedom for each node
    :return: volume, areas, forces, deformations
    """

    lb = 0.1 / len(f)
    l = ActiveMembers[:, 2]

    # defining variables
    a = cvx.Variable(len(ActiveMembers), name='a', nonneg=True)  # Cross-section areas
    p = [cvx.Variable(len(ActiveMembers), name='p_lc' + str(k)) for k in range(len(f))]  # Axial forces
    q = [cvx.Variable(len(ActiveMembers), name='q_lc' + str(k)) for k in range(len(f))]  # Element elastic energy
    beta = cvx.Variable(len(f), nonneg=True)  # Dual variables of lower bound on case weighting = slacks of eqn

    # setting objective function
    obj = cvx.Minimize(cvx.sum(l @ a) - cvx.sum(lb * beta))

    # setting equilibrium constraints
    B = calcBi(Nodes, ActiveMembers, dof)  # reduced
    equilibCon = [B.transpose() @ q[k] + f[k] == 0 for k in range(len(f))]

    # setting conic constraints
    coneConFlat = []
    invLen = [2 / li for li in l]
    for k in range(len(f)):
        xVals = cvx.vstack([2 * q[k], cvx.multiply(invLen, a) - p[k]])
        tVals = cvx.multiply(invLen, a) + p[k]
        coneConFlat.append(cvx.SOC(tVals, xVals))

    # setting loadcase total energy constraints
    sumCon = [cvx.sum(p[k]) + beta[k] == 0.5 for k in range(len(f))]

    # Solving problem
    prob = cvx.Problem(obj, sumCon + equilibCon + coneConFlat)
    volume = prob.solve(cvx.MOSEK, verbose=False)

    # getting solved values
    deformations = np.array([equilibCon[k].dual_value for k in range(len(f))])
    areas = a.value
    forces = np.array([q[k].value for k in range(len(f))])

    return volume, areas, forces, deformations

# Check dual (kinematic) violation
def stopViolation(Nodes, Members, deformations):
    """
    findes the members that are inactive but should be and makes them active. If there is none,
    return true to end the optimization

    :param Nodes: Array of nodes
    :param Members: Array of members: [node 1 index, node 2 index, length, Active Member Boolean]
    :param deformations: Array of deformations
    :return: If the program is done optimizing
    """

    # check only inactive members for efficiency
    InactiveMemberIndices = np.where(Members[:, 3] == False)[0]

    # if all members active terminate process
    if len(InactiveMemberIndices) == 0: return True

    InactiveMembers = Members[InactiveMemberIndices]
    MemberLengths = InactiveMembers[:, 2]
    B = calcBi(Nodes, InactiveMembers)  # un-reduced version

    #Values >1 violate
    violation = np.zeros(len(InactiveMembers))
    for k, deformation in enumerate(deformations):  # Calculate per-load-case
        gk = B.dot(deformation)
        violation += [0.5 * gk[i] ** 2 / (MemberLengths[i] ** 2) for i in range(len(InactiveMembers))]

    ViolatedInactiveMemberIndices = np.where(violation > 1.0001)[0]

    # indices of most violated elements within ViolatedInactiveMemberIndices
    SortedViolatedInactiveMemberIndices = np.flipud(np.argsort(violation[ViolatedInactiveMemberIndices]))

    # number of elements to be added, max 30% of current size
    numMembersToAdd = ceil(min(float(len(SortedViolatedInactiveMemberIndices)), (len(Members) - len(InactiveMembers)) * 0.3))

    for i in range(numMembersToAdd):
        Members[InactiveMemberIndices[ViolatedInactiveMemberIndices[SortedViolatedInactiveMemberIndices[i]]]][3] = True

    print('        adding', numMembersToAdd, ' of ', len(ViolatedInactiveMemberIndices), ' violated elements')

    # if no added elements, then terminate process
    return numMembersToAdd == 0

def plotTruss(Nodes, ActiveMembers, a, threshold):
    """
    Plots the truss

    :param Nodes: Array of nodes
    :param ActiveMembers: Array of active members: [node 1 index, node 2 index, length, Active Member Boolean]
    :param a: Area assigned to the member
    :param threshold: Required area of the member to be printed
    """

    plt.figure()
    plt.clf()
    plt.axis('equal')
    plt.draw()

    ThinknessMult = 0.08  # line thikness multiplyer

    # ploting members
    for i in [i for i in range(len(a)) if a[i] >= threshold]:

        pos = Nodes[ActiveMembers[i, [0, 1]].astype(int), :]
        plt.plot(pos[:, 0], pos[:, 1], linewidth=np.sqrt(a[i] * ThinknessMult), solid_capstyle='round')

    # ploting nodes
    #plt.plot(Nodes[:, 0], Nodes[:, 1], 'k.')
    plt.show()

def calcBi(Nodes, Members, dof=None):
    """
    Calculates the equilibrium matrix. In reduced form if dof information is given

    :param Nodes: Array of nodes
    :param Members: Array of members ether active or inactive [node 1 index, node 2 index, length, Active Member Boolean]
    :param dof: Array Degrese of freedom for each node
    :return: the equilibrium matrix
    """

    if dof is None: dof = []

    m, n1, n2 = len(Members), Members[:, 0].astype(int), Members[:, 1].astype(int)

    l, dx, dy = Members[:, 2], Nodes[n2, 0] - Nodes[n1, 0], Nodes[n2, 1] - Nodes[n1, 1]

    # If information on fixed degrees of freedom is provided
    if len(dof):
        d0, d1, d2, d3 = dof[n1 * 2], dof[n1 * 2 + 1], dof[n2 * 2], dof[n2 * 2 + 1]
        v = np.concatenate((-dx / l * d0, -dy / l * d1, dx / l * d2, dy / l * d3))

    else:
        v = np.concatenate((-dx / l, -dy / l, dx / l, dy / l))

    r = np.concatenate((n1 * 2, n1 * 2 + 1, n2 * 2, n2 * 2 + 1))

    c = np.concatenate((np.arange(m), np.arange(m), np.arange(m), np.arange(m)))

    return sparse.coo_matrix((v, (c, r)), shape=(m, len(Nodes) * 2))

def OptimizeTruss(zoneNodes, loadCasses, supports, nodeSpasing =None, nodes = None):

    start = time()

    zone = Polygon(zoneNodes)

    if nodes is None:
        Nodes, numColumbs, numRows = createNodeGrid(zone, nodeSpasing)

    else:
        Nodes = np.array(nodes)

    print(Nodes)

    # Load and support conditions
    f, dof = assignLoadsAndSupports(Nodes, loadCasses, supports)

    # Create the initial structure
    Members = createInitialStructure(Nodes, zone)

    print('Nodes: %d Members: %d' % (len(Nodes), len(Members)))

    plotTruss(Nodes, Members, np.ones((len(Members))) * 30, 0)
    # vol, a, q, u = solveOptimumProblem(Nodes, Members, f, dof)
    # plotTruss(Nodes, Members, a, max(a) * 1e-3)

    # Start the main member adding loop
    itr, vol, ActiveMembers, a, q = 0, 0, 0, 0, 0
    for itr in range(1, 5):

        Members = simplify(Nodes, Members, f, dof)

        ActiveMembers = Members[Members[:, 3] == True]# Only take members that are active

        vol, a, q, u = solveOptimumProblem(Nodes, ActiveMembers, f, dof)

        print("Itr: %d, volume: %f, active members: %d" % (itr, vol, len(ActiveMembers)))

        plotTruss(Nodes, ActiveMembers, a, max(a) * 1e-3)# plot interation for truss
        #plotTruss(Nodes, ActiveMembers, np.ones((len(a)))*30, 0)

        if stopViolation(Nodes, Members, u):
            break  # if no elements are added, terminate process

    end = time()
    print('Time taken ', end - start)

    plotTruss(Nodes, ActiveMembers, a, max(a) * 1e-3)

    return [vol, a, end - start, itr, len(ActiveMembers), len(Members)]


if __name__ == '__main__':
    ZoneNodes = [[0, 0], [0, 10], [10, 10], [10, 0]]
    Supports = [[0,0,True, True],[0,10, True, True]]
    LoadCasses = [[[10,0,0, -5]],[[6,5,-5,0]]]
    MaxSpassing = [5, 5]

    # ZoneNodes = [[0, 0], [12, 0], [12, 5], [78, 5], [78, 12], [162, 12], [162,5], [228,5], [228,0], [240, 0], [240,12], [288, 12], [288,28], [0,28]]
    # Supports = [[0, 0, True, True], [2, 0, True, True], [4, 0, True, True], [6, 0, True, True], [8, 0, True, True], [10, 0, True, True], [12, 0, True, True], [228, 0, True, True], [230, 0, True, True], [232, 0, True, True], [234, 0, True, True], [236, 0, True, True], [238, 0, True, True], [240, 0, True, True]]
    #
    # for i in range(len(Supports)):
    #     ZoneNodes[i] = [ZoneNodes[i][0]/4, ZoneNodes[i][1]/4]
    #
    # for i in range(len(Supports)):
    #     Supports[i] = [Supports[i][0]/4, Supports[i][1]/4, True, True]
    #
    # LoadCasses = [[[10, 28/4, 0, -1/7], [20, 28/4, 0, -1/7], [30, 28/4, 0, -1/7], [40, 28/4, 0, -1/7], [50, 28/4, 0, -1/7], [60, 28/4, 0, -1/7], [72, 28/4, 0, -1/7]]]
    # #for i in range(7, 50):
    #     LoadCasses[0].append([i, 28/4, 0, -1/108])
    # for i in range(60, 72):
    #     LoadCasses[0].append([i, 28/4, 0, -1/108])

    # ZoneNodes = [[0, 0], [0, 5], [20, 5], [20, 0]]
    # Supports = [[0,0,True, True],[20,0, True, True]]
    # LoadCasses = [[[5,5,0, -1], [15,5,0,-1]]]


    #MaxSpassing = [1, 1]

    Vol, Time, A, NumIter, NumActive, NumPotential = OptimizeTruss(
        zoneNodes = ZoneNodes, # the corners of the polygon defining where the truss can be
        nodeSpasing = MaxSpassing, # maximum spasing between nodes in x and y directions
        loadCasses = LoadCasses, # loads [[[x cord, y cord, x load, y load]]]
        supports = Supports)  # nodes [[x cord, y cord, x fixarity, y fixarity]]
    plt.show()
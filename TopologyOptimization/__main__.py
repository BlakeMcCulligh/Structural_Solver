import math
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

    TotalLoad = 0
    for case in loadCasses:
        for load in case:
            TotalLoad += abs(load[2]) + abs(load[3])
    if TotalLoad == 0:
        TotalLoad = 1

    f = [np.zeros(len(Nodes) * 2)] * len(loadCasses)
    for j, loadCase in enumerate(loadCasses):
        for load in loadCase:
            for i, nd in enumerate(Nodes):
                if nd[0] == load[0] and nd[1] == load[1]:
                    f[j][2*i] = load[2]/TotalLoad
                    f[j][2*i+1] = load[3]/TotalLoad

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
                seg = [] if convex else LineString([Nodes[i], Nodes[j]])
                if convex or zone.contains(seg) or zone.boundary.contains(seg):
                    if np.sqrt(dx ** 2 + dy ** 2) <= 42:
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

    y = cvx.Variable(len(ActiveMembers), boolean=True)


    ActveMCon = [a - 0.01 <= 1000 * y, a >= 0, cvx.sum(y) <= 75]

    obj = cvx.Minimize((cvx.sum(l @ a) - cvx.sum(lb * beta)))
    #sum(l[:] > 6)
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

    #memberOverlaps = find_overlaps(ActiveMembers[:,[0,1]], Nodes)

    #overlapCon = [(memberOverlaps @ a) >= 0.05]

    # Solving problem
    print("optimizing")
    prob = cvx.Problem(obj, sumCon + equilibCon + coneConFlat + ActveMCon)
    #prob = cvx.Problem(obj, sumCon + equilibCon + coneConFlat + overlapCon)
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
    #plt.show()

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


def find_overlaps(members, Nodes):
    overlaps = []

    def get_orientation(p, q, r):
        val = (q[1] - p[1]) * (r[0] - q[0]) - (q[0] - p[0]) * (r[1] - q[1])
        return 0 if val == 0 else (1 if val > 0 else 2)

    def on_segment(p, q, r):
        return min(p[0], r[0]) <= q[0] <= max(p[0], r[0]) and \
            min(p[1], r[1]) <= q[1] <= max(p[1], r[1])

    def do_intersect(s1, s2):
        p1, q1 = Nodes[int(s1[0])], Nodes[int(s1[1])]
        p2, q2 = Nodes[int(s2[0])], Nodes[int(s2[1])]
        o1, o2 = get_orientation(p1, q1, p2), get_orientation(p1, q1, q2)
        o3, o4 = get_orientation(p2, q2, p1), get_orientation(p2, q2, q1)

        if o1 != o2 and o3 != o4: return True
        if o1 == 0 and on_segment(p1, p2, q1): return True
        if o2 == 0 and on_segment(p1, q2, q1): return True
        if o3 == 0 and on_segment(p2, p1, q2): return True
        if o4 == 0 and on_segment(p2, q1, q2): return True
        return False

    overlaps = np.zeros((len(members), len(members)))

    for i in range(len(members)):
        for j in range(len(members)):
            if i != j and do_intersect(members[i], members[j]):
                overlaps[i,j] = 1
                overlaps[j,i] = 1
    return overlaps


def is_on_line_segment(points, p1, p2, tolerance=1e-6):
    """
    Checks which points from a list lie on the line segment defined by p1 and p2.

    Args:
        points (np.ndarray): An array of shape (N, 2) or (N, 3) for N points.
        p1 (np.ndarray): The start point of the line segment.
        p2 (np.ndarray): The end point of the line segment.
        tolerance (float): A tolerance for floating-point comparisons.

    Returns:
        np.ndarray: A boolean array indicating if each point is on the segment.
    """

    p1 = np.array(p1)
    p2 = np.array(p2)

    # Vector from p1 to the point
    vec_p1_to_points = points - p1
    # Vector from p1 to p2 (the segment vector)
    vec_segment = p2 - p1

    # 1. Check if the points are collinear (cross product is near zero)
    # For 2D points, the "cross product" is a scalar value
    if points.shape[1] == 2:
        cross_products = np.cross(vec_segment, vec_p1_to_points)
        is_collinear = np.abs(cross_products) < tolerance
    # For 3D points, the cross product is a vector, check its magnitude
    elif points.shape[1] == 3:
        cross_products = np.cross(np.tile(vec_segment, (points.shape[0], 1)), vec_p1_to_points)
        # Check if the magnitude of each resulting vector is near zero
        is_collinear = np.linalg.norm(cross_products, axis=1) < tolerance
    else:
        raise ValueError("Points must be 2D or 3D.")

    # 2. Check if the points are between p1 and p2 (dot product condition)
    # The dot product of vec_p1_to_points and vec_segment should be > 0 and < length of segment squared
    dot_products = np.sum(vec_p1_to_points * vec_segment, axis=1)
    segment_len_sq = np.sum(vec_segment * vec_segment)

    # Ensure dot product is positive (past p1) and less than segment length squared (before p2)
    is_within_bounds = (dot_products >= -tolerance) & (dot_products <= segment_len_sq + tolerance)

    # A point is on the segment if both conditions are true
    return is_collinear & is_within_bounds


def sort_nodes_along_line(nodes, start_point, end_point):
    """
    Sorts a list of nodes based on their position along the line
    defined by start_point and end_point.

    Args:
        nodes (list of tuples/arrays): The points to sort (e.g., [(x1, y1), ...]).
        start_point (tuple/array): The starting point of the reference line.
        end_point (tuple/array): The ending point of the reference line.

    Returns:
        list of tuples/arrays: The sorted nodes.
    """

    # 1. Calculate the line's direction vector (V)
    # V points from start to end
    direction_vector = end_point - start_point

    # 2. Normalize the direction vector to a unit vector (optional, but clean)
    # This simplifies the projection calculation slightly
    unit_direction = direction_vector / np.linalg.norm(direction_vector)

    # 3. Calculate projection of each node relative to the start point
    # We look at the vector from start_point to the node (P - Start)
    vectors_to_nodes = nodes - start_point

    # The scalar projection is the dot product of (P - Start) and the unit direction
    # This gives us a single scalar value for each node, which we can sort by
    projections = np.dot(vectors_to_nodes, unit_direction)

    # 4. Sort the nodes based on these projection values
    # np.argsort returns the indices that would sort the array
    sorted_indices = np.argsort(projections)

    # Use the sorted indices to reorder the original nodes list
    sorted_nodes = nodes[sorted_indices]

    return sorted_nodes

def constructLoadCases(loadCasses, maxLength, nodes):

    disLoadCasses = []

    for i, loadCase in enumerate(loadCasses):
        disLoadCasse = []
        for j, load in enumerate(loadCase):
            if load[4] is not None and not math.isnan(float(load[4])):
                print(load[4])
                disLoadCasse.append([load[4], load[5], load[6], load[7], load[8], load[9]])
        disLoadCasses.append(disLoadCasse)

    disLoadCassesSplit = []
    for i, disLoadCasse in enumerate(disLoadCasses):
        disLoadCasseSplit = []
        for j, disLoad in enumerate(disLoadCasse):
            p1 = np.array([disLoad[0], disLoad[1]])
            p2 = np.array([disLoad[2], disLoad[3]])
            nodes = np.array(nodes)
            nodesOnSegment = nodes[is_on_line_segment(nodes, p1, p2)]

            nodesOnSegment = sort_nodes_along_line(nodesOnSegment, p1, p2)
            nodesOnSegment = np.unique(nodesOnSegment, axis=0)

            PosibleConbinations = []
            n = len(nodesOnSegment)

            # Backtracking to find subsets of intervals that cover the range
            def backtrack(start_idx, current_combination):
                if start_idx == n - 1:
                    PosibleConbinations.append(list(current_combination))
                    return

                for end_idx in range(start_idx + 1, n):

                    length = np.linalg.norm(nodesOnSegment[start_idx] - nodesOnSegment[end_idx])

                    if length <= maxLength:
                        interval = (nodesOnSegment[start_idx], nodesOnSegment[end_idx], disLoad[4], disLoad[5])
                        backtrack(end_idx, current_combination + [interval])

            backtrack(0, [])
            disLoadCasseSplit.append(PosibleConbinations[0])
        disLoadCassesSplit.append(list(itertools.product(*disLoadCasseSplit)))
    disLoadCassesSplit = list(itertools.product(*disLoadCassesSplit))

    StructureCases = [None] * len(disLoadCassesSplit[0])
    for i, disLoadCasseSplit in enumerate(disLoadCassesSplit):

        for j, disLoadCaseLayout in enumerate(disLoadCasseSplit):

            if i == 0:
                StructureCases[j] = [None] * len(disLoadCassesSplit)

            loads = []
            for k, disLoadSegment in enumerate(disLoadCaseLayout):
                length = math.dist(disLoadSegment[0], disLoadSegment[1])
                Rx, Ry = disLoadSegment[2] * length / 2, disLoadSegment[3] * length / 2
                loads.append([disLoadSegment[0][0], disLoadSegment[0][1], Rx, Ry])
                loads.append([disLoadSegment[1][0], disLoadSegment[1][1], Rx, Ry])

            StructureCases[j][i] = loads

    for StructureCase in StructureCases:
        for i, loadCase in enumerate(loadCasses):
            for load in loadCase:
                if load[0] is not None and not math.isnan(float(load[0])):
                    StructureCase[i] = StructureCase[i] + [[load[0], load[1], load[2], load[3]]]

    return StructureCases

def OptimizeTruss(zoneNodes, loadCasses, supports, nodeSpasing =None, nodes = None, maxLength = None):

    maxLength = 12

    start = time()

    zone = Polygon(zoneNodes)

    if nodes is None:
        Nodes, numColumbs, numRows = createNodeGrid(zone, nodeSpasing)
    else:
        Nodes = np.array(nodes)

    StructureCases = constructLoadCases(loadCasses, maxLength, nodes)

    Members = createInitialStructure(Nodes, zone)
    #plotTruss(Nodes, Members, np.ones((len(Members))) * 30, 0)
    for loadCasses in StructureCases:
        if loadCasses is not None:

            f, dof = assignLoadsAndSupports(Nodes, loadCasses, supports)
            vol, a, q, u = solveOptimumProblem(Nodes, Members, f, dof)

            if a is not None:
                plotTruss(Nodes, Members, a, max(a) * 1e-3)
            else:
                print("Truss Is not stable")

        else:
            print("No Load Case Found")

    plt.show()

    #
    # # Load and support conditions
    # f, dof = assignLoadsAndSupports(Nodes, loadCasses, supports)
    #
    # # Create the initial structure
    # Members = createInitialStructure(Nodes, zone)
    #
    # print('Nodes: %d Members: %d' % (len(Nodes), len(Members)))
    #
    # plotTruss(Nodes, Members, np.ones((len(Members))) * 30, 0)
    # vol, a, q, u = solveOptimumProblem(Nodes, Members, f, dof)
    # plotTruss(Nodes, Members, a, max(a) * 1e-3)
    #
    # # Start the main member adding loop
    # itr, vol, ActiveMembers, a, q = 0, 0, 0, 0, 0
    # for itr in range(1, 5):
    #
    #     Members = simplify(Nodes, Members, f, dof)
    #
    #     ActiveMembers = Members[Members[:, 3] == True]# Only take members that are active
    #
    #     vol, a, q, u = solveOptimumProblem(Nodes, ActiveMembers, f, dof)
    #
    #     print("Itr: %d, volume: %f, active members: %d" % (itr, vol, len(ActiveMembers)))
    #
    #     plotTruss(Nodes, ActiveMembers, a, max(a) * 1e-3)# plot interation for truss
    #     #plotTruss(Nodes, ActiveMembers, np.ones((len(a)))*30, 0)
    #
    #     if stopViolation(Nodes, Members, u):
    #         break  # if no elements are added, terminate process
    #
    # end = time()
    # print('Time taken ', end - start)
    #
    # plotTruss(Nodes, ActiveMembers, a, max(a) * 1e-3)
    #
    # return [vol, a, end - start, itr, len(ActiveMembers), len(Members)]

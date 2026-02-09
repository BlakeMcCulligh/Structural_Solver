import math
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

    return Members

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
    obj = cvx.Minimize((cvx.sum(l @ a) - cvx.sum(lb * beta)))

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
    print("optimizing")
    prob = cvx.Problem(obj, sumCon + equilibCon + coneConFlat)
    volume = prob.solve(cvx.MOSEK, verbose=False)

    # getting solved values
    print("getting results")
    deformations = np.array([equilibCon[k].dual_value for k in range(len(f))])
    areas = a.value
    forces = np.array([q[k].value for k in range(len(f))])

    return volume, areas, forces, deformations

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

def getNodesOnSegment(nodes,p1, p2):
    p1 = np.array(p1)
    p2 = np.array(p2)

    nodesOnSegment = nodes[is_on_line_segment(nodes, p1, p2)]
    nodesOnSegment = sort_nodes_along_line(nodesOnSegment, p1, p2)
    nodesOnSegment = np.unique(nodesOnSegment, axis=0)

    return nodesOnSegment


def constructLoadCases(loadCasses, maxLength, nodes, supports):
    nodes = np.array(nodes)
    #print(loadCasses)

    #TODO ONLY FIST LODE CASE ACTIVE
    loads = loadCasses[0]

    disLoads = []
    for i, load in enumerate(loads):
        if load[4] is not None and not math.isnan(float(load[4])):
            disLoads.append([load[4], load[5], load[6], load[7], load[8], load[9]])

    # disLoads: list of loads, loads: [x1,y1, x2,y2, fx,fy]
    #print("disLoads: ", disLoads)

    posibleCombinationsOverall = []
    for i, load in enumerate(disLoads):
        p1, p2 = [load[0], load[1]], [load[2], load[3]]

        nodesOnSegment = getNodesOnSegment(nodes,p1, p2)

        posibleCombinationsPerDisLoad = []
        n = len(nodesOnSegment)

        # Backtracking to find subsets of intervals that cover the range
        def backtrack(start_idx, current_combination):
            if start_idx == n - 1:
                posibleCombinationsPerDisLoad.append(list(current_combination))
                return

            for end_idx in range(start_idx + 1, n):
                if np.linalg.norm(nodesOnSegment[start_idx] - nodesOnSegment[end_idx]) <= maxLength:
                    interval = (nodesOnSegment[start_idx], nodesOnSegment[end_idx], load[4], load[5])
                    backtrack(end_idx, current_combination + [interval])

        backtrack(0, [])

        # posible combinations: list of structure layouts, Structure layouts: list of loads, loads: [[x1,y1], [x2,y2], fx,fy]
        #print("posibleCombinationsPerDisLoad: ", posibleCombinationsPerDisLoad)
        posibleCombinationsOverall.append(posibleCombinationsPerDisLoad)

    # combining posible structure layouts for each distributed load
    posibleCombinationsOverall = list(itertools.product(*posibleCombinationsOverall))
    for i, posibleCombinationIn in enumerate(posibleCombinationsOverall):
        p = []
        for newloads in posibleCombinationIn:
            for newload in newloads:
                p.append(newload)
        posibleCombinationsOverall[i] = p

    # posibleCombinationsOverall: list of structure layouts, Structure layouts: list of loads, loads: [[x1,y1], [x2,y2], fx,fy]
    #print("posibleCombinationsOverall: ", posibleCombinationsOverall)

    # converting distributed loads into point loads
    posibleStructures = []
    for i, disLoads in enumerate(posibleCombinationsOverall):
        newloads = []
        for j, disLoad in enumerate(disLoads):
            length = math.dist(disLoad[0], disLoad[1])
            Rx, Ry = disLoad[2] * length / 2, disLoad[3] * length / 2
            newloads.append([float(disLoad[0][0]), float(disLoad[0][1]), Rx, Ry])
            newloads.append([float(disLoad[1][0]), float(disLoad[1][1]), Rx, Ry])
        posibleStructures.append(newloads)

    #print("posibleStructures 1: ", posibleStructures)

    # adding point loads
    for i, newloads in enumerate(posibleStructures):
        for load in loads:
            if load[0] is not None and not math.isnan(float(load[0])):
                posibleStructures[i] = posibleStructures[i] + [[load[0], load[1], load[2], load[3]]]

    #print("posibleStructures 2: ", posibleStructures)

    for i, loads in enumerate(posibleStructures):
        newLoads = []
        for load in loads:
            loadFound = False
            for newLoad in newLoads:
                if float(load[0]) == float(newLoad[0]) and float(load[1]) == float(newLoad[1]):
                    newLoad[2] += load[2]
                    newLoad[3] += load[3]
                    loadFound = True
                    break

            if not loadFound:
                for support in supports:
                    if float(support[0]) == float(load[0]) and float(support[1]) == float(load[1]):
                        loadFound = True

            if not loadFound:
                newLoads.append(load)
        posibleStructures[i] = [newLoads] # TODO brakets added as only one load case is posible curently

    #print("posibleStructures 3: ", posibleStructures)

    return posibleStructures

def OptimizeTruss(zoneNodes, loadCasses, supports, nodeSpasing =None, nodes = None, maxLength = None):

    maxLength = 42

    start = time()

    zone = Polygon(zoneNodes)

    if nodes is None:
        Nodes, numColumbs, numRows = createNodeGrid(zone, nodeSpasing)
    else:
        Nodes = np.array(nodes)

    StructureCases = constructLoadCases(loadCasses, maxLength, nodes, supports)

    Members = createInitialStructure(Nodes, zone)
    #plotTruss(Nodes, Members, np.ones((len(Members))) * 30, 0)
    for loadCasses in StructureCases:
        if loadCasses is not None:
            #print(loadCasses)

            f, dof = assignLoadsAndSupports(Nodes, loadCasses, supports)
            vol, a, q, u = solveOptimumProblem(Nodes, Members, f, dof)
            if a is not None:
                #print("Solved")
                plotTruss(Nodes, Members, a, max(a) * 1e-3)
            else:
                print("Truss Is not stable")

        else:
            print("No Load Case Found")
    plt.show()
    print("finished")
    print("Time taken: ", time() - start)


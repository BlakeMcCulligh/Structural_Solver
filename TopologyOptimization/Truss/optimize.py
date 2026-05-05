import math
import itertools

from scipy.spatial import Delaunay
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString, Polygon
from time import time

from CrossSectionOptimization.Truss.optimize import OptimizeTrussCrossSections
from OpeningAndSaving.Saving import saveTrussTopologyOptimizationExcel

def createNodeGrid(zone, nodeSpasing):
    """
    Crates the grid of printNodes to be used in the optimization

    :param zone: A shapely polygon that the truss must be within
    :param nodeSpasing: How far apart the printNodes should be in the location and y directions
    :return: array of printNodes, number of columbes of printNodes, number of rows of printNodes
    """

    minX, minY, maxX, maxY = zone.bounds
    dx, dy = maxX - minX, maxY - minY
    numColumbs, numRows = dx / nodeSpasing[0], dy / nodeSpasing[1]

    xv, yv = np.meshgrid(range(int(numColumbs) + 1), range(int(numRows) + 1))

    pts = [Point(xv.flat[i], yv.flat[i]) for i in range(xv.size)]

    # only add printNodes that are inside the zone
    Nodes = np.array([[pt.x, pt.y] for pt in pts if zone.intersects(pt)])

    return Nodes, numColumbs, numRows

def createInitialStructure(Nodes, zone):
    """
    Creates the list of members and findes what ones should be initaly active

    :param Nodes: Array of printNodes
    :param zone: A shapely polygon that the truss must be within
    :return: Array of members: [node 1 i, node 2 i, length]
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
                Members.append([line[0], line[1], np.sqrt(dx ** 2 + dy ** 2)])

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
                        Members.append([i, j, np.sqrt(dx ** 2 + dy ** 2)])

    Members = np.array(Members)

    return Members

def is_on_line_segment(points, p1, p2, tolerance=1e-6):
    """
    Checks which points from a list lie on the node segment defined by w1 and p2.

    Args:
        points (np.ndarray): An array of shape (N, 2) or (N, 3) for N points.
        p1 (np.ndarray): The start point of the node segment.
        p2 (np.ndarray): The end point of the node segment.
        tolerance (float): A tolerance for floating-point comparisons.

    Returns:
        np.ndarray: A boolean array indicating if each point is on the segment.
    """

    p1 = np.array(p1)
    p2 = np.array(p2)

    # Vector from w1 to the point
    vec_p1_to_points = points - p1
    # Vector from w1 to p2 (the segment vector)
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

    # 2. Check if the points are between w1 and p2 (dot product condition)
    # The dot product of vec_p1_to_points and vec_segment should be > 0 and < length of segment squared
    dot_products = np.sum(vec_p1_to_points * vec_segment, axis=1)
    segment_len_sq = np.sum(vec_segment * vec_segment)

    # Ensure dot product is positive (past w1) and less than segment length squared (before p2)
    is_within_bounds = (dot_products >= -tolerance) & (dot_products <= segment_len_sq + tolerance)

    # A point is on the segment if both conditions are true
    return is_collinear & is_within_bounds


def sort_nodes_along_line(nodes, start_point, end_point):
    """
    Sorts a list of printNodes based on their position along the node
    defined by start_point and end_point.

    Args:
        nodes (list of tuples/arrays): The points to sort (e.g., [(x1, y1), ...]).
        start_point (tuple/array): The starting point of the reference node.
        end_point (tuple/array): The ending point of the reference node.

    Returns:
        list of tuples/arrays: The sorted printNodes.
    """

    # 1. Calculate the node's direction vector (V)
    # V points from start to end
    direction_vector = end_point - start_point

    # 2. Normalize the direction vector to a unit vector (optional, but clean)
    # This simplifies the projection calculation slightly
    unit_direction = direction_vector / np.linalg.norm(direction_vector)

    # 3. Calculate projection of each node relative to the start point
    # We look at the vector from start_point to the node (Magnatude - Start)
    vectors_to_nodes = nodes - start_point

    # The scalar projection is the dot product of (Magnatude - Start) and the unit direction
    # This gives us a single scalar value for each node, which we can sort by
    projections = np.dot(vectors_to_nodes, unit_direction)

    # 4. Sort the printNodes based on these projection values
    # np.argsort returns the indices that would sort the array
    sorted_indices = np.argsort(projections)

    # Use the sorted indices to reorder the original printNodes list
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

    #TODO ONLY FIST LODE CASE CURRENTLY WORKS
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

def OptimizeTruss(filePath, zoneNodes, loadCasses, supports, nodeSpasing =None, nodes = None, maxLength = None):

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
    Vol = []
    A = []
    Q = []
    U = []
    for loadCasses in StructureCases:
        if loadCasses is not None:
            vol, a, q, u = OptimizeTrussCrossSections(Nodes, Members, loadCasses, supports)
            Vol.append(vol)
            A.append(a)
            Q.append(q)
            U.append(u)
        else:
            print("No Load Case Found")
    plt.show()
    print("finished")
    print("Time taken: ", time() - start)

    BestIndex = Vol.index(min(Vol))

    loadCasses = StructureCases[BestIndex]
    areas = A[BestIndex]
    deflections = U[BestIndex]
    forces = Q[BestIndex]
    volume = Vol[BestIndex]


    saveTrussTopologyOptimizationExcel(filePath, nodes, Members, loadCasses, supports, areas, deflections, forces, volume)


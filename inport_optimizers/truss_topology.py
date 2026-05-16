"""
Optimize the topology of a 2D truss.
"""

import math
import itertools

from scipy.spatial import Delaunay
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString, Polygon
from time import time

from inport_optimizers.truss_cross_section import optimize_truss_cross_sections
from opening_saving.saving import save_truss_topology_optimization_excel

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

def optimize_truss(file_path, zone_nodes, load_cases, supports, node_spacing =None, nodes = None, MAX_LENGTH = None):
    """
    Optimizes the topology of a 2D truss.

    :param file_path: File path to the Excel file to save to.
    :param zone_nodes: Nodes that define the polygone that all members must be in.
    :param load_cases: # TODO
    :param supports: Nodes that are supported and in what directions. Shape: (# Supported nodes, 4: [x,y,sx,sy])
    :param node_spacing: How far apart each node column and row should be. Is not needed if nodes are defined.
    :param nodes:  Locations of the nodes of the truss. Shape: (# nodes, 2). Is not needed if node_spacing is defined.
    :param MAX_LENGTH: Maximum length of a member. OPTIONAL
    """

    MAX_LENGTH = 42

    start = time()

    zone = Polygon(zone_nodes)

    if nodes is None:
        nodes, num_columns, num_rows = _create_node_grid(zone, node_spacing)
    else:
        nodes = np.array(nodes)

    structure_cases = _construct_load_cases(load_cases, MAX_LENGTH, nodes, supports)

    members = _create_initial_structure(nodes, zone)

    #_plot_truss(nodes, members, np.ones((len(members))) * 30, 0)

    Vol = []
    A = []
    Q = []
    U = []
    for load_cases in structure_cases:
        if load_cases is not None:
            vol, a, q, u = optimize_truss_cross_sections(nodes, members, load_cases, supports)
            Vol.append(vol)
            A.append(a)
            Q.append(q)
            U.append(u)
        else:
            print("No Load Case Found")

    plt.show()
    # print("finished")
    # print("Time taken: ", time() - start)

    best_index = Vol.index(min(Vol))

    load_cases = structure_cases[best_index]
    areas = A[best_index]
    deflections = U[best_index]
    forces = Q[best_index]
    volume = Vol[best_index]

    save_truss_topology_optimization_excel(file_path, nodes, members, load_cases, supports, areas, deflections, forces, volume)

def _create_node_grid(zone, node_spacing):
    """
    Crates the grid of nodes to be used in the optimization.

    :param zone: A shapely polygon that the truss must be within.
    :param node_spacing: How far apart the nodes should be in the x and y directions.
    :return: Array of nodes, number of columns of nodes, number of rows of nodes.
    """

    min_x, min_y, max_x, max_y = zone.bounds
    dx, dy = max_x - min_x, max_y - min_y
    num_columns, num_rows = dx / node_spacing[0], dy / node_spacing[1]

    xv, yv = np.meshgrid(range(int(num_columns) + 1), range(int(num_rows) + 1))

    pts = [Point(xv.flat[i], yv.flat[i]) for i in range(xv.size)]

    # only add nodes that are inside the zone
    nodes = np.array([[pt.x, pt.y] for pt in pts if zone.intersects(pt)])

    return nodes, num_columns, num_rows

def _create_initial_structure(nodes, zone):
    """
    Creates the list of members.

    :param nodes: Locations of the nodes of the truss. Shape: (# nodes, 2). Is not needed if node_spacing is defined.
    :param zone: A shapely polygon that the truss must be within.
    :return: members: Indices of the nodes that the members run between. Shape: (# members, 3: [node i, node j, length])
    """

    convex = True if zone.convex_hull.area == zone.area else False

    members = []

    tri = Delaunay(nodes)
    for simplex in tri.simplices:
        for i in range(3):
            line = tuple(sorted((simplex[i], simplex[(i + 1) % 3])))
            dx, dy = abs(nodes[line[0]][0] - nodes[line[1]][0]), abs(nodes[line[0]][1] - nodes[line[1]][1])
            seg = [] if convex else LineString([nodes[line[0]], nodes[line[1]]])

            if convex or zone.contains(seg) or zone.boundary.contains(seg):
                members.append([line[0], line[1], np.sqrt(dx ** 2 + dy ** 2)])

    members = np.array(members)
    A1 = np.zeros((len(nodes), len(nodes)))

    A1[members[:, 0].astype(int) , members[:, 1].astype(int) ] = 1
    for k in range(5):
        A1 = A1 @ A1
        np.fill_diagonal(A1, 0)
        A1[A1 > 0] = 1
        A1[A1 != 1] = 0

        A1 = A1 + A1.T
        A1[A1 > 0] = 1

    members = []
    for i in range(len(nodes)):
        for j in range(len(nodes)):
            if i != j and A1[i,j] == 1:
                dx, dy = abs(nodes[j][0] - nodes[i][0]), abs(nodes[j][1] - nodes[i][1])
                seg = [] if convex else LineString([nodes[i], nodes[j]])
                if convex or zone.contains(seg) or zone.boundary.contains(seg):
                    if np.sqrt(dx ** 2 + dy ** 2) <= 42:
                        members.append([i, j, np.sqrt(dx ** 2 + dy ** 2)])

    members = np.array(members)

    return members

def _construct_load_cases(load_cases, max_length, nodes, supports):
    """
    # TODO

    :param load_cases: # TODO
    :param max_length: Maximum length of a member.
    :param nodes: Locations of the nodes of the truss. Shape: (# nodes, 2).
    :param supports: Nodes that are supported and in what directions. Shape: (# Supported nodes, 4: [x,y,sx,sy])
    :return: #TODO
    """

    nodes = np.array(nodes)

    #TODO ONLY FIST LOAD CASE CURRENTLY WORKS
    loads = load_cases[0]

    disLoads = []
    for i, load in enumerate(loads):
        if load[4] is not None and not math.isnan(float(load[4])):
            disLoads.append([load[4], load[5], load[6], load[7], load[8], load[9]])

    # disLoads: list of loads, loads: [x1,y1, x2,y2, fx,fy]
    #print("disLoads: ", disLoads)

    possibleCombinationsOverall = []
    for i, load in enumerate(disLoads):
        p1, p2 = [load[0], load[1]], [load[2], load[3]]

        nodesOnSegment = _get_nodes_on_segment(nodes, p1, p2)

        possibleCombinationsPerDisLoad = []
        n = len(nodesOnSegment)

        # Backtracking to find subsets of intervals that cover the range
        def backtrack(start_idx, current_combination):
            if start_idx == n - 1:
                possibleCombinationsPerDisLoad.append(list(current_combination))
                return

            for end_idx in range(start_idx + 1, n):
                if np.linalg.norm(nodesOnSegment[start_idx] - nodesOnSegment[end_idx]) <= max_length:
                    interval = (nodesOnSegment[start_idx], nodesOnSegment[end_idx], load[4], load[5])
                    backtrack(end_idx, current_combination + [interval])

        backtrack(0, [])

        # possible combinations: list of structure layouts, Structure layouts: list of loads, loads: [[x1,y1], [x2,y2], fx,fy]
        #print("possibleCombinationsPerDisLoad: ", possibleCombinationsPerDisLoad)
        possibleCombinationsOverall.append(possibleCombinationsPerDisLoad)

    # combining possible structure layouts for each distributed load
    possibleCombinationsOverall = list(itertools.product(*possibleCombinationsOverall))
    for i, possibleCombinationIn in enumerate(possibleCombinationsOverall):
        p = []
        for new_loads in possibleCombinationIn:
            for new_load in new_loads:
                p.append(new_load)
        possibleCombinationsOverall[i] = p

    # possibleCombinationsOverall: list of structure layouts, Structure layouts: list of loads, loads: [[x1,y1], [x2,y2], fx,fy]
    #print("possibleCombinationsOverall: ", possibleCombinationsOverall)

    # converting distributed loads into point loads
    possibleStructures = []
    for i, disLoads in enumerate(possibleCombinationsOverall):
        new_loads = []
        for j, disLoad in enumerate(disLoads):
            length = math.dist(disLoad[0], disLoad[1])
            Rx, Ry = disLoad[2] * length / 2, disLoad[3] * length / 2
            new_loads.append([float(disLoad[0][0]), float(disLoad[0][1]), Rx, Ry])
            new_loads.append([float(disLoad[1][0]), float(disLoad[1][1]), Rx, Ry])
        possibleStructures.append(new_loads)

    #print("possibleStructures 1: ", possibleStructures)

    # adding point loads
    for i, new_loads in enumerate(possibleStructures):
        for load in loads:
            if load[0] is not None and not math.isnan(float(load[0])):
                possibleStructures[i] = possibleStructures[i] + [[load[0], load[1], load[2], load[3]]]

    #print("possibleStructures 2: ", possibleStructures)

    for i, loads in enumerate(possibleStructures):
        new_loads = []
        for load in loads:
            loadFound = False
            for new_load in new_loads:
                if float(load[0]) == float(new_load[0]) and float(load[1]) == float(new_load[1]):
                    new_load[2] += load[2]
                    new_load[3] += load[3]
                    loadFound = True
                    break

            if not loadFound:
                for support in supports:
                    if float(support[0]) == float(load[0]) and float(support[1]) == float(load[1]):
                        loadFound = True

            if not loadFound:
                new_loads.append(load)
        possibleStructures[i] = [new_loads] # TODO brackets added as only one load case is possible currently

    #print("possibleStructures 3: ", possibleStructures)

    return possibleStructures


def _get_nodes_on_segment(nodes, n1, n2):
    """
    Gets what nodes are on the segment defined by n1 and n2.

    :param nodes: Locations of the nodes of the truss. Shape: (# nodes, 2).
    :param n1: The starting node of the segment.
    :param n2: The ending node of the segment.
    :return: Ordered ndarray of nodes on the segment. Shape: (# nodes, 2).
    """

    n1 = np.array(n1)
    n2 = np.array(n2)

    nodes_on_segment = nodes[_is_on_line_segment(nodes, n1, n2)]
    nodes_on_segment = _sort_nodes_along_line(nodes_on_segment, n1, n2)
    nodes_on_segment = np.unique(nodes_on_segment, axis=0)

    return nodes_on_segment


def _is_on_line_segment(nodes, n1, n2, tolerance=1e-6):
    """
    Checks which points from a list lie on the node segment defined by w1 and n2.

    :param nodes: Locations of the nodes of the truss. Shape: (# nodes, 2).
    :param n1: The starting node of the segment.
    :param n2: The ending node of the segment.
    :param tolerance: A tolerance for floating-point comparisons.
    :return: A boolean array indicating if each point is on the segment.
    """

    n1 = np.array(n1)
    n2 = np.array(n2)

    vec_p1_to_points = nodes - n1
    vec_segment = n2 - n1

    # check if the nodes are collinear (cross product is near zero)
    cross_products = np.cross(vec_segment, vec_p1_to_points)
    is_collinear = np.abs(cross_products) < tolerance

    # Check if the nodes are between w1 and n2: dot product > 0 and < length of segment squared
    dot_products = np.sum(vec_p1_to_points * vec_segment, axis=1)
    segment_len_sq = np.sum(vec_segment * vec_segment)
    is_within_bounds = (dot_products >= -tolerance) & (dot_products <= segment_len_sq + tolerance)

    return is_collinear & is_within_bounds

def _sort_nodes_along_line(nodes, start_point, end_point):
    """
    Sorts a list of nodes based on their position along the node
    defined by start_point and end_point.

    :param nodes: Locations of the nodes of the truss. Shape: (# nodes, 2).
    :param start_point: The starting node of the segment.
    :param end_point: The ending node of the segment.
    :return: The sorted nodes. Shape: (# nodes, 2)
    """

    direction_vector = end_point - start_point
    unit_direction = direction_vector / np.linalg.norm(direction_vector)

    # calculate projection of each node relative to the start point
    vectors_to_nodes = nodes - start_point
    projections = np.dot(vectors_to_nodes, unit_direction)

    # sort the nodes based on these projection values
    sorted_indices = np.argsort(projections)
    sorted_nodes = nodes[sorted_indices]

    return sorted_nodes



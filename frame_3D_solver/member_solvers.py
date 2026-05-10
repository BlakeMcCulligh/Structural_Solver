"""
Holds functions to find forces and deflections within a member.
"""

import math
import numpy as np

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

def segment_member(members: np.ndarray, members_L: np.ndarray, members_cross_section_props:  np.ndarray,
                   materials: np.ndarray, point_loads: list, dist_loads: list, f_array: np.ndarray,
                   fer_array:  np.ndarray, d_array:  np.ndarray, num_m: int, num_c: int):
    """
    Segments all the members into sections where the loading of the member changes.

    :param members: ndarray. [i_node, j_node, material_index, set_cross_section_props]
    :param members_L: ndarray. length of each member.
    :param members_cross_section_props: ndarray. [A, Iy, Iz, J]
    :param materials: list. [E, G, nu, rho, fy]
    :param point_loads: list. Point loads applied to each member.
                        shape: (# members, # Casses, # loads in case on member, 7: (x, Px, Py, Pz, Mx, My, Mz))
    :param dist_loads: list. Distributed loads applied to each member.
                       shape: (# members, # Casses, varys,
                       # loads in case on member: (x1, x2, wx1, wx2, wy1, wy2, wz1, wz2))
    :param f_array: ndarray. Array of the local forces acting at the ends of the members.
                    shape: (# casses, # members, 12)
    :param fer_array: ndarray. local fixed end reaction vector for the members. shape: (# members, # Casses, 12, 1)
    :param d_array: ndarray. Array of the local end deflections for each member. shape: (# casses, # members, 12)
    :param num_m: Number of members.
    :param num_c: Number of load casses.
    :return:
        seg: list. Holds basic data to do with the segment.
             shape: (# Members, # Casses, # segments: varys, 5: [seg_x1, seg_x2, seg_EIz, seg_EIy, seg_EA])
        seg_internal_loads: list. Holds the end reaction forces of each segment.
                            shape: (# Members, # Casses, # Segments: varys, 3, 2: (P, M))
        seg_dist_loads: list. Holds the distributed loads on the segment.
                        shape: (# Members, # Casses, # segments: varys, 3, 2: (w1, w2))
        seg_thata: list. Holds the slope at the start of the segment.
                   shape: (# Members, # Casses, # Segments: varys, 3)
        seg_delta: list. Holds the deflection at the start of the segment.
                   shape: (# Members, # Casses, # Segments: varys, 3)
    """

    seg = []
    seg_internal_loads = []
    seg_dist_loads = []
    seg_thata = []
    seg_delta = []

    for m_index in range(num_m):

        seg_sub = []
        seg_sub_internal_loads = []
        seg_sub_dist_loads = []
        seg_sub_thata = []
        seg_sub_delta = []

        for c_index in range(num_c):

            # Create a list of discontinuity locations
            disconts = [0, float(members_L[m_index])]
            for load in point_loads[m_index][c_index]:
                disconts.append(load[0])
            for load in dist_loads[m_index][c_index]:
                disconts.append(load[0])  # Distributed load start locations
                disconts.append(load[1])  # Distributed load end locations
            disconts = sorted(set(disconts)) # Sort the list and eliminate duplicate values

            E = materials[members[m_index,2],0]
            A = members_cross_section_props[m_index,0]
            Iy = members_cross_section_props[m_index,1]
            Iz = members_cross_section_props[m_index,2]
            J = members_cross_section_props[m_index,3]
            L = members_L[m_index]

            seg_x1 = []
            seg_x2 = []
            seg_EIz = []
            seg_EIy = []
            seg_EA = []

            for index in range(len(disconts)-1):

                seg_x1.append(disconts[index])
                seg_x2.append(disconts[index+1])
                seg_EIz.append(E * Iz)
                seg_EIy.append(E * Iy)
                seg_EA.append(E * A)

            d = d_array[c_index, m_index]
            f = f_array[c_index, m_index]
            fer = fer_array[m_index, c_index]

            m1z = f[5,0]
            m2z = f[11,0]
            m1y = -f[4,0]
            m2y = -f[10,0]

            fem1z = fer[5,0]
            fem2z = fer[11,0]
            fem1y = -fer[4,0]
            fem2y = -fer[10,0]

            delta1y = d[1]
            delta2y = d[7]
            delta1z = d[2]
            delta2z = d[8]

            seg_deltaY1, seg_deltaZ1 = [delta1z], [delta1y]

            seg_thetaY1 = [-1 / 3 * ((m1y - fem1y) * L / (E * Iy) - (m2y - fem2y) * L
                                     / (2 * E * Iy) + 3 * (delta2z - delta1z) / L)]
            seg_thetaZ1 = [ 1 / 3 * ((m1z - fem1z) * L / (E * Iz) - (m2z - fem2z) * L
                                     / (2 * E * Iz) + 3 * (delta2y - delta1y) / L)]

            seg_deltaX1 = [d[0]]
            seg_L = []
            seg_dist_load_x, seg_dist_load_y, seg_dist_load_z = [], [], []
            seg_internal_load_x, seg_internal_load_y, seg_internal_load_z = [], [], [] # P1, X_T1 # V1, M1 # V1, M1

            for i in range(len(seg_x1)):
                x = seg_x1[i]
                seg_L.append(seg_x2[i] - seg_x1[i])

                seg_internal_load_x.append([f[0, 0], f[3, 0]])  # P1, X_T1
                seg_internal_load_y.append([f[2, 0], f[4, 0] + f[2, 0] * x])  # V1, M1
                seg_internal_load_z.append([f[1, 0], f[5, 0] - f[1, 0] * x])  # V1, M1

                seg_dist_load_z.append([0,0]) # w1, w2
                seg_dist_load_y.append([0,0]) # w1, w2
                seg_dist_load_x.append([0,0]) # w1, n2

                if i > 0:
                    seg_thetaZ1.append(seg_thetaZ1[i-1]
                                       - (-seg_internal_load_z[i-1][0] * seg_L[i-1] **2 / 2
                                            - seg_dist_load_z[i-1][0] * seg_L[i-1] ** 3/6
                                            + seg_L[i-1] * seg_internal_load_z[i-1][1]
                                            + seg_L[i-1] ** 4* (seg_dist_load_z[i-1][0] - seg_dist_load_z[i-1][1])
                                               /(24 * seg_L[i-1])) / seg_EIz[i-1])

                    seg_deltaZ1.append((seg_deltaZ1[i-1] + seg_thetaZ1[i-1] * seg_L[i-1]
                                        + seg_internal_load_z[i-1][0]*seg_L[i-1]**3/(6*seg_EIz[i-1])
                                        + seg_dist_load_z[i-1][0]*x**4/(24*seg_EIz[i-1])
                                        + seg_L[i-1]**2*(-seg_internal_load_z[i-1][1])/(2*seg_EIz[i-1])
                                        + x**5*(-seg_dist_load_z[i-1][0] + seg_dist_load_z[i-1][1])
                                            /(120*seg_EIz[i-1]*seg_L[i-1])))

                    seg_thetaY1.append(seg_thetaY1[i-1]
                                       + (-seg_internal_load_y[i-1][0]*seg_L[i-1]**2/2
                                          - seg_dist_load_y[i-1][0]*seg_L[i-1]**3/6
                                          + seg_L[i-1]*(-seg_internal_load_y[i-1][1])
                                          + seg_L[i-1]**4*(seg_dist_load_y[i-1][0] - seg_dist_load_y[i-1][1])
                                             /(24*seg_L[i-1]))/seg_EIy[i-1])

                    seg_deltaY1.append((seg_deltaY1[i-1] - seg_thetaY1[i-1]*seg_L[i-1]
                                        + seg_internal_load_y[i-1][0]*seg_L[i-1]**3/(6*seg_EIy[i-1])
                                        + seg_dist_load_y[i-1][0]*seg_L[i-1]**4/(24*seg_EIy[i-1])
                                        - seg_L[i-1]**2*(-seg_dist_load_y[i-1][1])/(2*seg_EIy[i-1])
                                        - seg_L[i-1]**5*(seg_dist_load_y[i-1][0] - seg_dist_load_y[i-1][1])
                                            /(120*seg_EIy[i-1]*seg_L[i-1])))

                    seg_deltaX1.append(seg_deltaX1[i-1]
                                       - 1/seg_EA[i-1]*(seg_internal_load_x[i-1][0]*seg_L[i-1]
                                                        + seg_dist_load_x[i-1][0]*seg_L[i-1]**2/2
                                                        + (seg_dist_load_x[i-1][1]
                                                           - seg_dist_load_x[i-1][0])*seg_L[i-1]**3/(6*seg_L[i-1])))

                for point_load in point_loads[m_index][c_index]:

                    if round(point_load[0], 10) <= round(x, 10):

                        seg_internal_load_x[i][0] += point_load[1]
                        seg_internal_load_y[i][0] += point_load[2]
                        seg_internal_load_y[i][1] -= point_load[2] * (x - point_load[0])
                        seg_internal_load_z[i][0] += point_load[3]
                        seg_internal_load_z[i][1] -= point_load[3] * (x - point_load[0])
                        seg_internal_load_x[i][1] += point_load[4]
                        seg_internal_load_y[i][1] += point_load[5]
                        seg_internal_load_z[i][1] += point_load[6]

                for dist_load in dist_loads[m_index][c_index]:

                    if round(dist_load[0], 10) <= round(x, 10):

                        if round(dist_load[1], 10) > round(x, 10):

                            #X Direction
                            seg_dist_load_x[i][0] += ((dist_load[3] - dist_load[2]) / (dist_load[1] - dist_load[0])
                                                      * (x - dist_load[0]) + dist_load[2])

                            seg_dist_load_x[i][1] += ((dist_load[3] - dist_load[2]) / (dist_load[1] - dist_load[0])
                                                      * (seg_x2[i] - dist_load[0]) + dist_load[2])

                            seg_internal_load_x[i][0] += ((dist_load[2]
                                                          + (dist_load[2] + (dist_load[3] - dist_load[2])
                                                             / (dist_load[1] - dist_load[0]) * (x - dist_load[0])))
                                                          / 2 * (x - dist_load[0]))

                            # Y Direction
                            seg_dist_load_z[i][0] += ((dist_load[5] - dist_load[4]) / (dist_load[1] - dist_load[0])
                                                      * (x - dist_load[0]) + dist_load[4])

                            seg_dist_load_z[i][1] += ((dist_load[5] - dist_load[4]) / (dist_load[1] - dist_load[0])
                                                      * (seg_x2[i] - dist_load[0]) + dist_load[4])

                            w2_int = (dist_load[4] + (dist_load[5] - dist_load[4]) / (dist_load[1] - dist_load[0])
                                      * (x - dist_load[0]))

                            seg_internal_load_z[i][0] += (dist_load[4] + w2_int) / 2 * (x - dist_load[0])

                            seg_internal_load_z[i][1] -= ((dist_load[0] - x)
                                                          * (2 * dist_load[4] * dist_load[0] - 3 * dist_load[4] * x
                                                             + dist_load[4] * x + w2_int * dist_load[0]
                                                             - 3 * w2_int * x + 2 * w2_int * x) / 6)

                            # Z Direction
                            seg_dist_load_y[i][0] += ((dist_load[7] - dist_load[6]) / (dist_load[1] - dist_load[0])
                                                      * (x - dist_load[0]) + dist_load[6])

                            seg_dist_load_y[i][1] += ((dist_load[7] - dist_load[6]) / (dist_load[1] - dist_load[0])
                                                      * (seg_x2[i] - dist_load[0]) + dist_load[6])

                            w2_int = (dist_load[3][6] + (dist_load[7] - dist_load[6]) / (dist_load[1] - dist_load[0])
                                      * (x - dist_load[0]))

                            seg_internal_load_y[i][0] += (dist_load[6] + w2_int) / 2 * (x - dist_load[0])

                            seg_internal_load_y[i][1] += ((dist_load[0] - x)
                                                          * (2 * dist_load[6] * dist_load[0] - 3 * dist_load[6] * x
                                                             + dist_load[6] * x + w2_int * dist_load[0]
                                                             - 3 * w2_int * x + 2 * w2_int * x) / 6)

                        else:

                            #X Direction
                            seg_internal_load_x[i][0] += ((dist_load[2] + dist_load[3])
                                                          / 2 * (dist_load[1] - dist_load[0]))

                            # Y direction
                            seg_internal_load_z[i][0] += ((dist_load[4] + dist_load[5])
                                                          / 2 * (dist_load[1] - dist_load[0]))

                            seg_internal_load_z[i][1] -= ((dist_load[0] - dist_load[1])
                                                          * (2 * dist_load[4] * dist_load[0] - 3 * dist_load[4] * x
                                                             + dist_load[4] * dist_load[1]
                                                             + dist_load[5] * dist_load[0]
                                                             - 3 * dist_load[5] * x
                                                             + 2 * dist_load[5] * dist_load[1]) / 6)

                            # z Direction
                            seg_internal_load_y[i][0] += ((dist_load[6] + dist_load[7])
                                                          / 2 * (dist_load[1] - dist_load[0]))

                            seg_internal_load_y[i][1] += ((dist_load[0] - dist_load[1])
                                                          * (2 * dist_load[6] * dist_load[0] - 3 * dist_load[6] * x
                                                             + dist_load[6] * dist_load[1]
                                                             + dist_load[7] * dist_load[0] - 3 * dist_load[7] * x
                                                             + 2 * dist_load[7] * dist_load[1]) / 6)

            seg_sub.append([seg_x1, seg_x2, seg_EIz, seg_EIy, seg_EA])
            seg_sub_internal_loads.append([seg_internal_load_x, seg_internal_load_y, seg_internal_load_z]) # [P1, T1], [Vy1, My1], [Vz1, Mz1]
            seg_sub_dist_loads.append([seg_dist_load_x, seg_dist_load_y, seg_dist_load_z]) # [wx1, wx2], [wy1, wy2], [wz1, wz2]
            seg_sub_thata.append([seg_thetaY1, seg_thetaZ1])
            seg_sub_delta.append([seg_deltaX1, seg_deltaY1, seg_deltaZ1])

        seg.append(seg_sub)
        seg_internal_loads.append(seg_sub_internal_loads)
        seg_dist_loads.append(seg_sub_dist_loads)
        seg_thata.append(seg_sub_thata)
        seg_delta.append(seg_sub_delta)

    return seg, seg_internal_loads, seg_dist_loads, seg_thata, seg_delta

def _extrem_finder(casses: list, m_index: int, seg: list, seg_internal_loads: list, seg_dist_loads: list,
                   abs_function: callable, abs_function_direction: str, poi_function: callable, direction: int,
                   sign: int | None = None, combo_indexs: list | None = None):
    """
    Finds the most extrem value and governing load case for the given member and direction.

    :param casses: list. Load Case Indexes.
    :param m_index: int. Index of the member to find values for.
    :param seg: list. Holds basic data to do with the segment.
                shape: (# Members, # Casses, # segments: varys, 5: [seg_x1, seg_x2, seg_EIz, seg_EIy, seg_EA])
    :param seg_internal_loads: list. Holds the end reaction forces of each segment.
                               shape: (# Members, # Casses, # Segments: varys, 3, 2: [P, M])
    :param seg_dist_loads: list. Holds the distributed loads on the segment.
                           shape: (# Members, # Casses, # segments: varys, 3, 2: (w1, w2))
    :param abs_function: callable. Finds the min or max of a list.
    :param abs_function_direction: String. Is the min or the max being found. ("min", "max")
    :param poi_function: callable. Returns a list of values that could be the maximum or minimum.
    :param direction: int. Direction to find the extrem for. (0=x, 1=y, 2=z)
    :param sign: int. Some function require subtractions insted of additions for cirten directions
                 so passing a -1 converts them to subtraction.
    :param combo_indexs: list. Indeces of combos to be checked.
    :return: Extreme Value, Governing Load Case
    """

    if combo_indexs is None: combo_indexs = casses

    global_, governing_combo = None, None

    seg, seg_internal_loads, seg_dist_loads = seg[m_index], seg_internal_loads[m_index], seg_dist_loads[m_index]

    for combo_index in combo_indexs:

        abs_ = []

        L = np.array(seg[combo_index][1]) - np.array(seg[combo_index][0])

        for seg_index in range(len(seg)):

            abs_.append(abs_function(
                poi_function(seg_dist_loads[combo_index][direction][seg_index][0],
                             seg_dist_loads[combo_index][direction][seg_index][1],
                             seg_internal_loads[combo_index][direction][seg_index][0],
                             seg_internal_loads[combo_index][direction][seg_index][1], L[combo_index], sign)))

        abs_ = abs_function(abs_)

        if abs_function_direction == "max":
            if global_ is None or abs_ > global_: global_, governing_combo = abs_, combo_index

        elif abs_function_direction == "min":
            if global_ is None or abs_ < global_: global_, governing_combo = abs_, combo_index

    return global_, governing_combo


""" --------------- SHEAR --------------- """


def shear(x: float, m_index: int, combo_index: int, members_L: np.ndarray,
          seg: list, seg_internal_loads: list, seg_dist_loads: list):
    """
    Findes the shear at location x.

    :param x: float. Loaction to find the shear at.
    :param m_index: int. Index of the member to find values for.
    :param combo_index: int. Load case to find the shear for.
    :param members_L: ndarray. length of each member.
    :param seg: list. Holds basic data to do with the segment.
                shape: (# Members, # Casses, # segments: varys, 5: [seg_x1, seg_x2, seg_EIz, seg_EIy, seg_EA])
    :param seg_internal_loads: list. Holds the end reaction forces of each segment.
                               shape: (# Members, # Casses, # Segments: varys, 3, 2: [P, M])
    :param seg_dist_loads: list. Holds the distributed loads on the segment.
                           shape: (# Members, # Casses, # segments: varys, 3, 2: (w1, w2))
    :return: Shear Y, Shear X
    """

    seg = seg[m_index][combo_index]
    seg_internal_loads = seg_internal_loads[m_index][combo_index]
    seg_dist_loads = seg_dist_loads[m_index][combo_index]

    for i in range(len(seg[0])):

        if round(seg[0][i], 10) <= round(x, 10) < round(seg[1][i], 10):
            VY = (seg_internal_loads[1][i][0] + seg_dist_loads[1][i][0] * (x - seg[0][i])
                  + (x - seg[0][i]) ** 2 * (-seg_dist_loads[1][i][0] + seg_dist_loads[1][i][1])
                  / (2 * (seg[1][i] - seg[0][i])))

            VZ = (seg_internal_loads[2][i][0] + seg_dist_loads[2][i][0] * (x - seg[0][i])
                  + (x - seg[0][i]) ** 2 * (-seg_dist_loads[2][i][0] + seg_dist_loads[2][i][1])
                  / (2 * (seg[1][i] - seg[0][i])))

            return VY, VZ

    if math.isclose(x, members_L[m_index]):

        last_index = len(seg[0]) - 1

        VY = (seg_internal_loads[1][last_index][0] + seg_dist_loads[1][last_index][0] * (x - seg[0][last_index])
              + (x - seg[0][last_index]) ** 2 * (-seg_dist_loads[1][last_index][0] + seg_dist_loads[1][last_index][1])
              / (2 * (seg[1][last_index] - seg[0][last_index])))

        VZ = (seg_internal_loads[2][last_index][0] + seg_dist_loads[2][last_index][0] * (x - seg[0][last_index])
              + (x - seg[0][last_index]) ** 2 * (-seg_dist_loads[2][last_index][0] + seg_dist_loads[2][last_index][1])
              / (2 * (seg[1][last_index] - seg[0][last_index])))

        return VY, VZ

    return 0, 0

def _shear_loc_of_intrest(w1: float, w2: float, L: float):
    """
    Returns the shear location of intrest for the segment.

    :param w1: float. Distributed load at the start of the segment.
    :param w2: float. Distributed load at the end of the segment.
    :param L: float. Length of the segment.
    :return: x: float. Posible location of extreme shear.
    """

    if w1 - w2 == 0: x1 = 0
    else: x1 = w1 * L / (w1 - w2)
    if round(x1, 10) < 0 or round(x1, 10) > round(L, 10): x1 = 0
    return x1

def _seg_V_POI(w1: float, w2: float, V1: float, M_1: float, L: float, sign):
    """
    Finds the shear at all the locations where it could be a maximum or minimum.

    :param w1: float. Distributed load at the start of the segment.
    :param w2: float. Distributed load at the end of the segment.
    :param V1: float. Shear at the start of the segment.
    :param M_1: float. Moment at the start of the segment.
    :param L: float. Length of the segment.
    :param sign: NOT USED
    :return: list of posible extreme shears.
    """

    x = _shear_loc_of_intrest(w1, w2, L)
    shear1 = V1 + w1 * x + x ** 2 * (-w1 + w2) / (2 * L)
    shear2 = V1
    shear3 = V1 + w1 * L + L ** 2 * (-w1 + w2) / (2 * L)
    return [shear1, shear2, shear3]

def max_shear(casses, m_index, seg, seg_internal_loads, seg_dist_loads, combo_indexs = None):
    """
    Finds the maximum shear for the member and its governing load case.

    :param casses: list. Load Case Indexes.
    :param m_index: int. Index of the member to find values for.
    :param seg: list. Holds basic data to do with the segment.
                shape: (# Members, # Casses, # segments: varys, 5: [seg_x1, seg_x2, seg_EIz, seg_EIy, seg_EA])
    :param seg_internal_loads: list. Holds the end reaction forces of each segment.
                               shape: (# Members, # Casses, # Segments: varys, 3, 2: [P, M])
    :param seg_dist_loads: list. Holds the distributed loads on the segment.
                           shape: (# Members, # Casses, # segments: varys, 3, 2: (w1, w2))
    :param combo_indexs: int. Load case to find the shear for. OPTIONAL
    :return: Shear Y, Shear Z, Governing Case Y, Governing Case Z
    """

    Y, yC = _extrem_finder(casses, m_index, seg, seg_internal_loads, seg_dist_loads, max, "max",
                           _seg_V_POI, 1, combo_indexs=combo_indexs)
    Z, zC = _extrem_finder(casses, m_index, seg, seg_internal_loads, seg_dist_loads, max, "max",
                           _seg_V_POI, 2, combo_indexs=combo_indexs)

    return Y, Z, yC, zC

def min_shear(casses, m_index, seg, seg_internal_loads, seg_dist_loads, combo_indexs = None):
    """
    Finds the minimum shear for the member and its governing load case.

    :param casses: list. Load Case Indexes.
    :param m_index: int. Index of the member to find values for.
    :param seg: list. Holds basic data to do with the segment.
                shape: (# Members, # Casses, # segments: varys, 5: [seg_x1, seg_x2, seg_EIz, seg_EIy, seg_EA])
    :param seg_internal_loads: list. Holds the end reaction forces of each segment.
                                shape: (# Members, # Casses, # Segments: varys, 3, 2: [P, M])
    :param seg_dist_loads: list. Holds the distributed loads on the segment.
                            shape: (# Members, # Casses, # segments: varys, 3, 2: (w1, w2))
    :param combo_indexs: int. Load case to find the shear for. OPTIONAL
    :return: Shear Y, Shear Z, Governing Case Y, Governing Case Z
    """

    Y, yC = _extrem_finder(casses, m_index, seg, seg_internal_loads, seg_dist_loads, min, "min",
                           _seg_V_POI, 1, combo_indexs=combo_indexs)
    Z, zC = _extrem_finder(casses, m_index, seg, seg_internal_loads, seg_dist_loads, min, "min",
                           _seg_V_POI, 2, combo_indexs=combo_indexs)

    return Y, Z, yC, zC


""" --------------- MOMENT --------------- """


def _moment_calc(M1: float, V1: float, w1: float, w2: float, L: float, x: float, sign: int):
    """
    Calculates the moment at location x.

    :param M1: float. Moment at the start of the segment.
    :param V1: float. Shear at the start of the segment.
    :param w1: float. Distributed load at the start of the segment.
    :param w2: float. Distributed load at the end of the segment.
    :param L: float. Length of the segment.
    :param x: float. Location to get the moment at.
    :param sign: int. Used for converting values to negative.
    :return: Moment at location x.
    """

    return sign*M1 - V1*x - w1*x**2/2 - x**3*(-w1 + w2)/(6*L)

def moment(x: float, m_index: int, combo_index: int, members_L: np.ndarray,
           seg: list, seg_internal_loads: list, seg_dist_loads: list):
    """
    Calculates the moment at location x.

    :param x: float. Location to get the moment at.
    :param m_index: int. Index of the member.
    :param combo_index: int. Index of the load case.
    :param members_L: ndarray. length of each member.
    :param seg: list. Holds basic data to do with the segment.
                shape: (# Members, # Casses, # segments: varys, 5: [seg_x1, seg_x2, seg_EIz, seg_EIy, seg_EA])
    :param seg_internal_loads: list. Holds the end reaction forces of each segment.
                               shape: (# Members, # Casses, # Segments: varys, 3, 2: [P, M])
    :param seg_dist_loads: list. Holds the distributed loads on the segment.
                           shape: (# Members, # Casses, # segments: varys, 3, 2: (w1, w2))
    :return: Moment Y, Moment Z
    """

    seg = seg[m_index][combo_index]
    seg_internal_loads = seg_internal_loads[m_index][combo_index]
    seg_dist_loads = seg_dist_loads[m_index][combo_index]

    for i in range(len(seg[0])):
        if round(seg[0][i], 10) <= round(x, 10) < round(seg[1][i], 10):

            MZ = _moment_calc(seg_internal_loads[2][i][1], seg_internal_loads[2][i][0], seg_dist_loads[2][i][0],
                              seg_dist_loads[2][i][1], (seg[1][i] - seg[0][i]), (x - seg[0][i]), 1)

            MY = _moment_calc(seg_internal_loads[1][i][1], seg_internal_loads[1][i][0], seg_dist_loads[1][i][0],
                              seg_dist_loads[1][i][1], (seg[1][i] - seg[0][i]), (x - seg[0][i]), -1)
            return MZ, MY

    if math.isclose(x, members_L[m_index]):
        i = len(seg[0]) - 1

        MZ = _moment_calc(seg_internal_loads[2][i][1], seg_internal_loads[2][i][0], seg_dist_loads[2][i][0],
                          seg_dist_loads[2][i][1], (seg[1][i] - seg[0][i]), (x - seg[0][i]), 1)

        MY = _moment_calc(seg_internal_loads[1][i][1], seg_internal_loads[1][i][0], seg_dist_loads[1][i][0],
                          seg_dist_loads[1][i][1], (seg[1][i] - seg[0][i]), (x - seg[0][i]), -1)
        return MY, MZ
    return 0, 0

def _moment_loc_of_interist(w1: float, w2: float, V1: float, L: float):
    """
    Finds points where the moment could be at its maximum or minimum.

    :param w1: float. Distributed load at the start of the segment.
    :param w2: float. Distributed load at the end of the segment.
    :param V1: float. Shear at the start of the segment.
    :param L: float. Length of the segment.
    :return: x1, x2
    """

    a = -(w2 - w1) / (2 * L)
    b = -w1
    c = -V1

    if a == 0:
        if b != 0: x1 = -c / b
        else: x1 = 0
        x2 = 0

    elif b ** 2 - 4 * a * c < 0:
        x1, x2 = 0, 0

    else:
        x1 = (-b + (b ** 2 - 4 * a * c) ** 0.5) / (2 * a)
        x2 = (-b - (b ** 2 - 4 * a * c) ** 0.5) / (2 * a)

    if round(x1, 10) < 0 or round(x1, 10) > round(L, 10): x1 = 0
    if round(x2, 10) < 0 or round(x2, 10) > round(L, 10): x2 = 0

    return x1, x2

def _seg_M_POI(w1: float, w2: float, V1: float, M_1: float, L: float, sign: int):
    """
    Finds the moment at all the locations where it could be a maximum or minimum.

    :param w1: float. Distributed load at the start of the segment.
    :param w2: float. Distributed load at the end of the segment.
    :param V1: float. Shear at the start of the segment.
    :param M_1: float. Moment at the start of the segment.
    :param L: float. Length of the segment.
    :param sign: int. Used for converting values to negative.
    :return: list of posible extreme moments.
    """

    x1, x2 = _moment_loc_of_interist(w1, w2, V1, L)
    M1 = _moment_calc(M_1, V1, w1, w2, L, x1, sign)
    M2 = _moment_calc(M_1, V1, w1, w2, L, x2, sign)
    M3 = _moment_calc(M_1, V1, w1, w2, L, 0, sign)
    M4 = _moment_calc(M_1, V1, w1, w2, L, L, sign)
    return [M1, M2, M3, M4]

def max_moment(casses: list, m_index: int, seg: list, seg_internal_loads: list, seg_dist_loads: list,
               combo_indexs: int | None = None):
    """
    Finds the maximum moment for the member and its governing load case.

    :param casses: list. Load Case Indexes.
    :param m_index: int. Index of the member to find values for.
    :param seg: list. Holds basic data to do with the segment.
                shape: (# Members, # Casses, # segments: varys, 5: [seg_x1, seg_x2, seg_EIz, seg_EIy, seg_EA])
    :param seg_internal_loads: list. Holds the end reaction forces of each segment.
                                shape: (# Members, # Casses, # Segments: varys, 3, 2: [P, M])
    :param seg_dist_loads: list. Holds the distributed loads on the segment.
                            shape: (# Members, # Casses, # segments: varys, 3, 2: (w1, w2))
    :param combo_indexs: int. Load case to find the moment for. OPTIONAL
    :return: Moment Y, Moment Z, Governing Case Y, Governing Case Z
    """

    Y, yC = _extrem_finder(casses, m_index, seg, seg_internal_loads, seg_dist_loads, max, "max",
                           _seg_M_POI, 1, sign=-1, combo_indexs= combo_indexs)
    Z, zC = _extrem_finder(casses, m_index, seg, seg_internal_loads, seg_dist_loads, max, "max",
                           _seg_M_POI, 2, sign=1, combo_indexs=combo_indexs)

    return Y, Z, yC, zC

def min_moment(casses: list, m_index: int, seg: list, seg_internal_loads: list, seg_dist_loads: list,
               combo_indexs: int | None = None):
    """
    Finds the minimum moment for the member and its governing load case.

    :param casses: list. Load Case Indexes.
    :param m_index: int. Index of the member to find values for.
    :param seg: list. Holds basic data to do with the segment.
                shape: (# Members, # Casses, # segments: varys, 5: [seg_x1, seg_x2, seg_EIz, seg_EIy, seg_EA])
    :param seg_internal_loads: list. Holds the end reaction forces of each segment.
                                shape: (# Members, # Casses, # Segments: varys, 3, 2: [P, M])
    :param seg_dist_loads: list. Holds the distributed loads on the segment.
                            shape: (# Members, # Casses, # segments: varys, 3, 2: (w1, w2))
    :param combo_indexs: int. Load case to find the moment for. OPTIONAL
    :return: Moment Y, Moment Z, Governing Case Y, Governing Case Z
    """

    Y, yC = _extrem_finder(casses, m_index, seg, seg_internal_loads, seg_dist_loads, min, "min",
                           _seg_M_POI, 1, sign=-1, combo_indexs=combo_indexs)
    Z, zC = _extrem_finder(casses, m_index, seg, seg_internal_loads, seg_dist_loads, min, "min",
                           _seg_M_POI, 2, sign=1, combo_indexs=combo_indexs)

    return Y, Z, yC, zC


""" --------------- TORQUE --------------- """


def torque(x: float, m_index: int, combo_index: int, members_L: np.ndarray, seg: list, seg_internal_loads: list):
    """
    Calculates the torque at location x.

    :param x: float. Location to get the torque at.
    :param m_index: int. Index of the member.
    :param combo_index: int. Index of the load case.
    :param members_L: ndarray. length of each member.
    :param seg: list. Holds basic data to do with the segment.
                shape: (# Members, # Casses, # segments: varys, 5: [seg_x1, seg_x2, seg_EIz, seg_EIy, seg_EA])
    :param seg_internal_loads: list. Holds the end reaction forces of each segment.
                               shape: (# Members, # Casses, # Segments: varys, 3, 2: [P, M])
    :return: Torque
    """

    seg = seg[m_index][combo_index]
    seg_internal_loads = seg_internal_loads[m_index][combo_index]

    for i in range(len(seg[0])):
        if round(seg[0][i], 10) <= round(x, 10) < round(seg[1][i], 10):
            return seg_internal_loads[0][i][1]

    if math.isclose(x, members_L[m_index]):
        i = len(seg[0]) - 1
        return seg_internal_loads[0][i][1]

    return 0

def _seg_T_POI(w1, w2, P1, M1, L, sign):
    """
    Calculates the torque where it can be a maximum or minimum.

    :param w1: Not Used
    :param w2: Not Used
    :param P1: Not Used
    :param M1: Moment at the start of the member.
    :param L: Not Used
    :param sign: Not Used
    :return: List of posible max or min torques
    """

    return [M1]

def max_tourque(casses, m_index, seg, seg_internal_loads, seg_dist_loads, combo_indexs = None):
    """
    Finds the maximum torque for the member and its governing load case.

    :param casses: list. Load Case Indexes.
    :param m_index: int. Index of the member to find values for.
    :param seg: list. Holds basic data to do with the segment.
                shape: (# Members, # Casses, # segments: varys, 5: [seg_x1, seg_x2, seg_EIz, seg_EIy, seg_EA])
    :param seg_internal_loads: list. Holds the end reaction forces of each segment.
                               shape: (# Members, # Casses, # Segments: varys, 3, 2: [P, M])
    :param seg_dist_loads: list. Holds the distributed loads on the segment.
                           shape: (# Members, # Casses, # segments: varys, 3, 2: (w1, w2))
    :param combo_indexs: int. Load case to find the torque for. OPTIONAL
    :return: Torque, Governing Case
    """

    return _extrem_finder(casses, m_index, seg, seg_internal_loads, seg_dist_loads, max, "max",
                          _seg_T_POI, 0, combo_indexs=combo_indexs)

def min_tourque(casses, m_index, seg, seg_internal_loads, seg_dist_loads, combo_indexs = None):
    """
    Finds the minimum torque for the member and its governing load case.

    :param casses: list. Load Case Indexes.
    :param m_index: int. Index of the member to find values for.
    :param seg: list. Holds basic data to do with the segment.
                shape: (# Members, # Casses, # segments: varys, 5: [seg_x1, seg_x2, seg_EIz, seg_EIy, seg_EA])
    :param seg_internal_loads: list. Holds the end reaction forces of each segment.
                                shape: (# Members, # Casses, # Segments: varys, 3, 2: [P, M])
    :param seg_dist_loads: list. Holds the distributed loads on the segment.
                            shape: (# Members, # Casses, # segments: varys, 3, 2: (w1, w2))
    :param combo_indexs: int. Load case to find the torque for. OPTIONAL
    :return: Torque, Governing Case
    """

    return _extrem_finder(casses, m_index, seg, seg_internal_loads, seg_dist_loads, min, "min",
                          _seg_T_POI, 0, combo_indexs=combo_indexs)


""" --------------- AXIAL --------------- """


def _axial_calc(p1: float, p2: float, P1: float, L: float, x: float):
    """
    Calculates the axial force at location x.

    :param p1: float. Distributed axial force at the start of the segment.
    :param p2: float. Distributed axial force at the end of the segment.
    :param P1: float. Axial force at the start of the segment.
    :param L: float. Length of the segment.
    :param x: float. Location to solve for the axial force.
    :return: float. Axial Force
    """
    return P1 + (p2 - p1)/(2*L)*x**2 + p1*x

def axial(x: float, m_index: int, combo_index: int, members_L: np.ndarray,
          seg: list, seg_internal_loads: list, seg_dist_loads: list):
    """
    Calculates the axial force at location x.

    :param x: float. Location to get the axial force at.
    :param m_index: int. Index of the member.
    :param combo_index: int. Index of the load case.
    :param members_L: ndarray. length of each member.
    :param seg: list. Holds basic data to do with the segment.
                shape: (# Members, # Casses, # segments: varys, 5: [seg_x1, seg_x2, seg_EIz, seg_EIy, seg_EA])
    :param seg_internal_loads: list. Holds the end reaction forces of each segment.
                                shape: (# Members, # Casses, # Segments: varys, 3, 2: [P, M])
    :param seg_dist_loads: list. Holds the distributed loads on the segment.
                            shape: (# Members, # Casses, # segments: varys, 3, 2: (w1, w2))
    :return: Axial Force
    """

    seg = seg[m_index][combo_index]
    seg_internal_loads = seg_internal_loads[m_index][combo_index]
    seg_dist_loads = seg_dist_loads[m_index][combo_index]

    for i in range(len(seg[0])):
        if round(seg[0][i], 10) <= round(x, 10) < round(seg[1][i], 10):

            return _axial_calc(seg_dist_loads[0][1][0], seg_dist_loads[0][1][1], seg_internal_loads[0][i][0],
                               (seg[1][i] - seg[0][i]), (x - seg[0][i]))

    if math.isclose(x, members_L[m_index]):
        i = len(seg[0]) - 1

        return _axial_calc(seg_dist_loads[0][1][0], seg_dist_loads[0][1][1], seg_internal_loads[0][i][0],
                           (seg[1][i] - seg[0][i]), (x - seg[0][i]))
    return 0

def _axial_loc_of_interist(p1: float, p2: float, L: float):
    """
    Finds points where the axial force could be at its maximum or minimum.

    :param p1: float. Distributed load at the start of the segment.
    :param p2: float. Distributed load at the end of the segment.
    :param L: float. Length of the segment.
    :return: x
    """

    if p1 - p2 != 0:
        x1 = L * p1 / (p1 - p2)
    else:
        x1 = 0
    if round(x1, 10) < 0 or round(x1, 10) > round(L, 10):
        x1 = 0
    return x1

def _seg_P_POI(w1: float, w2: float, P1: float, M1: float, L: float, sign: int):
    """
    Finds the axial force at all the locations where it could be a maximum or minimum.

    :param w1: float. Distributed load at the start of the segment.
    :param w2: float. Distributed load at the end of the segment.
    :param P1: float. Axial force at the start of the segment.
    :param M1: float. Moment at the start of the segment.
    :param L: float. Length of the segment.
    :param sign: int. Used for converting values to negative.
    :return: list of posible extreme axial forces.
    """

    x1 = _axial_loc_of_interist(w1, w2, L)
    P_1 = _axial_calc(w1, w2, P1, L, x1)
    P_2 = _axial_calc(w1, w2, P1, L, 0)
    P_3 = _axial_calc(w1, w2, P1, L, L)
    return [P_1, P_2, P_3]

def max_axial(casses: list, m_index: int, seg: list, seg_internal_loads: list, seg_dist_loads: list,
              combo_indexs: int | None = None):
    """
    Finds the maximum axial force for the member and its governing load case.

    :param casses: list. Load Case Indexes.
    :param m_index: int. Index of the member to find values for.
    :param seg: list. Holds basic data to do with the segment.
                shape: (# Members, # Casses, # segments: varys, 5: [seg_x1, seg_x2, seg_EIz, seg_EIy, seg_EA])
    :param seg_internal_loads: list. Holds the end reaction forces of each segment.
                                shape: (# Members, # Casses, # Segments: varys, 3, 2: [P, M])
    :param seg_dist_loads: list. Holds the distributed loads on the segment.
                            shape: (# Members, # Casses, # segments: varys, 3, 2: (w1, w2))
    :param combo_indexs: int. Load case to find the axial force for. OPTIONAL
    :return: Axial Force, Governing Case
    """

    return _extrem_finder(casses, m_index, seg, seg_internal_loads, seg_dist_loads, max, "max",
                          _seg_P_POI, 0, combo_indexs)

def min_axial(casses: list, m_index: int, seg: list, seg_internal_loads: list, seg_dist_loads: list,
              combo_indexs: int | None = None):
    """
    Finds the minimum axial force for the member and its governing load case.

    :param casses: list. Load Case Indexes.
    :param m_index: int. Index of the member to find values for.
    :param seg: list. Holds basic data to do with the segment.
                shape: (# Members, # Casses, # segments: varys, 5: [seg_x1, seg_x2, seg_EIz, seg_EIy, seg_EA])
    :param seg_internal_loads: list. Holds the end reaction forces of each segment.
                                shape: (# Members, # Casses, # Segments: varys, 3, 2: [P, M])
    :param seg_dist_loads: list. Holds the distributed loads on the segment.
                            shape: (# Members, # Casses, # segments: varys, 3, 2: (w1, w2))
    :param combo_indexs: int. Load case to find the axial force for. OPTIONAL
    :return: Axial Force, Governing Case
    """

    return _extrem_finder(casses, m_index, seg, seg_internal_loads, seg_dist_loads, min, "min",
                          _seg_P_POI, 0, combo_indexs)


""" --------------- DEFLECTION --------------- """


def _axial_deflection_calc(delta_x1: float, EA: float, P1: float, w1: float, w2: float, L: float, x: float):
    """
    Calculates the axial deflection at location x.

    :param delta_x1: float. X direction deflection at the start of the segment.
    :param EA: float. Youngs Modulus times Area for the segment.
    :param P1: float. Axial force at the start of the segment.
    :param w1: float. Distributed axial force at the start of the segment.
    :param w2: float. Fistributed axial force at the end of the segment.
    :param L: float. Length of the segment.
    :param x: float. Location to find the axial deflection for.
    :return: deflction X
    """

    return delta_x1 - 1/EA*(P1 * x + w1 * x ** 2 / 2 + (w2 - w1) * x ** 3 / (6 * L))

def _deflection_calc(delta1: float, theta1: float, V1: float, EI: float, w1: float, w2: float, M1: float, L: float,
                     sign: int, x: float):
    """
    Calculates the deflection at location x.

    :param delta1: float. Deflection at the start of the segment.
    :param theta1: float. Slope at the start of the segment.
    :param V1: float. Shear at the start of the segment.
    :param EI: float. Youngs Modulus times Moment of Ineta for the segment.
    :param w1: float. Distributed load at the start of the segment.
    :param w2: float. Distributed load at the end of the segment.
    :param M1: float. Moment at the start of the segment.
    :param L: float. Length of the segment.
    :param sign: int. Used for converting values to negative.
    :param x: Float. Location to find the deflection for.
    :return: deflection
    """

    return (delta1 + sign * theta1*x + V1*x**3/(6 * EI) + w1*x**4/(24 * EI) + sign*x**2*(-M1)/(2 * EI)
            + sign*x**5*(sign* (-w1) + sign*w2)/(120 * EI * L))

def deflection(x: float, m_index: int, combo_index: int, members_L: np.ndarray, seg: list, seg_internal_loads: list,
               seg_dist_loads: list, seg_thata: list, seg_delta: list):
    """
    Finds the deflection of the member at location x.

    :param x: float. Location to find the deflection for.
    :param m_index: int. Index of the member to find values for.
    :param combo_index: int. Load case to find the deflection for.
    :param members_L: ndarray. length of each member.
    :param seg: list. Holds basic data to do with the segment.
                shape: (# Members, # Casses, # segments: varys, 5: [seg_x1, seg_x2, seg_EIz, seg_EIy, seg_EA])
    :param seg_internal_loads: list. Holds the end reaction forces of each segment.
                                shape: (# Members, # Casses, # Segments: varys, 3, 2: [P, M])
    :param seg_dist_loads: list. Holds the distributed loads on the segment.
                            shape: (# Members, # Casses, # segments: varys, 3, 2: (w1, w2))
    :param seg_thata: list. Holds the slope at the start of the segment.
               shape: (# Members, # Casses, # Segments: varys, 3)
    :param seg_delta: list. Holds the deflection at the start of the segment.
               shape: (# Members, # Casses, # Segments: varys, 3)
    :return: Deflection X, Deflection Y, Deflection Z
    """

    seg = seg[m_index][combo_index]
    seg_internal_loads = seg_internal_loads[m_index][combo_index]
    seg_dist_loads = seg_dist_loads[m_index][combo_index]
    seg_thata = seg_thata[m_index][combo_index]
    seg_delta = seg_delta[m_index][combo_index]

    for i in range(len(seg[0])):
        if round(seg[0][i], 10) <= round(x, 10) < round(seg[1][i], 10):

            DX = _axial_deflection_calc(seg_delta[0][i], seg[4][i], seg_internal_loads[0][i][0],
                                        seg_dist_loads[0][i][0], seg_dist_loads[0][i][1],
                                        (seg[1][i] - seg[0][i]), (x - seg[0][i]))

            DY = _deflection_calc(seg_delta[1][i], seg_thata[0][i], seg_internal_loads[1][i][0], seg[3][i],
                                  seg_dist_loads[1][i][0], seg_dist_loads[1][i][1], seg_internal_loads[1][i][1],
                                  (seg[1][i] - seg[0][i]), -1, (x - seg[0][i]))

            DZ = _deflection_calc(seg_delta[2][i], seg_thata[1][i], seg_internal_loads[2][i][0], seg[2][i],
                                  seg_dist_loads[2][i][0], seg_dist_loads[2][i][1], seg_internal_loads[2][i][1],
                                  (seg[1][i] - seg[0][i]), 1, (x - seg[0][i]))
            return [DX, DY, DZ]

    if math.isclose(x, members_L[m_index]):
        i = len(seg[0]) - 1

        DX = _axial_deflection_calc(seg_delta[0][i], seg[4][i], seg_internal_loads[0][i][0], seg_dist_loads[0][i][0],
                                    seg_dist_loads[0][i][1], (seg[1][i] - seg[0][i]), (x - seg[0][i]))

        DY = _deflection_calc(seg_delta[1][i], seg_thata[0][i], seg_internal_loads[1][i][0], seg[3][i],
                              seg_dist_loads[1][i][0], seg_dist_loads[1][i][1], seg_internal_loads[1][i][1],
                              (seg[1][i] - seg[0][i]), -1, (x - seg[0][i]))

        DZ = _deflection_calc(seg_delta[2][i], seg_thata[1][i], seg_internal_loads[2][i][0], seg[2][i],
                              seg_dist_loads[2][i][0], seg_dist_loads[2][i][1], seg_internal_loads[2][i][1],
                              (seg[1][i] - seg[0][i]), 1, (x - seg[0][i]))
        return [DX, DY, DZ]

    return [0,0,0]

def _difflection_extrem_finder(casses: list, m_index: int, seg: list, seg_internal_loads: list, seg_dist_loads: list,
                               seg_thata: list, seg_delta: list, abs_function_direction: str,
                               members_L: np.ndarray, combo_indexs: list | None = None):
    """
    Finds the most extrem value and governing load case for the given member and direction.

    :param casses: list. Load Case Indexes.
    :param m_index: int. Index of the member to find values for.
    :param seg: list. Holds basic data to do with the segment.
                shape: (# Members, # Casses, # segments: varys, 5: [seg_x1, seg_x2, seg_EIz, seg_EIy, seg_EA])
    :param seg_internal_loads: list. Holds the end reaction forces of each segment.
                               shape: (# Members, # Casses, # Segments: varys, 3, 2: [P, M])
    :param seg_dist_loads: list. Holds the distributed loads on the segment.
                           shape: (# Members, # Casses, # segments: varys, 3, 2: (w1, w2))
    :param seg_thata: list. Holds the slope at the start of the segment.
                      shape: (# Members, # Casses, # Segments: varys, 3)
    :param seg_delta: list. Holds the deflection at the start of the segment.
                      shape: (# Members, # Casses, # Segments: varys, 3)
    :param abs_function_direction: String. Is the min or the max being found. ("min", "max")
    :param members_L: ndarray. Length of each member.
    :param combo_indexs: list. Indeces of combos to be checked.
    :return: Deflection, Governing Casse
    """

    if combo_indexs is None: combo_indexs = casses

    global_, governing_combo = [None, None, None], [None, None, None]

    L = float(members_L[m_index])

    for comboINDEX in combo_indexs:

        dmax = deflection(0, m_index, comboINDEX, members_L, seg, seg_internal_loads, seg_dist_loads,
                          seg_thata, seg_delta)

        for i in range(100):

            d = deflection(L * i / 99, m_index, comboINDEX, members_L, seg, seg_internal_loads,
                           seg_dist_loads, seg_thata, seg_delta)

            for j in range(3):

                if abs_function_direction == "max":
                    if d[j] > dmax[j]:
                        dmax[j] = d[j]
                elif abs_function_direction == "min":
                    if d[j] < dmax[j]:
                        dmax[j] = d[j]
        j: int
        for j in range(3):

            if global_[j] is None:
                # noinspection PyTypeChecker
                global_[j] = dmax[j]
                governing_combo[j] = comboINDEX

            else:

                if abs_function_direction == "max":
                    # noinspection PyTypeChecker
                    if dmax[j] > global_[j]:
                        # noinspection PyTypeChecker
                        global_[j] = dmax[j]
                        governing_combo[j] = comboINDEX

                elif abs_function_direction == "min":
                    # noinspection PyTypeChecker
                    if dmax[j] < global_[j]:
                        # noinspection PyTypeChecker
                        global_[j] = dmax[j]
                        governing_combo[j] = comboINDEX

    return global_, governing_combo

def max_difflection(casses: list, members_L: np.ndarray, m_index: int, seg: list, seg_internal_loads: list,
                    seg_dist_loads: list, seg_thata: list, seg_delta: list, combo_indexs: list | None = None):
    """
    Finds the maximum deflection for the member and its governing load case.

    :param casses: list. Load Case Indexes.
    :param members_L: ndarray. Length of each member.
    :param m_index: int. Index of the member to find values for.
    :param seg: list. Holds basic data to do with the segment.
                shape: (# Members, # Casses, # segments: varys, 5: [seg_x1, seg_x2, seg_EIz, seg_EIy, seg_EA])
    :param seg_internal_loads: list. Holds the end reaction forces of each segment.
                               shape: (# Members, # Casses, # Segments: varys, 3, 2: [P, M])
    :param seg_dist_loads: list. Holds the distributed loads on the segment.
                           shape: (# Members, # Casses, # segments: varys, 3, 2: (w1, w2))
    :param seg_thata: list. Holds the slope at the start of the segment.
                      shape: (# Members, # Casses, # Segments: varys, 3)
    :param seg_delta: list. Holds the deflection at the start of the segment.
                      shape: (# Members, # Casses, # Segments: varys, 3)
    :param combo_indexs: list. Load case to find the axial force for. OPTIONAL
    :return: Deflection, Governing Casse
    """

    return _difflection_extrem_finder(casses, m_index, seg, seg_internal_loads, seg_dist_loads, seg_thata,
                                      seg_delta, "max", members_L, combo_indexs)

def min_difflection(casses: list, members_L: np.ndarray, m_index: int, seg: list, seg_internal_loads: list,
                    seg_dist_loads: list, seg_thata: list, seg_delta: list, combo_indexs: list | None = None):
    """
    Finds the minimum deflection for the member and its governing load case.

    :param casses: list. Load Case Indexes.
    :param members_L: ndarray. Length of each member.
    :param m_index: int. Index of the member to find values for.
    :param seg: list. Holds basic data to do with the segment.
                shape: (# Members, # Casses, # segments: varys, 5: [seg_x1, seg_x2, seg_EIz, seg_EIy, seg_EA])
    :param seg_internal_loads: list. Holds the end reaction forces of each segment.
                                shape: (# Members, # Casses, # Segments: varys, 3, 2: [P, M])
    :param seg_dist_loads: list. Holds the distributed loads on the segment.
                            shape: (# Members, # Casses, # segments: varys, 3, 2: (w1, w2))
    :param seg_thata: list. Holds the slope at the start of the segment.
                        shape: (# Members, # Casses, # Segments: varys, 3)
    :param seg_delta: list. Holds the deflection at the start of the segment.
                        shape: (# Members, # Casses, # Segments: varys, 3)
    :param combo_indexs: list. Load case to find the axial force for. OPTIONAL
    :return: Deflection, Governing Casse
    """

    return _difflection_extrem_finder(casses, m_index, seg, seg_internal_loads, seg_dist_loads, seg_thata, seg_delta,
                                     "min", members_L, combo_indexs)
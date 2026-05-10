"""
Helper functions used for the analysis of a 3D frame.
"""

import numpy as np

import frame3DSolver.fixedEndReactionsCalculaters as ferCalc
import frame3DSolver.memberForceSolvers as MFSolvers

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

def part_D(nodes_support):
    """
    Builds lists of unreleased and released degree of freedom.

    :param nodes_support: ndarray. What nodes are supported in what DOFs.
    :return:
        D_unknown: ndarray. A ndarray of the indices for the released DOFs.
        D_known: ndarray. A ndarray of the indices for the unreleased DOFs.
    """
    n_nodes = nodes_support.shape[0]
    dof_indices = np.arange(n_nodes * 6).reshape(n_nodes, 6)
    D_unknown = dof_indices[~nodes_support].ravel()
    D_known = dof_indices[nodes_support].ravel()
    return D_unknown, D_known

def prep_members(nodes_cord, members, members_releases, members_point_loads, members_dist_loads, num_c):
    """
    Prepars the members for an analysis. This function can be run at the start of an optimization and
    does not need to be rerun is only cross-sections are changed.


    :param nodes_cord: ndarray. Cordinates of the nodes [X, Y, Z]
    :param members:  list. [i_node, j_node, material_index, setCrossSectionProps]
    :param members_releases: list. What directions are relesed Bool.
                             [Dxi, Dyi, Dzi, Rxi, Ryi, Rzi, Dxj, Dyj, Dzj, Rxj, Ryj, Rzj]
    :param members_point_loads: list. Point loads applyed to members.
                                [case, [[x], [[Px, Py, Pz, Mx, My, Mz]]]]
    :param members_dist_loads: list. Distributed loads applyed to members.
                               [case, [[[x1, x2], [wx1, wx2, wy1, wy2, wz1, wz2]]]]
    :param num_c: int. Number of casses.
    :return:
        DOFs: ndarray. A 3D ndarray of the DOF indices for each member. shape: (# members, 12)
        L: ndarray. A ndarray of the lengths of the members. shape: (# members)
        member_unrelesed_DOFs: list: 3D list of unreleased DOFs for each member.
                            shape: (# members, varys: # unrelesed DOFs)
        member_relesed_DOFs: list: 3D list of released DOFs for each member. shape: (# members, varys: # released DOFs)
        T: ndarray: 4D ndarray of the transformation matrices for each member. shape: (# members, 12, 12)
    """

    L = _get_L(nodes_cord, members)
    members = members[:, :2]
    DOFs = _build_dof_vector(members)
    member_unrelesed_DOFs, member_relesed_DOFs = _member_part_D(members_releases)
    T = _get_member_t(nodes_cord, members, L)
    point_loads = _assemble_point_loads(members_point_loads, num_c)
    dist_loads = _assemble_dist_loads(members_dist_loads, num_c)
    return DOFs, L, member_unrelesed_DOFs, member_relesed_DOFs, T, point_loads, dist_loads

def _get_L(nodes_cord, members):
    """
    Builds an array of the lengths of the members.

    :param nodes_cord: ndarray. Cordinates of the nodes [X, Y, Z]
    :param members: list. [i_node, j_node, material_index, setCrossSectionProps]
    :return:
        L: ndarray. A ndarray of the lengths of the members.
    """

    dx = nodes_cord[members[:, 1], 0] - nodes_cord[members[:, 0], 0]
    dy = nodes_cord[members[:, 1], 1] - nodes_cord[members[:, 0], 1]
    dz = nodes_cord[members[:, 1], 2] - nodes_cord[members[:, 0], 2]
    return np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

def _build_dof_vector(nodes: np.ndarray):
    """
    Returns the flattened list of global DOF indices for the supplied nodes.

    :param nodes: ndarray. A 3D ndarray of the indices of the nodes to get the DOF indices for. shape: (# members, 2)
    :return: ndarray. A 3D ndarray of the DOF indices for each member. shape: (# members, 12)
    """
    dofs = nodes[:, :, None] * 6 + np.arange(6)
    return dofs.reshape(nodes.shape[0], -1)

def _member_part_D(releases: np.ndarray):
    """
    Builds lists of unreleased and released degree of freedom indices for each member.

    :param releases: ndarray. 3D ndarray of booleans of wether a DOF is released. true = released, false = unreleased.
                     shape: (# members, 12)
    :return:
        unreleased: list: 3D list of unreleased DOFs for each member. shape: (# members, varys: # unrelesed DOFs)
        released: list: 3D list of released DOFs for each member. shape: (# members, varys: # released DOFs)
    """

    idx = np.arange(12)
    unreleased = [idx[~r] for r in releases]
    released   = [idx[r]  for r in releases]
    return unreleased, released

def _get_member_t(nodes_cord, members, L: np.ndarray):
    """
    Builds an array of the transformation matrices for each member.

    :param nodes_cord: ndarray. Cordinates of the nodes [X, Y, Z]
    :param members: list. [i_node, j_node, material_index, setCrossSectionProps]
    :param L: ndarray. A ndarray of the lengths of the members. shape: (# members)
    :return: ndarray: 4D array of the transformation matrices for each member. shape: (# members, 12, 12)
    """

    i = members[:, 0]
    j = members[:, 1]
    Xi, Yi, Zi = nodes_cord[i].T
    Xj, Yj, Zj = nodes_cord[j].T
    # Local x-axis (direction cosines)
    x = np.stack([(Xj - Xi)/L, (Yj - Yi)/L, (Zj - Zi)/L], axis=1)
    # Initialize y, z
    n = len(members)
    y = np.zeros((n, 3))
    z = np.zeros((n, 3))
    # Masks
    vertical   = np.isclose(Xi, Xj) & np.isclose(Zi, Zj)
    horizontal = np.isclose(Yi, Yj) & ~vertical
    general    = ~(vertical | horizontal)
    # Vertical members
    up = Yj > Yi
    y[vertical & up]  = [-1, 0, 0]
    y[vertical & ~up] = [ 1, 0, 0]
    z[vertical]       = [0, 0, 1]
    # Horizontal members
    y[horizontal] = [0, 1, 0]
    z[horizontal] = np.cross(x[horizontal], y[horizontal])
    z[horizontal] /= np.linalg.norm(z[horizontal], axis=1)[:, None]
    # General members
    proj = np.stack([Xj - Xi, np.zeros(n), Zj - Zi], axis=1)
    z_temp1 = np.cross(proj, x)
    z_temp2 = np.cross(x, proj)
    cond = (Yj > Yi)[:, None]
    z[general] = np.where(cond[general], z_temp1[general], z_temp2[general])
    z[general] /= np.linalg.norm(z[general], axis=1)[:, None]
    y[general] = np.cross(z[general], x[general])
    y[general] /= np.linalg.norm(y[general], axis=1)[:, None]
    # Direction cosine matrices
    dirCos = np.stack([x, y, z], axis=1)  # (n, 3, 3)
    # Build transformation matrices
    T = np.zeros((n, 12, 12))
    for k in range(4): T[:, k*3:(k+1)*3, k*3:(k+1)*3] = dirCos
    return T


def get_global_fixed_end_reaction_vector(nodes_cord, members_dof, T: np.ndarray, D_unknown: np.ndarray,
                                         D_known: np.ndarray, fer_condensed_array, num_m: int, num_c: int):
    """
    Builds array of the global fixed end reaction vector for both fixed and released DOFs.

    :param nodes_cord: ndarray. Cordinates of the nodes [X, Y, Z]
    :param members_dof: ndarray. A 3D ndarray of the DOF indices for each member. shape: (# members, 12)
    :param T: ndarray. 4D array of the transformation matrices for each member. shape: (# members, 12, 12)
    :param D_unknown: ndarray. A ndarray of the indices for the released DOFs.
    :param D_known: ndarray. A ndarray of the indices for the unreleased DOFs.
    :param fer_condensed_array: ndarray. Condensed local fixed end reaction vector. shape: (# members, # Casses, 12, 1)
    :param num_m: int. Number of members.
    :param num_c: int. Number of cases.
    :return:
        FER1: ndarray. An array of the released DOFs fixed end reactions.
        FER2: ndarray. An array of the unreased DOFs fixed end reactions.
    """

    FER_array = _get_FER(fer_condensed_array, T)
    FER_stack = np.array([[np.asarray(FER_array[i][j], dtype=float).reshape(12)
                           for j in range(num_c)] for i in range(num_m)])  # shape: (# members, # Casses, 12)
    members_dof = np.asarray(members_dof)  # (num_m, 12)
    nDOF = len(nodes_cord) * 6
    FER_val = np.zeros((num_c, nDOF))
    dof_flat = members_dof.ravel()
    for j in range(num_c):np.add.at(FER_val[j], dof_flat, FER_stack[:, j, :].ravel())
    FER_val = FER_val[..., None]
    FER1, FER2 = _part_force_vector(FER_val, D_unknown, D_known)
    return FER1, FER2

def _part_force_vector(FER_val: np.ndarray, D_unknown: np.ndarray, D_known: np.ndarray):
    """
    Partitions force vectors into 2 vectors one of releced and one unrelesed DOFs.

    :param FER_val: ndarray. A ndarray of forces to be partioned.
    :param D_unknown: ndarray. A ndarray of the indices for the released DOFs.
    :param D_known: ndarray. A ndarray of the indices for the unreleased DOFs.
    :return:
        FER1: ndarray. An array of the released DOFs forces.
        FER2: ndarray. An array of the unreased DOFs forces.
    """

    FER_val = np.asarray(FER_val)  # shape: (num_c, nDOF, 1)
    D_unknown = np.asarray(D_unknown, dtype=int)
    D_known = np.asarray(D_known, dtype=int)
    FER1 = FER_val[:, D_unknown, :]
    FER2 = FER_val[:, D_known, :]
    return FER1, FER2

def get_k_local_array(materials, members, members_CrossSectionProps, L: np.ndarray, log):
    """
    Builds an array of the local stiffness matrices for each member.

    :param materials: list. [E, G, nu, rho, fy]
    :param members: list. [i_node, j_node, material_index, setCrossSectionProps]
    :param members_CrossSectionProps: ndarray. [A, Iy, Iz, J]
    :param L: ndarray. A ndarray of the lengths of the members. shape: (# members)
    :return: 4D array of the local stiffness matrices for each member. shape: (# members, 12, 12)
    """

    E = materials[members[:, 2], 0]
    G = materials[members[:, 2], 1]
    A = members_CrossSectionProps[:, 0]
    Iy = members_CrossSectionProps[:, 1]
    Iz = members_CrossSectionProps[:, 2]
    J = members_CrossSectionProps[:, 3]

    n = len(L)
    # Precompute terms
    AE_L   = A * E / L
    EIy_L  = E * Iy / L
    EIz_L  = E * Iz / L
    EIy_L2 = E * Iy / L**2
    EIz_L2 = E * Iz / L**2
    EIy_L3 = E * Iy / L**3
    EIz_L3 = E * Iz / L**3
    GJ_L   = G * J / L
    # Allocate
    k = np.zeros((n, 12, 12))
    # Fill matrix (vectorized indexing)
    # Row 0 / 6
    k[:, 0, 0] = AE_L
    k[:, 0, 6] = -AE_L
    k[:, 6, 0] = -AE_L
    k[:, 6, 6] = AE_L
    # Row 1 / 7
    k[:, 1, 1] = 12 * EIz_L3
    k[:, 1, 5] = 6 * EIz_L2
    k[:, 1, 7] = -12 * EIz_L3
    k[:, 1, 11] = 6 * EIz_L2
    k[:, 7, 1] = -12 * EIz_L3
    k[:, 7, 5] = -6 * EIz_L2
    k[:, 7, 7] = 12 * EIz_L3
    k[:, 7, 11] = -6 * EIz_L2
    # Row 2 / 8
    k[:, 2, 2] = 12 * EIy_L3
    k[:, 2, 4] = -6 * EIy_L2
    k[:, 2, 8] = -12 * EIy_L3
    k[:, 2, 10] = -6 * EIy_L2
    k[:, 8, 2] = -12 * EIy_L3
    k[:, 8, 4] = 6 * EIy_L2
    k[:, 8, 8] = 12 * EIy_L3
    k[:, 8, 10] = 6 * EIy_L2
    # Row 3 / 9
    k[:, 3, 3] = GJ_L
    k[:, 3, 9] = -GJ_L
    k[:, 9, 3] = -GJ_L
    k[:, 9, 9] = GJ_L
    # Row 4 / 10
    k[:, 4, 2] = -6 * EIy_L2
    k[:, 4, 4] = 4 * EIy_L
    k[:, 4, 8] = 6 * EIy_L2
    k[:, 4, 10] = 2 * EIy_L
    k[:, 10, 2] = -6 * EIy_L2
    k[:, 10, 4] = 2 * EIy_L
    k[:, 10, 8] = 6 * EIy_L2
    k[:, 10, 10] = 4 * EIy_L
    # Row 5 / 11
    k[:, 5, 1] = 6 * EIz_L2
    k[:, 5, 5] = 4 * EIz_L
    k[:, 5, 7] = -6 * EIz_L2
    k[:, 5, 11] = 2 * EIz_L
    k[:, 11, 1] = 6 * EIz_L2
    k[:, 11, 5] = 2 * EIz_L
    k[:, 11, 7] = -6 * EIz_L2
    k[:, 11, 11] = 4 * EIz_L

    if log:
        print("local k: ", k)
    return k

def member_part_k_array(k_array: np.ndarray, R_unrelesed_array: list or np.ndarray,
                        R_relesed_array: list or np.ndarray, numM: int):
    """
    Partitions the stiffness matrices into sub-matrices based on unreleased and released degree of freedom indices.

    :param k_array: ndarray. 4D array of the stiffness matrices for each member. shape: (# members, 12, 12)
    :param R_unrelesed_array: list. 3D list of unreleased DOFs for each member.
                              shape: (# members, varys: # unrelesed DOFs)
    :param R_relesed_array: list. 3D list of released DOFs for each member. shape: (# members, varys: # released DOFs)
    :param numM: int. Number of members.
    :return:  4 lists of sub-matrices based on the degrese of freedom.
              general shape: (# members, ndarray(varys, varys))
    """

    k11_array = []
    k12_array = []
    k21_array = []
    k22_array = []
    for i in range(numM):
        k = np.array(k_array[i])
        R_unrelesed = np.array(R_unrelesed_array[i], dtype=int)
        R_relesed = np.array(R_relesed_array[i], dtype=int)
        k11_array.append(k[R_unrelesed, :][:, R_unrelesed])
        k12_array.append(k[R_unrelesed, :][:, R_relesed])
        k21_array.append(k[R_relesed, :][:, R_unrelesed])
        k22_array.append(k[R_relesed, :][:, R_relesed])
    return k11_array, k12_array, k21_array, k22_array

def _get_FER(fer: np.ndarray, T: np.ndarray):
    """
    Returns the global fixed end reaction vectors for all members and casses

    :param fer: ndarry. The array of local fer vectors to be transformed. shape: (# members, # Casses, 12, 1)
    :param T: ndarray. 4D array of the transformation matrices for each member. shape: (# members, 12, 12)
    :return: global fixed end reaction vectors for all members and casses
    """
    invT = T.transpose(0, 2, 1) # T⁻¹ = Tᵀ for orthogonal transformation matrices
    FER = invT[:, None] @ fer
    return FER

def _assemble_point_loads(members_point_loads, num_c):
    """
    Assembles the point load vectors for all members and casses

    :param members_point_loads:  list. Point loads applyed to members. [case, [[x], [[Px, Py, Pz, Mx, My, Mz]]]]
    :param num_c: int. Number of casses
    :return: list. Point loads applied to each member.
             shape: (# members, # Casses, # loads in case on member, 7: (x, Px, Py, Pz, Mx, My, Mz))
    """

    new_point_loads = []
    for i, load_casses in enumerate(members_point_loads):
        new_loads  = [None] * num_c
        for loads in load_casses:
            new_load  = []
            for j in range(len(loads[1][0])):
                new_load.append([loads[1][0][j]] + loads[1][1][j])
            new_loads[loads[0]] = new_load
        for j in range(len(new_loads)):
            if new_loads[j] is None:
                # noinspection PyTypeChecker
                new_loads[j] = []
        new_point_loads.append(new_loads)
    return new_point_loads

def _assemble_dist_loads(members_dist_loads, num_c):
    """
    Assembles the distributed load vectors for all members and casses

    :param members_dist_loads: list. Distributed loads applyed to members.
                              [case, [[[x1, x2], [wx1, wx2, wy1, wy2, wz1, wz2]]]]
    :param num_c: int. Number of casses
    :return: list. Distributed loads applied to each member.
             shape: (# members, # Casses, varys, # loads in case on member: (x1, x2, wx1, wx2, wy1, wy2, wz1, wz2))
    """

    new_dist_loads = []
    for i, load_casses in enumerate(members_dist_loads):
        new_loads = [None] * num_c
        for loads in load_casses:
            new_load = []
            for j in range(len(loads[1])):
                new_load.append(loads[1][j][0] + loads[1][j][1])
            new_loads[loads[0]] = new_load
        for j in range(len(new_loads)):
            if new_loads[j] is None:
                # noinspection PyTypeChecker
                new_loads[j] = []
        new_dist_loads.append(new_loads)
    return new_dist_loads

def get_member_fer_unc(members_L, point_loads: list, dist_loads: list, num_m: int, num_c: int):
    """
    Returns the local fixed end reaction vector for the members, ignoring the effects of end releases.

    :param members_L: ndarray. length of each member.
    :param point_loads: list. Point loads applied to each member.
                         shape: (# members, # Casses, # loads in case on member, 7: (x, Px, Py, Pz, Mx, My, Mz))
    :param dist_loads: list. Distributed loads applied to each member.
                       shape: (# members, # Casses, varys,
                       # loads in case on member: (x1, x2, wx1, wx2, wy1, wy2, wz1, wz2))
    :param num_m: int. Number of members.
    :param num_c: int. Number of Casses.
    :return: ndarray. local fixed end reaction vector for the members. shape: (# members, # Casses, 12, 1)
    """
    fer_unc = np.zeros((num_m, num_c, 12, 1))
    fer_unc = np.add(fer_unc, _get_fixed_end_reactions_point_load_array(members_L, point_loads))
    fer_unc = np.add(fer_unc, _get_fixed_end_reactions_dist_load_array(members_L, dist_loads))
    return fer_unc

def _get_fixed_end_reactions_point_load_array(members_L, point_loads: list):
    """
    Builds an array of the fixed end reacctions for point loads.

    :param members_L: ndarray. length of each member.
    :param point_loads: list. Point loads applied to each member.
                        shape: (# members, # Casses, varys: # loads in case on member, 7: (x, Px, Py, Pz, Mx, My, Mz))
    :return: ndarray. fixed end reactions for point loads. shape: (# members, # Casses, 12, 1)
    """

    members_L = np.asarray(members_L)
    reactions = []
    for i, loadCases in enumerate(point_loads):
        L = members_L[i]
        reactionsMember = []
        for loads in loadCases:
            if loads is None or len(loads) == 0:
                reactionsMember.append(np.zeros((12, 1)))
                continue
            loads = np.asarray(loads)
            x  = loads[:, 0]
            Fx = loads[:, 1]
            Fy = loads[:, 2]
            Fz = loads[:, 3]
            Mx = loads[:, 4]
            My = loads[:, 5]
            Mz = loads[:, 6]
            R = np.zeros((12,))
            if np.any(Fx): R += ferCalc.pointLoadX_batch(Fx, x, L).sum(axis=0)
            if np.any(Fy): R += ferCalc.pointLoadY_batch(Fy, x, L).sum(axis=0)
            if np.any(Fz): R += ferCalc.pointLoadZ_batch(Fz, x, L).sum(axis=0)
            if np.any(Mx): R += ferCalc.momentX_batch(Mx, x, L).sum(axis=0)
            if np.any(My): R += ferCalc.momentY_batch(My, x, L).sum(axis=0)
            if np.any(Mz): R += ferCalc.momentZ_batch(Mz, x, L).sum(axis=0)
            reactionsMember.append(R.reshape(12, 1))
        reactions.append(reactionsMember)
    return reactions

def _get_fixed_end_reactions_dist_load_array(members_L, dist_loads):
    """
    Builds an array of the fixed end reacctions for distributed loads.

    :param members_L: ndarray. length of each member.
    :param dist_loads: list. Distributed loads applied to each member.
                       shape: (# members, # Casses, varys,
                       varys: # loads in case on member: (x1, x2, wx1, wx2, wy1, wy2, wz1, wz2))
    :return: ndarray. fixed end reactions for distributed loads. shape: (# members, # Casses, 12, 1)
    """

    reactions = []
    for i, load_casses in enumerate(dist_loads):
        L = members_L[i]
        reactions_member = []
        for j, loads in enumerate(load_casses):
            reactions = np.zeros((12, 1))
            for k, load in enumerate(loads):
                if not (load[2] == 0 and load[3] == 0):
                    reactions = reactions + ferCalc.distributedLoadX([load[2],load[3]], [load[0], load[1]], L)
                if not (load[4] == 0 and load[5] == 0):
                    reactions = reactions + ferCalc.distributedLoadY([load[4],load[5]], [load[0], load[1]], L)
                if not (load[6] == 0 and load[7] == 0):
                    reactions = reactions + ferCalc.distributedLoadZ([load[6],load[7]], [load[0], load[1]], L)
            reactions_member.append(reactions)
        reactions.append(reactions_member)
    return reactions

def member_part_fer(fer: np.ndarray, R_unrelesed: np.ndarray or list, R_relesed: np.ndarray or list,
                    num_m: int, num_c: int):
    """
    Partition fer into 2 sub matrices based on unreleased and released degree of freedom.

    :param fer: ndarray. local fixed end reaction vector for the members. shape: (# members, # Casses, 12, 1)
    :param R_unrelesed: ndarray or list. list: 3D list of unreleased DOFs for each member.
                        shape: (# members, varys: # unrelesed DOFs)
    :param R_relesed: ndarray or list. list: 3D list of released DOFs for each member.
                      shape: (# members, varys: # released DOFs)
    :param num_m: int. Number of members.
    :param num_c: int. Number of Casses.
    :return:
        fer1: list. fixed end reactions for each member for the unrelesed DOFs.
              shape: (# members, # Casses, varys: # unrelesed DOFs, 1)
        fer2: list. fixed end reactions for eacch member for the released FODs.
              shape: (# members, # Casses, varys: # relesed DOFs, 1)
    """

    fer1_array = []
    fer2_array = []
    for i in range(num_m):
        fer1_sub_array = []
        fer2_sub_array = []
        for j in range(num_c):
            fer1 = fer[i,j,R_unrelesed[i]]
            fer2 = fer[i,j,R_relesed[i]]
            fer1_sub_array.append(fer1)
            fer2_sub_array.append(fer2)
        fer1_array.append(fer1_sub_array)
        fer2_array.append(fer2_sub_array)
    return fer1_array, fer2_array

def get_fer(members_releases, k12: list, k22: list, fer1: np.ndarray, fer2: np.ndarray, num_m: int, num_c: int):
    """
    Returns the condensed local fixed end reaction vector for all members and casses.

    :param members_releases: list. What directions are relesed Bool.
                             [Dxi, Dyi, Dzi, Rxi, Ryi, Rzi, Dxj, Dyj, Dzj, Rxj, Ryj, Rzj]
    :param k12: list. sub-matrix of the member local stiffness matrix. shape: (# members, varys)
    :param k22: list. sub-matrix of the member local stiffness matrix. shape: (# members, varys)
    :param fer1: ndarray. fixed end reactions for each member for the unrelesed DOFs.
                 shape: (# members, # Casses, varys: # unrelesed DOFs, 1)
    :param fer2: ndarray. fixed end reactions for eacch member for the released FODs.
                 shape: (# members, # Casses, varys: # relesed DOFs, 1)
    :param num_m: int. Number of members.
    :param num_c: int. Number of casses.
    :return: ndarray. Condensed local fixed end reaction vector. shape: (# members, # Casses, 12, 1)
    """

    fer_condenced_array = []
    for i in range(num_m):
        fer_condensed_sub_array = []
        for j in range(num_c):
            fer_1 = np.array(fer1[i][j])
            fer_2 = np.array(fer2[i][j])
            k_12 = np.array(k12[i])
            k_22 = np.array(k22[i])
            fer_condensed = (fer_1 - k_12 @ k_22 @ fer_2)
            a = 0
            for DOF in members_releases[i]:
                if DOF:
                    fer_condensed = np.insert(fer_condensed, a, 0, axis=0)
                a += 1
            fer_condensed_sub_array.append(fer_condensed)
        fer_condenced_array.append(fer_condensed_sub_array)
    return fer_condenced_array

def get_parted_global_nodal_force_vector(nodes_loads, casses, D_unknown, D_known, num_n: int):
    """
    Builds the pertitioned global nodal force vector.

    :param nodes_loads: list. [case, [Px, Py, Pz, Mx, My, Mz]]
    :param casses: list. Load Case Indexes.
    :param D_unknown: ndarray. A ndarray of the indices for the released DOFs.
    :param D_known: ndarray. A ndarray of the indices for the unreleased DOFs.
    :param num_n: int. Number of nodes
    :return: ndarray. Global nodal force vector partitioned into 2 vectors.
    """

    p_array = []
    for load_case_i in casses:
        p = np.zeros((num_n * 6, 1))
        for i in range(num_n):
            if np.size(nodes_loads[load_case_i]) != 0:
                local = np.array(nodes_loads[load_case_i][1])
            else:
                local = np.zeros(6, dtype=float)
            dofs = np.empty(len([i]) * 6, dtype=np.int64)
            local_ = np.arange(6, dtype=np.int64)
            for j, node in enumerate([i]):
                start = j * 6
                dofs[start:start + 6] = node * 6 + local_
            p[dofs, 0] += local
        p_array.append(p)
    return _part_force_vector(np.array(p_array), D_unknown, D_known)

def k_member_make_global(k_local_array, T_array):
    """
    Converts member stiffness matrices from the local cordinate system to the global cordinate system.

    :param k_local_array: ndarray. 4D array of the local stiffness matrices for each member. shape: (# members, 12, 12)
    :param T_array: ndarray: 4D array of the transformation matrices for each member. shape: (# members, 12, 12)
    :return: ndarray: 4D array of the global stiffness matrices for each member. shape: (# members, 12, 12)
    """

    k_member_global = []
    for i in range(len(T_array)):
        k_member_global.append(np.linalg.inv(T_array[i]) @ k_local_array[i] @ T_array[i])
    return k_member_global

def get_K_global(members_dof, k_global: np.ndarray, num_n: int, num_m:int):
    """
    Assembles the global stiffness matrix for the frame.

    :param members_dof: ndarray. A 3D ndarray of the DOF indices for each member. shape: (# members, 12)
    :param k_global: ndarray. 4D array of the global stiffness matrices for each member. shape: (# members, 12, 12)
    :param num_n: int. Number of nodes.
    :param num_m: int. Number of casses.
    :return: ndarray. Global stiffness matrix. shape: (# nodes * 6, # nodes * 6)
    """

    K = np.zeros((num_n * 6, num_n * 6))
    for i in range(num_m):
        K[np.ix_(members_dof[i], members_dof[i])] += np.asarray(k_global[i], dtype=float)
    return K

def partition_K_gloabl(K_global: np.ndarray, R_unrelesed: np.ndarray, R_relesed: np.ndarray):
    """
    Partitiones the global stiffness matrix with respect to the frames supports.

    :param K_global: ndarray. Global stiffness matrix. shape: (# nodes * 6, # nodes * 6)
    :param R_unrelesed: ndarray. A ndarray of the indices for the unreleased DOFs. (not supported)
    :param R_relesed: ndarray. A ndarray of the indices for the unreleased DOFs. (supported)
    :return: ndarray. 4 sub-matrices of the global stiffness matrix with respect to the frames supports.
    """

    K11 = (K_global[R_unrelesed, :][:, R_unrelesed])
    K12 = (K_global[R_unrelesed, :][:, R_relesed])
    K21 = (K_global[R_relesed, :][:, R_unrelesed])
    K22 = (K_global[R_relesed, :][:, R_relesed])
    return K11, K12, K21, K22

def get_D(K11: np.ndarray, K12:np.ndarray, P1_array:np.ndarray, FER1_array:np.ndarray, index_unsupported:np.ndarray,
          index_supported: np.ndarray, num_n: int, num_c:int, log):
    """
    Solves for the nodal displacement for the frame.

    :param K11: ndarray. Sub-matrices of the global stiffness matrix.
    :param K12: ndarray. Sub-matrices of the global stiffness matrix.
    :param P1_array: ndarray. Vector of the nodel loads.
    :param FER1_array: ndarray. An array of the released DOFs fixed end reactions.
    :param index_unsupported: ndarray. Indeces of the unsupported DOFs.
    :param index_supported: ndarray. Indeces of the supported DOFs.
    :param num_n: int. Number of nodes.
    :param num_c: int. Number of casses.
    :return:
        D: ndarray. Array of the nodal displacements.
        DX: ndarray. Array of the nodel displacements in the X direction.
        DY: ndarray. Array of the nodel displacements in the Y direction.
        DZ: ndarray. Array of the nodel displacements in the Z direction.
        RX: ndarray. Array of the nodel rotation displacements in the X direction.
        RY: ndarray. Array of the nodel rotation displacements in the Y direction.
        RZ: ndarray. Array of the nodel rotation displacements in the Z direction.
    """

    D1_array = []
    D2 = np.zeros((len(index_supported), 1))
    for i in range(num_c):
        P1 = np.array(P1_array[i])
        FER1 = np.array(FER1_array[i])
        if K11.shape == (0, 0):
            D1_array.append([])
        else:
            try:
                D1_array.append(np.linalg.solve(K11, np.subtract(np.subtract(P1, FER1), np.matmul(K12, D2))))
            except:
                raise Exception('The stiffness matrix is singular, which implies rigid body motion. '
                                'The structure is unstable. Aborting analysis.')
    D_array = _assemble_D_array(D1_array, D2, index_supported, index_unsupported, num_n, num_c)

    if log:
        print("D_array: ", D_array)

    DX_array, DY_array, DZ_array, RX_array, RY_array, RZ_array = _get_node_direction_deflections(D_array, num_n, num_c)

    if log:
        print("DX_array: ", DX_array)
        print("DY_array: ", DY_array)
        print("DZ_array: ", DZ_array)
        print("RX_array: ", RX_array)
        print("RY_array: ", RY_array)
        print("RZ_array: ", RZ_array)

    return D_array, DX_array, DY_array, DZ_array, RX_array, RY_array, RZ_array

def _assemble_D_array(D1_array, D2, index_supported, index_unsupported, num_n: int, num_c: int):
    """

    :param D1_array: ndarray. Array of displecements for all unspported DOFs.
    :param D2: ndarray. Array of displacements for all supported DOFs.
    :param index_supported: ndarray. Indeces of the supported DOFs.
    :param index_unsupported: ndarray. Indeces of the unsupported DOFs.
    :param num_n: int: Number of nodes.
    :param num_c: intL Number of casses.
    :return: ndarray. Array of the nodal displacements.
    """

    D_array = []
    for a in range(num_c):
        D1 = np.array(D1_array[a])
        D = np.zeros((num_n * 6, 1))
        for j in range(num_n):
            for i in range(6):
                if j * 6 + i in index_supported:
                    D[(j * 6 + i, 0)] = D2[index_supported.tolist().index(j * 6 + i), 0]
                else:
                    D[(j * 6 + i, 0)] = D1[index_unsupported.tolist().index(j * 6 + i), 0]
        D_array.append(D)
    return D_array

def _get_node_direction_deflections(D_array, num_n: int, num_c: int):
    """
    Gets the nodal deflection in each direction.

    :param D_array: ndarray. Array of the nodal displacements.
    :param num_n: int: Number of nodes.
    :param num_c: int: Number of casses.
    :return:
        DX: ndarray. Array of the nodel displacements in the X direction.
        DY: ndarray. Array of the nodel displacements in the Y direction.
        DZ: ndarray. Array of the nodel displacements in the Z direction.
        RX: ndarray. Array of the nodel rotation displacements in the X direction.
        RY: ndarray. Array of the nodel rotation displacements in the Y direction.
        RZ: ndarray. Array of the nodel rotation displacements in the Z direction.
    """

    DX_array = []
    DY_array = []
    DZ_array = []
    RX_array = []
    RY_array = []
    RZ_array = []
    for i in range(num_c):
        DX_case = []
        DY_case = []
        DZ_case = []
        RX_case = []
        RY_case = []
        RZ_case = []
        for j in range(num_n):
            DX_case.append(D_array[i][j * 6 + 0, 0])
            DY_case.append(D_array[i][j * 6 + 1, 0])
            DZ_case.append(D_array[i][j * 6 + 2, 0])
            RX_case.append(D_array[i][j * 6 + 3, 0])
            RY_case.append(D_array[i][j * 6 + 4, 0])
            RZ_case.append(D_array[i][j * 6 + 5, 0])
        DX_array.append(DX_case)
        DY_array.append(DY_case)
        DZ_array.append(DZ_case)
        RX_array.append(RX_case)
        RY_array.append(RY_case)
        RZ_array.append(RZ_case)
    return (np.array(DX_array), np.array(DY_array), np.array(DZ_array), np.array(RX_array), np.array(RY_array),
            np.array(RZ_array))

def get_member_direction_deflections(members, DX_array: np.ndarray, DY_array: np.ndarray, DZ_array: np.ndarray,
                                     RX_array: np.ndarray, RY_array: np.ndarray, RZ_array: np.ndarray,
                                     num_m: int, num_c: int):
    """
    Gets the global end deflections for each member.

    :param members:  list. [i_node, j_node, material_index, setCrossSectionProps]
    :param DX_array: ndarray. Array of the nodel displacements in the X direction.
    :param DY_array: ndarray. Array of the nodel displacements in the Y direction.
    :param DZ_array: ndarray. Array of the nodel displacements in the Z direction.
    :param RX_array: ndarray. Array of the nodel rotation displacements in the X direction.
    :param RY_array: ndarray. Array of the nodel rotation displacements in the Y direction.
    :param RZ_array: ndarray. Array of the nodel rotation displacements in the Z direction.
    :param num_m: int. Number of members.
    :param num_c: int. Number of casses.
    :return: ndarray. Array of the global end deflections for each member. shape: (# casses, # members, 12)
    """

    D_member_array = []
    for i in range(num_c):
        D_sub = []
        for j in range(num_m):
            Dx1 = DX_array[i, members[j,0]]
            Dy1 = DY_array[i, members[j,0]]
            Dz1 = DZ_array[i, members[j,0]]
            Rx1 = RX_array[i, members[j,0]]
            Ry1 = RY_array[i, members[j,0]]
            Rz1 = RZ_array[i, members[j,0]]
            Dx2 = DX_array[i, members[j, 1]]
            Dy2 = DY_array[i, members[j, 1]]
            Dz2 = DZ_array[i, members[j, 1]]
            Rx2 = RX_array[i, members[j, 1]]
            Ry2 = RY_array[i, members[j, 1]]
            Rz2 = RZ_array[i, members[j, 1]]
            D_sub.append(np.array([Dx1, Dy1, Dz1, Rx1, Ry1, Rz1, Dx2, Dy2, Dz2, Rx2, Ry2, Rz2]))
        D_member_array.append(D_sub)
    return np.array(D_member_array)

def get_d(T_array: np.ndarray, D_array: np.ndarray, num_m: int, num_c: int):
    """
    Transforms the end deflection array for each member from the global cordinate system into the local
    cordinate system.

    :param T_array: ndarray. 4D ndarray of the transformation matrices for each member. shape: (# members, 12, 12)
    :param D_array: ndarray. Array of the global end deflections for each member. shape: (# casses, # members, 12)
    :param num_m: int. Number of members.
    :param num_c: int. Number of casses.
    :return: ndarray. Array of the local end deflections for each member. shape: (# casses, # members, 12)
    """

    d = []
    for i in range(num_c):
        d_sub = []
        for j in range(num_m):
            d_sub.append(T_array[j] @ D_array[i][j])
        d.append(d_sub)
    return np.array(d)

def get_F(T_array: np.ndarray, f_array: np.ndarray, num_m: int, num_c:int):
    """
    Gets the global forces acting at each end of the members.

    :param T_array: ndarray. 4D ndarray of the transformation matrices for each member. shape: (# members, 12, 12)
    :param f_array: ndarray. Array of the local forces acting at the ends of the members.
                    shape: (# casses, # members, 12)
    :param num_m: int. Number of members.
    :param num_c: int. Number of casses.
    :return: ndarray. Array of the global forces acting at the ends of the members. shape: (# casses, # members, 12)
    """

    F = []
    for i in range(num_c):
        F_sub = []
        for j in range(num_m):
            F_sub.append(np.linalg.inv(T_array[j]) @ f_array[i, j])
        F.append(F_sub)
    return np.array(F)

def get_f(k_array: np.ndarray, d_array: np.ndarray, fer_array: np.ndarray, num_m: int, num_c: int):
    """
    Gets the local forces acting at each end of the members.

    :param k_array: ndarray. 4D array of the local stiffness matrices for each member. shape: (# members, 12, 12)
    :param d_array: ndarray. Array of the local end deflections for each member. shape: (# casses, # members, 12)
    :param fer_array: ndarray. Condensed local fixed end reaction vector. shape: (# members, # Casses, 12, 1)
    :param num_m: int: Number of members.
    :param num_c: int: Number of casses.
    :return: ndarray. Array of the local forces acting at the ends of the members. shape: (# casses, # members, 12)
    """

    f = []
    for i in range(num_c):
        f_sub = []
        for j in range(num_m):
            f_sub.append(k_array[j] @ d_array[i, j] + fer_array[j, i])
        f.append(f_sub)
    return np.array(f)

def get_weight(materials, members, members_L, members_cross_section_props):
    """
    Gets the weight of each member of the frame.

    :param materials:  list. [E, G, nu, rho, fy]
    :param members: list. [i_node, j_node, material_index, setCrossSectionProps]
    :param members_L: ndarray. length of each member.
    :param members_cross_section_props: ndarray. [A, Iy, Iz, J]
    :return: ndarry. Array of the weights of each member of the frame.
    """

    return materials[members[:,2],3] * members_cross_section_props[:, 0] * members_L

def get_reactions(nodes_support, nodes_loads, members, members_releases, F_array, num_c:int, num_m:int, num_n:int):
    """
    Gets the reactions for the frame. If an DOF is not supported the reaction is 0.

    :param nodes_support: list. [support_DX, support_DY, support_DZ, support_RX, support_RY, support_RZ]
    :param nodes_loads:  list. [case, [Px, Py, Pz, Mx, My, Mz]]
    :param members: list. [i_node, j_node, material_index, setCrossSectionProps]
    :param members_releases: list. What directions are relesed Bool.
                             [Dxi, Dyi, Dzi, Rxi, Ryi, Rzi, Dxj, Dyj, Dzj, Rxj, Ryj, Rzj]
    :param F_array: ndarray. Array of the global forces acting at the ends of the members.
                    shape: (# casses, # members, 12)
    :param num_c: int: Number of casses.
    :param num_m: int: Number of members.
    :param num_n: int: Number of nodes.
    :return: ndarray. Array of the reactions at the supports of the frame. shape: (# casses, # nodes, 6)
    """

    reactions = []
    for i in range(num_c):
        reactions_sub = []
        for j in range(num_n):
            r = [0.0] * 6
            for s in range(6):
                if nodes_support[i][s]:
                    for k in range(num_m):
                        if members[k, 0] == j:
                            if not members_releases[k][s]: r[s] += F_array[i, k, s, 0]
                        elif members[k, 1] == j:
                            if not members_releases[k][s + 6]: r[s] += F_array[i, k, s + 6, 0]
                    for joint_load_case in nodes_loads[j]:
                        if joint_load_case[0] == i:
                            r[s] += joint_load_case[1][s]
            reactions_sub.append(r)
        reactions.append(reactions_sub)
    return reactions

#todo
def solve_internal_forces(members, members_L, members_cross_section_props, materials, casses, point_loads, dist_loads,
                          f_array, fer_array, d_array, num_m, num_c):
    segments = MFSolvers.segment_Member(members, members_L, members_cross_section_props, materials, point_loads,
                                        dist_loads, f_array, fer_array, d_array, num_m, num_c)

    seg, seg_internal_loads, seg_dist_loads, seg_thata, seg_delta = segments
    abs_F = []
    abs_M = []
    for mINDEX in range(num_m):
        FX_min, _ = MFSolvers.min_axial(casses, mINDEX, seg, seg_internal_loads, seg_dist_loads)
        FX_max, _ = MFSolvers.max_axial(casses, mINDEX, seg, seg_internal_loads, seg_dist_loads)
        FY_min, FZ_min, _, _ = MFSolvers.min_shear(casses, mINDEX, seg, seg_internal_loads, seg_dist_loads)
        FY_max, FZ_max, _, _ = MFSolvers.max_shear(casses, mINDEX, seg, seg_internal_loads, seg_dist_loads)
        MX_min, _ = MFSolvers.min_tourque(casses, mINDEX, seg, seg_internal_loads, seg_dist_loads)
        MX_max, _ = MFSolvers.max_tourque(casses, mINDEX, seg, seg_internal_loads, seg_dist_loads)
        MY_min, MZ_min, _, _ = MFSolvers.min_moment(casses, mINDEX, seg, seg_internal_loads, seg_dist_loads)
        MY_max, MZ_max, _, _ = MFSolvers.max_moment(casses, mINDEX, seg, seg_internal_loads, seg_dist_loads)
        abs_F.append([max(FX_max,abs(FX_min)), max(FY_max,abs(FY_min)), max(FZ_max,abs(FZ_min))])
        abs_M.append([max(MX_max, abs(MX_min)), max(MY_max, abs(MY_min)), max(MZ_max, abs(MZ_min))])
    return abs_F, abs_M



"""-----------------Optimization-----------------------------"""

from frame3DSolver.CrossSectionCalculaters import SquareHSS, TubeHSS, Angle, RectHSS
import scipy.optimize as opt

# noinspection PyUnusedLocal
def _cost(D, DX, DY, DZ, RX, RY, RZ, Weight, Reactions, InternalForces, cost_function):
    return eval(cost_function)

def get_cost(X, constants):
    """
    Gets the _cost for the current optimization variables.

    :param X: list. Optimization variables.
    :param constants: list. arguments needed to calculate the _cost: [cost_function, frame, member_group, member_group_type,
                     weight_needed, reactions_needed, internal_forces_needed, log]
    :return: float. _cost
    """

    (cost_function, frame, member_group, member_group_type,
     weight_needed, reactions_needed, internal_forces_needed, log) = constants

    cross_section_props = get_cross_section_props(X, member_group, member_group_type)
    if log: print("Variable cross section properties: ", cross_section_props)

    j = 0
    for i in range(len(frame.members_CrossSectionProps)):
        if not frame.members[i][3]:
            frame.members_CrossSectionProps[i] = cross_section_props[j]
            j += 1
    if log: print("Cross Seciton Properites: ", frame.members_CrossSectionProps)

    D, DX, DY, DZ, RX, RY, RZ, weight, reactions, internalForces = (
        frame.analysis_linear(weight_needed, reactions_needed, internal_forces_needed, log))

    return _cost(D, DX, DY, DZ, RX, RY, RZ, weight, reactions, internalForces, cost_function)

def chack_inputs(members, member_group: list, member_group_type: list):
    """
    Checks if the input lists for the optimization are valid.

    :param members: list. [i_node, j_node, material_index, setCrossSectionProps]
    :param member_group: list. list of indices of member groups for non set members to be assigned to.
    :param member_group_type: list. List of cross-section types for each member group to be assigned.
    """

    number_not_set_members = 0
    for member in members:
        if not member[3]:
            number_not_set_members += 1
    if len(member_group) != number_not_set_members:
        raise Exception('Number of not set members is not equal to number of members assigned to a group.')
    if max(member_group) + 1 != len(member_group_type):
        raise Exception('Number of member groups does not equal number of groups members are assigned to')

def get_num_varables(member_group_type: list):
    """
    Gets the number of variables to be in the optimization problem.

    :param member_group_type: list. List of cross-section types for each member group to be assigned.
    :return:
        numVariables: int, number of variables in the optimization.
    """

    num_varables = 0
    for t in member_group_type:
        if t == "Angle":
            num_varables += 3
        elif t == "RectHSS":
            num_varables += 3
        elif t == "SquareHSS":
            num_varables += 2
        elif t == "TubeHSS":
            num_varables += 2
        else:
            raise Exception('member groupe type is not a valid cross section type. '
                            'Chose one of the following: Angle, RectHSS, SquareHSS, TubeHSS')
    return num_varables

def get_bounds(lower_bound: list or float, upper_bound: list or float, num_varables: int):
    """
    Sets up the scipy optimization bounds.

    :param lower_bound: list or float. Lower bound on the optimization variables.
                        if list must be length of number of variables.
    :param upper_bound: list or float. Upper bound on the optimization variables.
                        if list must be length of number of variables.
    :param num_varables: int, Number of variables to be in the optimization problem.
    :return: Scipy optimization bounds.
    """

    if isinstance(lower_bound, float):
        lower_bound = [lower_bound] * num_varables
    if isinstance (upper_bound, float):
        upper_bound = [upper_bound] * num_varables
    return opt.Bounds(lower_bound, upper_bound)

def get_cross_section_props(X, member_groups, member_group_type):
    """
    Gets the cross-section properties of the cross-section being optimized from the optimization variables.

    :param X: list: Optimization variables.
    :param member_groups: list. list of indices of member groups for non set members to be assigned to.
    :param member_group_type: list. List of cross-section types for each member group to be assigned.
    :return: list: Cross-section properties for the optimization cross-sections.
             shape: (# cross-sections being optimized, 4
    """

    group_properties = []
    j = 0
    for t in member_group_type:
        if t == "Angle":
            A = Angle.getA(X[j], X[j + 1], X[j + 2])
            Iy = Angle.getIy(X[j], X[j + 1], X[j + 2])
            Iz = Angle.getIx(X[j], X[j + 1], X[j + 2])
            J = Angle.getJ(X[j], X[j + 1], X[j + 2])
            group_properties.append([A, Iy, Iz, J])
            j += 3
        elif t == "RectHSS":
            A = RectHSS.getA(X[j], X[j + 1], X[j + 2])
            Iy = RectHSS.getIy(X[j], X[j + 1], X[j + 2])
            Iz = RectHSS.getIx(X[j], X[j + 1], X[j + 2])
            J = RectHSS.getJ(X[j], X[j + 1], X[j + 2])
            group_properties.append([A, Iy, Iz, J])
            j += 3
        elif t == "SquareHSS":
            A = SquareHSS.getA(X[j], X[j + 1])
            I = SquareHSS.getI(X[j], X[j + 1])
            J = SquareHSS.getJ(X[j], X[j + 1])
            group_properties.append([A, I, I, J])
            j += 2
        elif t == "TubeHSS":
            A = TubeHSS.getA(X[j], X[j + 1])
            I = TubeHSS.getI(X[j], X[j + 1])
            J = TubeHSS.getJ(X[j], X[j + 1])
            group_properties.append([A, I, I, J])
            j += 2

    member_properties = []
    for member_g in member_groups:
        member_properties.append(group_properties[member_g])

    return member_properties
"""
Optimizes the cross-sections of a 2D truss.
"""

import numpy as np
import cvxpy as cvx
import matplotlib.pyplot as plt
from scipy import sparse

from opening_saving.saving import save_truss_cross_section_optimization_excel

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

def truss_main(file_path, nodes, members, load_cases, supports):
    """
    Prepares the inputs, and runs the optimization of the truss

    :param file_path: File path to the Excel file to save to.
    :param nodes: Locations of the nodes of the truss. Shape: (# nodes, 2)
    :param members: Indices of the nodes that the members run between. Shape: (# members, 2)
    :param load_cases: Load cases applied to the truss. Shape: (# Load Cases, # Loads per Case, 4: [x,y,fx,fy])
    :param supports: Nodes that are supported and in what directions. Shape: (# Supported nodes, 4: [x,y,sx,sy])
    """

    nodes = np.array(nodes)

    for i, m in enumerate(members):
        members[i].append(((nodes[m[1]][0] - nodes[m[0]][0])**2 + (nodes[m[1]][1] - nodes[m[0]][1])**2)**0.5)
    members = np.array(members)

    vol, A, q, u = optimize_truss_cross_sections(nodes, members, load_cases, supports)
    plt.show()

    save_truss_cross_section_optimization_excel(file_path, nodes, members, load_cases, supports, A, u, q, vol)

def optimize_truss_cross_sections(nodes, members, load_cases, supports):
    """
    Optimizes the cross-sections of a truss.

    :param nodes: Locations of the nodes of the truss. Shape: (# nodes, 2)
    :param members: Indices of the nodes that the members run between. Shape: (# members, 3: [node i, node j, length])
    :param load_cases: Load cases applied to the truss. Shape: (# Load Cases, # Loads per Case, 4: [x,y,fx,fy])
    :param supports: Nodes that are supported and in what directions. Shape: (# Supported nodes, 4: [x,y,sx,sy])
    :return: vol: Volume, a: areas, q: forces, u: deformations
        vol: Overall volume of the truss.
        a: Areas of each member. Shape: (# Members)
        q: Axial forces within each member. Shape: (# Members * # Load Cases)
        u: Deflections at each node. Shape: (2 * # Nodes * # Load Cases)
    """

    f, dof = _assign_loads_and_supports(nodes, load_cases, supports)
    vol, a, q, u = _solve_optimum_problem(nodes, members, f, dof)
    if a is not None:
        # print("Solved")
        _plot_truss(nodes, members, a, max(a) * 1e-3)
        return vol, a, q, u
    else:
        print("Truss Is not stable")
        return None, None, None, None

def _assign_loads_and_supports(nodes, load_cases, supports):
    """
    Makes the vectors for the forces and deflections to be used in stiffness matrices.

    :param nodes: Locations of the nodes of the truss. Shape: (# nodes, 2)
    :param load_cases: Load cases applied to the truss. Shape: (# Load Cases, # Loads per Case, 4: [x,y,fx,fy])
    :param supports: Nodes that are supported and in what directions. Shape: (# Supported nodes, 4: [x,y,sx,sy])
    :return: Vector of forces, Vector of deflections
    """

    total_load = 0
    for case in load_cases:
        for load in case:
            total_load += abs(load[2]) + abs(load[3])
    if total_load == 0:
        total_load = 1

    f = [np.zeros(len(nodes) * 2)] * len(load_cases)
    for j, load_case in enumerate(load_cases):
        for load in load_case:
            for i, nd in enumerate(nodes):
                if nd[0] == load[0] and nd[1] == load[1]:
                    f[j][2*i] = load[2]/total_load
                    f[j][2*i+1] = load[3]/total_load

    dof = np.ones((len(nodes), 2))
    for support in supports:
        for i, nd in enumerate(nodes):
            if nd[0] == support[0] and nd[1] == support[1]:
                if support[2] and support[3]: dof[i, :] = [0, 0]
                elif support[2]: dof[i, :] = [1, 0]
                elif support[3]: dof[i, :] = [0, 1]
                else: dof[i, :] = [1, 1]

    dof = np.array(dof).flatten()

    return f, dof

def _solve_optimum_problem(nodes, members, f, dof):
    """
    Solves for the optimum areas for all the active members

    :param nodes: Locations of the nodes of the truss. Shape: (# nodes, 2)
    :param members: Array of active members: [node 1 i, node 2 i, length]
    :param f: List of forces acting on each node for each load case
    :param dof: List of degrees of freedom for each node
    :return: volume, areas, forces, deformations
    """

    lb = 0.1 / len(f)
    l = members[:, 2]

    # defining variables
    a = cvx.Variable(len(members), name='a', nonneg=True)  # Cross-section areas
    p = [cvx.Variable(len(members), name='p_lc' + str(k)) for k in range(len(f))]  # Axial forces
    q = [cvx.Variable(len(members), name='q_lc' + str(k)) for k in range(len(f))]  # Element elastic energy
    beta = cvx.Variable(len(f), nonneg=True)  # Dual variables of lower bound on case weighting = slacks of eqn

    # setting objective function
    obj = cvx.Minimize((cvx.sum(l @ a) - cvx.sum(lb * beta)))

    # setting equilibrium constraints
    B = _calc_bi(nodes, members, dof)  # reduced
    equilibCon = [B.transpose() @ q[k] + f[k] == 0 for k in range(len(f))]

    # setting conic constraints
    coneConFlat = []
    invLen = [2 / li for li in l]
    for k in range(len(f)):
        xVals = cvx.vstack([2 * q[k], cvx.multiply(invLen, a) - p[k]])
        tVals = cvx.multiply(invLen, a) + p[k]
        coneConFlat.append(cvx.SOC(tVals, xVals))

    # setting load case total energy constraints
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

def _calc_bi(nodes, members, dof=None):
    """
    Calculates the equilibrium matrix. In reduced form if dof information is given.

    :param nodes: Locations of the nodes of the truss. Shape: (# nodes, 2)
    :param members: Indices of the nodes that the members run between. Shape: (# members, 3: [node i, node j, length])
    :param dof: Array degrees of freedom for each node.
    :return: Equilibrium Matrix
    """

    if dof is None: dof = []

    m, n1, n2 = len(members), members[:, 0].astype(int), members[:, 1].astype(int)

    l, dx, dy = members[:, 2], nodes[n2, 0] - nodes[n1, 0], nodes[n2, 1] - nodes[n1, 1]

    # If information on fixed degrees of freedom is provided
    if len(dof):
        d0, d1, d2, d3 = dof[n1 * 2], dof[n1 * 2 + 1], dof[n2 * 2], dof[n2 * 2 + 1]
        v = np.concatenate((-dx / l * d0, -dy / l * d1, dx / l * d2, dy / l * d3))

    else:
        v = np.concatenate((-dx / l, -dy / l, dx / l, dy / l))

    r = np.concatenate((n1 * 2, n1 * 2 + 1, n2 * 2, n2 * 2 + 1))

    c = np.concatenate((np.arange(m), np.arange(m), np.arange(m), np.arange(m)))

    return sparse.coo_matrix((v, (c, r)), shape=(m, len(nodes) * 2))

def _plot_truss(nodes, members, a, threshold):
    """
    Plots the truss.

    :param nodes: Locations of the nodes of the truss. Shape: (# nodes, 2)
    :param members: Indices of the nodes that the members run between. Shape: (# members, 3: [node i, node j, length])
    :param a: Areas of each member. Shape: (# Members)
    :param threshold: Required area of a member for it to be printed.
    """

    plt.figure()
    plt.clf()
    plt.axis('equal')
    plt.draw()

    THICKNESS_MULTIPLIER = 0.08

    for i in [i for i in range(len(a)) if a[i] >= threshold]:
        pos = nodes[members[i, [0, 1]].astype(int), :]
        plt.plot(pos[:, 0], pos[:, 1], linewidth=np.sqrt(a[i] * THICKNESS_MULTIPLIER), solid_capstyle='round')
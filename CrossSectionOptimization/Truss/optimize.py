import numpy as np
import cvxpy as cvx
import matplotlib.pyplot as plt
from scipy import sparse

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

def solveOptimumProblem(Nodes, Members, f, dof):
    """
    Solves for the optimum areas for all the acitve members

    :param Nodes: Array of nodes
    :param Members: Array of active members: [node 1 index, node 2 index, length]
    :param f: list of forces actcting on each node for each load case
    :param dof: list of degress of freedom for each node
    :return: volume, areas, forces, deformations
    """

    lb = 0.1 / len(f)
    l = Members[:, 2]

    # defining variables
    a = cvx.Variable(len(Members), name='a', nonneg=True)  # Cross-section areas
    p = [cvx.Variable(len(Members), name='p_lc' + str(k)) for k in range(len(f))]  # Axial forces
    q = [cvx.Variable(len(Members), name='q_lc' + str(k)) for k in range(len(f))]  # Element elastic energy
    beta = cvx.Variable(len(f), nonneg=True)  # Dual variables of lower bound on case weighting = slacks of eqn

    # setting objective function
    obj = cvx.Minimize((cvx.sum(l @ a) - cvx.sum(lb * beta)))

    # setting equilibrium constraints
    B = calcBi(Nodes, Members, dof)  # reduced
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

def calcBi(Nodes, Members, dof=None):
    """
    Calculates the equilibrium matrix. In reduced form if dof information is given

    :param Nodes: Array of nodes
    :param Members: Array of members ether active or inactive [node 1 index, node 2 index, length]
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

def plotTruss(Nodes, Members, a, threshold):
    """
    Plots the truss

    :param Nodes: Array of nodes
    :param Members: Array of members: [node 1 index, node 2 index, length]
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

        pos = Nodes[Members[i, [0, 1]].astype(int), :]
        plt.plot(pos[:, 0], pos[:, 1], linewidth=np.sqrt(a[i] * ThinknessMult), solid_capstyle='round')


def OptimizeTrussCrossSections(Nodes, Members, loadCasses, supports):
    """
    Optimizes the cross-sections of a truss

    :param Nodes: The locations of the nodes of the truss [[x1,y1],[x2,y2],...]
    :param Members: The nodes that the members run between [[n11,n21,L1],[n12,n22,L2],...]
    :param loadCasses: the load casses applyed to the truss [[[x1,y1,fx1,dy1],[x2,y2,fx2,dy2],...],...]
    :param supports: the nodes the supports are at and what directions are supported [[x1,y1,sx1,sy1],[x2,y2,sx2,sy2],...]
    :return: vol: Volume, a: areas, q: forces, u: deformations
    """
    print(Nodes)
    print(Members)
    print(loadCasses)
    print(supports)

    f, dof = assignLoadsAndSupports(Nodes, loadCasses, supports)
    vol, a, q, u = solveOptimumProblem(Nodes, Members, f, dof)
    if a is not None:
        # print("Solved")
        plotTruss(Nodes, Members, a, max(a) * 1e-3)
        return vol, a, q, u
    else:
        print("Truss Is not stable")
        return None, None, None, None

def main(nodes, members, loadCasses, supports):
    """
        Optimizes the cross-sections of a truss

        :param nodes: The locations of the nodes of the truss [[x1,y1],[x2,y2],...]
        :param members: The nodes that the members run between [[n11,n21],[n12,n22],...]
        :param loadCasses: the load casses applyed to the truss [[[x1,y1,fx1,fy1],[x2,y2,fx2,fy2],...],...]
        :param supports: the nodes the supports are at and what directions are supported [[x1,y1,sx1,sy1],[x2,y2,sx2,sy2],...]
    """

    Nodes = np.array(nodes)

    for i, m in enumerate(members):
        members[i].append(((nodes[m[1]][0] - nodes[m[0]][0])**2 + (nodes[m[1]][1] - nodes[m[0]][1])**2)**0.5)
    Members = np.array(members)

    vol, A, q, u = OptimizeTrussCrossSections(Nodes, Members, loadCasses, supports)
    plt.show()

    print("vol: ", vol)
    print("A: ", A)
    print("q: ", q)
    print("u: ", u)

import numpy as np
import scipy.optimize as opt
from scipy.optimize import NonlinearConstraint

from CrossSectionOptimization.Frame.FrameStiffnessMatrix2D import MemberStiffness, assembleGlobalStiffnessMatrix, getR, \
    solveGlobalStiffnessMatrix

nCand = 0
nMemb = 0
numChecked = 0

def getK0(E, A0, I0, nodes, members):
    nMemb = len(members)

    K0 = [0] * nMemb
    for j in range(nMemb):
        K0[j], _ = MemberStiffness(E, A0, I0,
                                             nodes[int(members[j][0])][0], nodes[int(members[j][0])][1],
                                             nodes[int(members[j][1])][0], nodes[int(members[j][1])][1])
    return K0

a=1
# 11a
def get_w_c(X, p):
    return X ** p

# 11b
def get_w_s(X, q):
    return X ** q

def get_w_d(X, s_eq):
    return X*(1+s_eq)/(1+s_eq*X)

def get_K_DMO(K_master, K0, X, p):
    wc = get_w_c(X, p).tolist()
    for i, row in enumerate(wc):
        for j, ele in enumerate(row):
            wc[i][j] = np.ones((6,6)) * ele

    wc = np.array(wc)
    return np.sum(((wc * (K_master - K0)) + K0), axis=0)

def get_m(X, s_eq, members, crossSections): # 7a
    a = crossSections[:,2] * crossSections[:,1]
    b = np.array([a]*len(members)).T
    return (get_w_d(X, s_eq) * b @ members[:,2]).sum()

def get_C(R, D):  # 7b
    return R.T @ D

def costFunction(m_X0, C_X0, alpha, C_X, m_X): # 9
    return alpha * (C_X/C_X0) + (1-alpha) * (m_X/m_X0)

def objective(X, iterationConstants, constants):


    p, s_eq = iterationConstants
    nodes, members, supports, crossSections, X0, K0, C_X0, K_master, T_master, alpha, R = constants

    X = X.reshape(len(crossSections),len(members))
    print("X: ", X)

    K_DMO = get_K_DMO(K_master, K0, X, p)
    K_global = assembleGlobalStiffnessMatrix(K_DMO, nodes, members)

    D, R = solveGlobalStiffnessMatrix(K_global, R, len(nodes), supports)

    C_X = get_C(R, D)
    m_X = get_m(X, s_eq, members, crossSections)

    m_X0 = get_m(X0, s_eq, members, crossSections)

    #print("Max Deflection: ", D)
    # print((max(D)*10**2))
    # print("Cost Function", costFunction(m_X0, C_X0, alpha, C_X, m_X))
    cost = costFunction(m_X0, C_X0, alpha, C_X, m_X) #+ (max(D)*10**2)

    return cost

#Constraints
def MemberSumConstraint(X):
    global numChecked
    X = X.reshape((nCand, nMemb))
    sumMember = np.max(X, axis=0)[numChecked]
    print(X)
    print(sumMember)
    # sumMember = sum(X[list(range(numChecked * nCand, (numChecked+1) * nCand))])
    if numChecked == nMemb - 1:
         numChecked = 0
    else:
         numChecked += 1
    return sumMember
# sum(X[list(range(i * nCand, (i+1) * nCand +1))])

def assembleConstraints():
    # constraints = []
    # for i in range(nMemb):
    #     constraints.append(opt.NonlinearConstraint(MemberSumConstraint, 0.999, 1.001))

    constraints = []
    for j in range(nMemb):
        def constraint_fun(x, col_idx=j):
            # Reshape x to a matrix internally to access columns
            x_matrix = x.reshape((nCand, nMemb))
            return np.sum(x_matrix[:, col_idx]) - 1.0

        #constraints.append({'type': 'eq', 'fun': constraint_fun})
        constraints.append(NonlinearConstraint(constraint_fun, -0.1, 0.1))

    # for i in range(nCand * nMemb):
    #     def constraint_fun_absalute(x, index = i):
    #         if x[index] == 0 or x[index] == 1:
    #             return 0
    #         else:
    #             return 1
    #     constraints.append({'type': 'eq', 'fun': constraint_fun_absalute})

    return constraints

# # Stress Constraints
# def stressConstrints(X, iterationConstants, constants):
#     p, s_eq = iterationConstants
#     nodes, members, supports, crossSections, X0, K0, C_X0, K_master, T_master, alpha, R = constants
#
#     X = X.reshape(len(crossSections), len(members))
#
#     K_DMO = get_K_DMO(K_master, K0, X, p)
#     K_global = assembleGlobalStiffnessMatrix(K_DMO, nodes, members)
#
#     D, R = solveGlobalStiffnessMatrix(K_global, R, len(nodes), supports)
#
#     F = getForces(members, D, K_DMO, T_master)
# # 10
# def sigma_VM(X, q, sigma_max):
#     return get_w_s(X, q) * sigma_max
# # 12
# def StressConstraint(X, q, sigma_max, sigma_Y):
#     return phi(X, q, sigma_max, sigma_Y) <= 0
# # 13
# def phi(X, q, sigma_max, sigma_Y):
#     return np.sum(sigma_VM(X, q, sigma_max) / sigma_Y - get_w_s(X, q), axis=0)

def FrameMain(nodes, members, loads, supports, crossSections, crossSectionNames, E, A0, I0, alpha):
    """

    :param nodes: [x,y]
    :param members: [n1, n2, L]
    :param loads:
    :param supports:
    :param crossSections: [Name, A, I, rho (density)]
    :param E:
    :param A0:
    :param I0:
    :param alpha:
    :return:
    """
    print("Starting optimization")
    global nCand
    global nMemb

    nCand = len(crossSections)
    nMemb = len(members)

    R = getR(nodes, loads)


    # Calculating all the stiffness matrices for every posible combination of member and cross-section
    K_master = [[0]*nMemb]*nCand
    T_master = [0]*nMemb
    for i in range(nCand):
        for j in range(nMemb):
            K_master[i][j], T_master[j] = MemberStiffness(E, crossSections[i][0], crossSections[i][1], nodes[int(members[j][0])][0], nodes[int(members[j][0])][1], nodes[int(members[j][1])][0], nodes[int(members[j][1])][1])

    # setting inital guess stuff
    x0 = 1/nCand

    # design varable array i: cross-section, j: member x_i,j = {0,1}. Is relaxed during solving to x_i,j = [0,1]
    X0 = np.array([[x0]*nMemb]*nCand)

    K0_DMO = getK0(E, A0, I0, nodes, members)

    K0_global = assembleGlobalStiffnessMatrix(K0_DMO, nodes, members)

    D0, R0 = solveGlobalStiffnessMatrix(K0_global, R, len(nodes), supports)
    C_X0 = get_C(R0, D0)

    constraints = assembleConstraints()

    # runnning iterations

    constants = [nodes, members, supports, crossSections, X0, np.array(K0_DMO), C_X0, K_master, T_master, alpha, R]

    p = 50
    r = 50
    s_eq = 50
    q = 50

    X = X0

    iteration = 0
    while True:
        iteration += 1
        print("Iteration: ", iteration)

        iterationConstants = [p, s_eq]

        X = X.reshape(len(crossSections)* len(members))

        #print("Starting point: ", X)
        OptimizeResult = opt.minimize(objective, X, args=(iterationConstants, constants), method='SLSQP', tol = 0.001, bounds=opt.Bounds([0]*len(X),[1]*len(X)), constraints=constraints)
        #integrality=[True]*len(X)
        #OptimizeResult = opt.differential_evolution(costFunction, bounds=opt.Bounds([0]*len(X),[1]*len(X)), args=(iterationConstants, constants),  constraints=constraints, disp=True)

        X = OptimizeResult.x
        X = X.reshape((len(crossSections), len(members)))
        print("Results: ", np.round(X,2))

        if all(np.max(X, axis=0) > 0.95):
            break
        else:
            # update penaltys
            p += 1
            q += 1.11
            r += 1
            s_eq += 1
            if q > p:
                print("adjust")
                q = p

            # if iteration%20 == 0:
            #     A_terget = np.sum((np.array([crossSections[:,0]]*nMemb).T * X), axis=0)
            #     I_target = np.sum((np.array([crossSections[:,1]]*nMemb).T * X), axis=0)
            #
            #     newCrossSections = []
            #
            #     for i in range(nMemb):
            #         target1 = np.array([[A_terget[i], I_target[i], 0]] * nCand)
            #
            #         mask = (crossSections >= target1).all(axis=1)
            #
            #         if not np.any(mask):
            #             print("No Cross-Sections Large Enough")
            #             newCrossSections.append(0)
            #         else:
            #
            #             valid_rows = crossSections[mask]
            #
            #             target2 = np.array([[A_terget[i], I_target[i], 0]] * len(valid_rows))
            #
            #             distances = np.sum(target2 - valid_rows, axis=1)
            #
            #             closest_valif_row_idx_in_valid_rows = np.argmin(abs(distances))
            #
            #             original_indices = np.where(mask)[0]
            #             original_row_idx = original_indices[closest_valif_row_idx_in_valid_rows]
            #
            #             newCrossSections.append(original_row_idx)
            #
            #     newX = np.zeros((nCand, nMemb))
            #     for i in range(len(newCrossSections)):
            #         newX[newCrossSections[i],i] = 1
            #
            #     print("newX: ", newX)
            #
            #     X = 3*X/4 + newX/4

    # finding what cross-section was selected for each member
    selectedCrossSectionIndex = np.argmax(X, axis=0)
    #print("Selected cross-section index: ", selectedCrossSectionIndex)
    selectedCrossSections = crossSectionNames[selectedCrossSectionIndex]
    print("Cross Sections Slected: ", selectedCrossSections)

    return selectedCrossSections






# Nodes = np.array([[0,0],[1,0],[1,1],[0,1]])
# Members = np.array([[0,3,1],[1,2,1],[2,3,1]])
# Loads = np.array([[2,0,-1,0],[3,0,-1,0]])
# Supports = np.array([[0,1,1,0],[1,1,1,0]])
# CrossSections = np.array([[5,3,1],[10,20,1 ]])
# crossSectionNames = np.array(["C1", "C2"])
# E_set = 200000
#
# A_0 = 0.01
# I_0 = 0.01
# Alpha = 0.2
#
# selectedSections = main(Nodes, Members, Loads, Supports, CrossSections, crossSectionNames, E_set, A_0, I_0, Alpha)
# print("Selected Sections: ", selectedSections)












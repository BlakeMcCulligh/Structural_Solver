import time
import numpy as np
from numpy.linalg import LinAlgError
import scipy.optimize as opt
from scipy.optimize import NonlinearConstraint

from CrossSectionOptimization.Frame.FrameStiffnessMatrix2D import solveGlobalStiffnessMatrix

numSingular = 0
numRun = 0

class UnderDefinedStructureError(Exception):
    pass

def getK_local(nMemb, E, A, I, L):
    K_local = []
    for i in range(nMemb):
        AEL = A[i] * E / L[i]
        EI = E * I[i]
        k_local = np.array([
            [AEL, 0, 0, -AEL, 0, 0],
            [0, 12 * EI / L[i] ** 3, 6 * EI / L[i] ** 2, 0, -12 * EI / L[i] ** 3, 6 * EI / L[i] ** 2],
            [0, 6 * EI / L[i] ** 2, 4 * EI / L[i], 0, -6 * EI / L[i] ** 2, 2 * EI / L[i]],
            [-AEL, 0, 0, AEL, 0, 0],
            [0, -12 * EI / L[i] ** 3, -6 * EI / L[i] ** 2, 0, 12 * EI / L[i] ** 3, -6 * EI / L[i] ** 2],
            [0, 6 * EI / L[i] ** 2, 2 * EI / L[i], 0, -6 * EI / L[i] ** 2, 4 * EI / L[i]]
        ])
        K_local.append(k_local)
    return K_local

def getT(nMemb, C, S):
    T = []
    for i in range(nMemb):
        c, s = C[i], S[i]
        T.append(np.array([
                [c, s, 0, 0, 0, 0],
                [-s, c, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, c, s, 0],
                [0, 0, 0, -s, c, 0],
                [0, 0, 0, 0, 0, 1]
            ]))
    return T

def assembaleGlobalStiffnessMatrix(N, M, T, K_local_e):
    K_global = np.zeros((len(N) * 3, len(N) * 3)).tolist()
    for index, m in enumerate(M):
        k_global_e = T[index].T @ K_local_e[index] @ T[index]
        global_dofs = [3 * m[0], 3 * m[0] + 1, 3 * m[0] + 2, 3 * m[1], 3 * m[1] + 1, 3 * m[1] + 2]
        for i in range(6):
            for j in range(6):
                K_global[int(global_dofs[i])][int(global_dofs[j])] += k_global_e[i, j]
    return np.array(K_global)

def getF(N, Loads):
    F  = np.zeros(len(N) * 3)
    for load in Loads:
        F[load[0]*3] = load[1]
        F[load[0]*3 + 1] = load[2]
    return F

def createBounds(nMemb, rangeA, rangeI):
    b = opt.Bounds([rangeA[0]]*nMemb + [rangeI[0]]*nMemb, [rangeA[1]]*nMemb + [rangeI[1]]*nMemb)
    return b

def createIConstraint(nMemb, AlowableMultiplyer):
    constI = []
    for i in range(nMemb):
        def IConst(X,indexA=i,indexI = i+nMemb):
            return X[indexI] / X[indexA]
        constI.append(NonlinearConstraint(IConst, AlowableMultiplyer[0], AlowableMultiplyer[1]))
    return constI

def calculateNumFailed(M, T, U_global, K, A, I, depth, Fy):
    U_global_e = []
    for i, m in enumerate(M):
        U_global_e.append([U_global[m[0]*3], U_global[m[0]*3+1], U_global[m[0]*3+2], U_global[m[1]*3], U_global[m[1]*3+1], U_global[m[1]*3+2]])
    Sigma = []
    for i in range(len(M)):
        U_local = T[i] @ U_global_e[i]
        F_members = K[i] @ U_local
        sigma_Axial = F_members[[0, 3]] / A[i]
        sigma_bending = F_members[[1, 4]] * depth[i] / I[i]
        sigma = sigma_Axial + sigma_bending
        Sigma.append(sigma)
    numFailed = np.sum(np.array(Sigma) > Fy)
    print("Number Failed: ", numFailed)
    return numFailed


def objectiveFunction(X, Constants):
    E, N, M, F, S, L, T, Fy, depthMult = Constants
    X = X.reshape(len(M), 2)
    A = X[:, 0]
    I = X[:, 1]
    global numRun
    numRun += 1
    nMemb = len(M)
    K_local = getK_local(nMemb, E, A, I, L)
    K = assembaleGlobalStiffnessMatrix(N, M, T, K_local)
    depth = I / A * depthMult
    try:
        U, F = solveGlobalStiffnessMatrix(K, F, len(N), S)
        numFailed = calculateNumFailed(M, T, U, K_local, A, I, depth, Fy)
        cost = (sum(A * L) + max(U)) * (numFailed + 1)
    except LinAlgError:
        print(" Sigular Matrix")
        cost = 9.99e15
        global numSingular
        numSingular += 1
        if numRun == numSingular:
            raise UnderDefinedStructureError("Structure is unstable, add more supports or remove member releces")
    print(cost)
    return cost

def optimize(E, N, M, Loads, S, Fy, depthMult):
    L = ((N[M[:,1]][:,0] - N[M[:,0]][:,0])**2 + (N[M[:,1]][:,1] - N[M[:,0]][:,1])**2)**0.5
    c = (N[M[:, 1]][:, 0] - N[M[:, 0]][:, 0]) / L
    s = (N[M[:, 1]][:, 1] - N[M[:, 0]][:, 1]) / L
    F = getF(N, Loads)
    T = getT(len(M), c, s)
    bounds = createBounds(len(M), [0, 100], [100, 10000])
    constraintI = createIConstraint(len(M), [2,50])
    Constants = [E, N, M, F, S, L, T, Fy, depthMult]
    startValue = [10, 200]
    X = np.array([startValue[0]] * len(M) + [startValue[1]] * len(M))
    print("Optimizing")
    startTime = time.time()
    OptimizeResult = opt.minimize(objectiveFunction, X, args=(Constants), method='SLSQP', tol=0.001, bounds = bounds, constraints=constraintI)
    print("Elaps Time: ", time.time()-startTime)
    # best time 0.012356996536254883
    print("Results: ", OptimizeResult.x)

# testing
nodes = np.array([[0,0],[5,0],[5,5],[0,5]])
members = np.array([[0,3],[3,2],[2,1]])
loads = np.array([[2,0,-10,0]])
supports = np.array([[0,1,1,0],[1,1,1,0]])
E_set = 200000
Fy_set = 1
DepthMult = 1000
optimize(E_set, nodes, members, loads, supports, Fy_set, DepthMult)

# TODO
# 1. add stress constraint system : DONE
# 2. add distributed loads
# 3. add releces
# 4. add A and I constraint system based on cross-section table, also add depth function
# 5. add member selection based on results
# 6. add load casses
# 7. add disply
# 8. add reading in and out
# 9. add member groups and difrenet cross-section tables for each



# ---------------------------------
# OLD CODE THAT IS NOT USED ANYMORE
# ---------------------------------

# def MemberStiffness(E, A, I, L, c, s):
#     #https://github.com/anastruct/anaStruct/blob/master/anastruct/fem/elements.py line 323 for how to do hinges
#     AEL = A*E/L
#     EI = E*I
#
#     k_local = np.array([
#         [ AEL,           0,          0, -AEL,           0,          0],
#         [   0,  12*EI/L**3,  6*EI/L**2,    0, -12*EI/L**3,  6*EI/L**2],
#         [   0,   6*EI/L**2,     4*EI/L,    0,  -6*EI/L**2,     2*EI/L],
#         [-AEL,           0,          0,  AEL,           0,          0],
#         [   0, -12*EI/L**3, -6*EI/L**2,    0,  12*EI/L**3, -6*EI/L**2],
#         [   0,   6*EI/L**2,     2*EI/L,    0,  -6*EI/L**2,     4*EI/L]
#     ])
#
#     T = np.array([
#         [ c, s, 0,  0, 0, 0],
#         [-s, c, 0,  0, 0, 0],
#         [ 0, 0, 1,  0, 0, 0],
#         [ 0, 0, 0,  c, s, 0],
#         [ 0, 0, 0, -s, c, 0],
#         [ 0, 0, 0,  0, 0, 1]
#     ])
#
#     return T.T @ k_local @ T

# def getStiffnessMatrix(E, N, M, A, I, L, c, s):
#     K_global = np.zeros((len(N)*3, len(N)*3)).tolist()
#
#     for index, m in enumerate(M):
#
#         K = MemberStiffness(E, A[index], I[index], L[index], c[index], s[index])
#
#         global_dofs = [3 * m[0], 3 * m[0] + 1, 3 * m[0] + 2, 3 * m[1], 3 * m[1] + 1, 3 * m[1] + 2]
#
#         for i in range(6):
#             for j in range(6):
#                 K_global[int(global_dofs[i])][int(global_dofs[j])] += K[i, j]
#
#     return np.array(K_global)
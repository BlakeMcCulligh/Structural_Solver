import numpy as np
import cvxpy as cvx
import matplotlib.pyplot as plt
from scipy import sparse
import sympy as sp

from CrossSectionOptimization.Frame.FrameStiffnessMatrix2D import solveGlobalStiffnessMatrix


def MemberStiffness(E, A, I, x1, y1, x2, y2):
    #https://github.com/anastruct/anaStruct/blob/master/anastruct/fem/elements.py line 323 for how to do hinges
    L = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
    c = (x2 - x1) / L
    s = (y2 - y1) / L

    k_local = np.array([
        [ A*E/L,            0,           0, -A*E/L,            0,           0],
        [     0,  12*E*I/L**3,  6*E*I/L**2,      0, -12*E*I/L**3,  6*E*I/L**2],
        [     0,   6*E*I/L**2,     4*E*I/L,      0,  -6*E*I/L**2,     2*E*I/L],
        [-A*E/L,            0,           0,  A*E/L,            0,           0],
        [     0, -12*E*I/L**3, -6*E*I/L**2,      0,  12*E*I/L**3, -6*E*I/L**2],
        [     0,   6*E*I/L**2,     2*E*I/L,      0,  -6*E*I/L**2,     4*E*I/L]
    ])

    T = np.array([
        [ c, s, 0,  0, 0, 0],
        [-s, c, 0,  0, 0, 0],
        [ 0, 0, 1,  0, 0, 0],
        [ 0, 0, 0,  c, s, 0],
        [ 0, 0, 0, -s, c, 0],
        [ 0, 0, 0,  0, 0, 1]
    ])

    return T.T @ k_local @ T

def getStiffnessMatrix(E, N, M, A, I):
    K_global = np.zeros((len(N)*3, len(N)*3)).tolist()

    for index, m in enumerate(M):
        x1, y1, x2, y2 = N[m[0]][0], N[m[0]][1], N[m[1]][0], N[m[1]][1]

        #K = MemberStiffness(E, A[index], I[index], x1, y1, x2, y2)
        K = MemberStiffness(E, A, I, x1, y1, x2, y2)

        dofs_i = [3 * m[0], 3 * m[0] + 1, 3 * m[0] + 2]
        dofs_j = [3 * m[1], 3 * m[1] + 1, 3 * m[1] + 2]
        global_dofs = dofs_i + dofs_j

        for i in range(6):
            for j in range(6):
                global_i = int(global_dofs[i])
                global_j = int(global_dofs[j])
                K_global[global_i][global_j] += K[i, j]

    return np.array(K_global)

def getF(N, Loads):

    F  = np.zeros(len(N) * 3)

    for load in Loads:
        F[load[0]*3] = load[1]
        F[load[0]*3 + 1] = load[2]
    return F


def optimize(E, N, M, Loads, S):

    L = ((N[M[:,1]][:,0] - N[M[:,0]][:,0])**2 + (N[M[:,1]][:,1] - N[M[:,0]][:,1])**2)**0.5

    F = getF(N, Loads)

    #A = cvx.Variable(len(M), name='A', nonneg=True)  # Cross-section areas
    #I = cvx .Variable(len(M), name='I', nonneg=True) # Moment of Inertia
    #U = cvx .Variable(len(N)*3, name='U', nonneg=True) # deformation vector

    A, I = sp.symbols('A, I')
    K = getStiffnessMatrix(E, N, M, A, I)
    print(solveGlobalStiffnessMatrix(K, F, len(N), S))

    # print(K[0][0])
    #
    # obj = cvx.Minimize(cvx.max(U)+cvx.sum(A@L))
    #
    # stiffnessCons = [F == K @ U]
    #
    # prob = cvx.Problem(obj, stiffnessCons)
    # volume = prob.solve(cvx.MOSEK, verbose=False)
    # print(volume)




# testing
nodes = np.array([[0,0],[5,0],[5,5],[0,5]])
members = np.array([[0,3],[3,2],[2,1]])
loads = np.array([[2,0,-1,0]])
supports = np.array([[0,1,1,0],[1,1,1,0]])
E_set = 200000
#print("Matrix: ", getStiffnessMatrix(1, nodes, members, 1, 1))
optimize(E_set, nodes, members, loads, supports)








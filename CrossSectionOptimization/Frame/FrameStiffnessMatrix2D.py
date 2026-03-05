import numpy as np


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

    return T.T @ k_local @ T, T

def assembleGlobalStiffnessMatrix(K_DMO, nodes, members):
    total_dof = len(nodes) * 3
    K_global = np.zeros((total_dof, total_dof)).tolist()

    for i, m in enumerate(members):

        dofs_i = [3 * m[0], 3 * m[0] + 1, 3 * m[0] + 2]
        dofs_j = [3 * m[1], 3 * m[1] + 1, 3 * m[1] + 2]
        global_dofs = dofs_i + dofs_j

        for j in range(6):
            for k in range(6):
                global_i = int(global_dofs[j])
                global_j = int(global_dofs[k])
                K_global[global_i, global_j] += K_DMO[i][j, k]

    return K_global

def getR(nodes, loads):

    R  = np.zeros(len(nodes) * 3)

    for load in loads:
        R[load[0]*3] = load[1]
        R[load[0]*3 + 1] = load[2]
    return R

def solveGlobalStiffnessMatrix(K, R, nNodes, supports):
    fixed_dofs = []
    for support in supports:
        if support[1] == 1:
            fixed_dofs.append(support[0]*3)
        if support[2] == 1:
            fixed_dofs.append(support[0] * 3 + 1)
        if support[3] == 1:
            fixed_dofs.append(support[0] * 3 + 2)



    free_dofs = list(set(range(nNodes*3)) - set(fixed_dofs))
    Kff = K[np.ix_(free_dofs, free_dofs)]
    Df = np.linalg.solve(Kff, R[free_dofs])
    D = np.zeros(nNodes*3)
    D[free_dofs] = Df

    R[fixed_dofs] = K[fixed_dofs,:]@D

    return D, R

def getForces(members, D, K_DMO, T_master):

    F = []
    for i, m in enumerate(members):
        D_global = [D[m[0]*3], D[m[0]*3+1], D[m[0]*3+2], D[m[1]*3], D[m[1]*3+1], D[m[1]*3+2]]

        D_local = T_master[i] @ D_global

        F.append( ( T_master[i] @ K_DMO[i] @ T_master[i].T ) @ D_local )

    return F





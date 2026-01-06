import numpy as np

def ElementStiffnessMatrixLocal(E, A, L):
    K_mPrime = np.zeros(6)

    K_mPrime[0][0] = 1
    K_mPrime[3][0] = -1
    K_mPrime[0][3] = -1
    K_mPrime[3][3] = 1

    mult = E*A/L
    K_mPrime = mult * K_mPrime

    return K_mPrime

def TransformationMatrox(L, N1, N2):

    Cx = (N2[0] - N1[0]) / L
    Cy = (N2[1] - N1[1]) / L
    Cz = (N2[2] - N1[2]) / L

    T_m =  np.zeros(6)

    T_m[0][0] = Cx
    T_m[3][3] = Cx

    T_m[0][1] = Cy
    T_m[3][4] = Cy

    T_m[0][2] = Cz
    T_m[3][5] = Cz

    return T_m

def EleStiffnessMatrixLocalToGlobal(K_mPrime, T_m):
    T_m_T = T_m.transpose()
    K_m = T_m_T @ K_mPrime @ T_m
    return K_m

def AssembleGlobalStiffnessMatrix(kel, elements, Nnodes):
    """
    Creates the global stiffness matrix

    :param kel: elements stiffness matrices
    :param elements: elements objects
    :param Nnodes: number of nodes
    :return: Global stiffness matrix
    """

    K = np.zeros(3*Nnodes)

    for i in range(len(elements)):

        eldofs = np.concatenate([np.arange(3 * elements[i].N1, 3 * (elements[i].N1 + 1)), np.arange(3 * elements[i].N2, 3 * (elements[i].N2 + 1))])

        K[np.ix_(eldofs, eldofs)] += kel[i]

    return K

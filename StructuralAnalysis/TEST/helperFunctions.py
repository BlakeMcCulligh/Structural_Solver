
import numpy as np
import math

from StructuralAnalysis.TEST.__main__ import Frame3D_T

def partD(model: Frame3D_T):

    D_unknown = []
    D_known = []
    D_known_val = []

    for i in range(len(model.nodes_cord)):
        sup = model.nodes_Support[i]
        enf = model.nodes_D_Enforced[i]

        for direc in range(6):
            if sup[direc] == False and enf[direc] is None:
                D_unknown.append(i * 6 + direc)
            elif enf[direc] is not None:
                D_known.append(i * 6 + direc)
                D_known_val.append(enf[direc])
            else:
                D_known.append(i * 6 + direc)
                D_known_val.append(0.0)

    D_known_val = np.array(D_known_val, ndmin=2).T

    return D_unknown, D_known, D_known_val

def prepMembers(model: Frame3D_T):
    DOFs = []
    Ls = []
    memberPartD = []
    T = []
    K = []
    for i in range(len(model.members)):
        DOFs.append(buildDOFVector([model.members[0], model.members[1]]))
        Ls.append(np.linalg.norm(model.nodes_cord[model.members[0]], model.nodes_cord[model.members[1]]))
        memberPartD.append(member_PartD(model.members_Releases[i]))
        T.append(getMemberT(model,i,Ls))
        if model.members[3]:
            K.append(getMeberK(model, i, T, Ls))
        else:
            K.append(None)

    return DOFs, Ls, memberPartD, T, K

def buildDOFVector(listNodes_INDEX):
    dofs = np.empty(len(listNodes_INDEX) * 6, dtype=np.int64)
    local = np.arange(6, dtype=np.int64)
    for i in listNodes_INDEX:
        start = i * 6
        dofs[start:start + 6] = i * 6 + local
    return dofs

def member_PartD(member_Releases):
    R_unrelesed = []
    R_relesed = []
    for i in range(12):
        if not member_Releases[i]:
            R_unrelesed.append(i)
        else:
            R_relesed.append(i)
    return R_unrelesed, R_relesed

def getMemberT(model: Frame3D_T, memberIndex, Ls):

    Xi = model.nodes_cord[model.members[memberIndex][0]][0]
    Xj = model.nodes_cord[model.members[memberIndex][1]][0]
    Yi = model.nodes_cord[model.members[memberIndex][0]][1]
    Yj = model.nodes_cord[model.members[memberIndex][1]][1]
    Zi = model.nodes_cord[model.members[memberIndex][0]][2]
    Zj = model.nodes_cord[model.members[memberIndex][1]][2]

    L = Ls[memberIndex]

    # Calculate the direction cosines for the local x-axis
    x = [(Xj - Xi) / L, (Yj - Yi) / L, (Zj - Zi) / L]

    # Vertical members
    if math.isclose(Xi, Xj) and math.isclose(Zi, Zj):
        if Yj > Yi:
            y = [-1, 0, 0]
            z = [0, 0, 1]
        else:
            y = [1, 0, 0]
            z = [0, 0, 1]

    # Horizontal members
    elif math.isclose(Yi, Yj):
        y = [0, 1, 0]
        z = np.cross(x, y)
        z = np.divide(z, (z[0] ** 2 + z[1] ** 2 + z[2] ** 2) ** 0.5)

    # Members neither vertical nor horizontal
    else:
        proj = [Xj - Xi, 0, Zj - Zi]

        if Yj > Yi:
            z = np.cross(proj, x)
        else:
            z = np.cross(x, proj)

        z = np.divide(z, (z[0] ** 2 + z[1] ** 2 + z[2] ** 2) ** 0.5)
        y = np.cross(z, x)
        y = np.divide(y, (y[0] ** 2 + y[1] ** 2 + y[2] ** 2) ** 0.5)

    # Create the direction cosines matrix
    dirCos = np.array([x, y, z])

    # Build the transformation matrix
    transMatrix = np.zeros((12, 12))
    transMatrix[0:3, 0:3] = dirCos
    transMatrix[3:6, 3:6] = dirCos
    transMatrix[6:9, 6:9] = dirCos
    transMatrix[9:12, 9:12] = dirCos

    return transMatrix

def getMeberK(model: Frame3D_T, memberIndex, T, Ls):
    return np.matmul(np.matmul(np.linalg.inv(T), getMemberk(model, memberIndex, Ls)), T)

def getMemberk(model: Frame3D_T, memberIndex, Ls):
    E = model.materials[model.members[3]][0]
    G = model.materials[model.members[3]][1]
    Iy = model.members_CrossSectionProps[1]
    Iz = model.members_CrossSectionProps[2]
    J = model.members_CrossSectionProps[3]
    A = model.members_CrossSectionProps[0]

    L = Ls[memberIndex]

    # Create the uncondensed local stiffness matrix
    k = np.array([[A * E / L, 0, 0, 0, 0, 0, -A * E / L, 0, 0, 0, 0, 0],
               [0, 12 * E * Iz / L ** 3, 0, 0, 0, 6 * E * Iz / L ** 2, 0, -12 * E * Iz / L ** 3, 0, 0, 0,
                6 * E * Iz / L ** 2],
               [0, 0, 12 * E * Iy / L ** 3, 0, -6 * E * Iy / L ** 2, 0, 0, 0, -12 * E * Iy / L ** 3, 0,
                -6 * E * Iy / L ** 2, 0],
               [0, 0, 0, G * J / L, 0, 0, 0, 0, 0, -G * J / L, 0, 0],
               [0, 0, -6 * E * Iy / L ** 2, 0, 4 * E * Iy / L, 0, 0, 0, 6 * E * Iy / L ** 2, 0, 2 * E * Iy / L, 0],
               [0, 6 * E * Iz / L ** 2, 0, 0, 0, 4 * E * Iz / L, 0, -6 * E * Iz / L ** 2, 0, 0, 0, 2 * E * Iz / L],
               [-A * E / L, 0, 0, 0, 0, 0, A * E / L, 0, 0, 0, 0, 0],
               [0, -12 * E * Iz / L ** 3, 0, 0, 0, -6 * E * Iz / L ** 2, 0, 12 * E * Iz / L ** 3, 0, 0, 0,
                -6 * E * Iz / L ** 2],
               [0, 0, -12 * E * Iy / L ** 3, 0, 6 * E * Iy / L ** 2, 0, 0, 0, 12 * E * Iy / L ** 3, 0,
                6 * E * Iy / L ** 2, 0],
               [0, 0, 0, -G * J / L, 0, 0, 0, 0, 0, G * J / L, 0, 0],
               [0, 0, -6 * E * Iy / L ** 2, 0, 2 * E * Iy / L, 0, 0, 0, 6 * E * Iy / L ** 2, 0, 4 * E * Iy / L, 0],
               [0, 6 * E * Iz / L ** 2, 0, 0, 0, 2 * E * Iz / L, 0, -6 * E * Iz / L ** 2, 0, 0, 0, 4 * E * Iz / L]])

    k11, k12, k21, k22 = memberPart(model, k, memberIndex)
    k = np.subtract(k11, np.matmul(np.matmul(k12, np.linalg.inv(k22)), k21))
    i = 0
    for DOF in model.members_Releases[memberIndex]:
        if DOF:
            k = np.insert(k, i, 0, axis=0)
            k = np.insert(k, i, 0, axis=1)
        i += 1
    return k

def memberPart(model: Frame3D_T, matrix, i):
    R_unrelesed, R_relesed = member_PartD(model.members_Releases[i])
    if matrix.shape[1] == 1:
        m1 = matrix[R_unrelesed, :]
        m2 = matrix[R_relesed, :]
        return m1, m2
    else:
        m11 = matrix[R_unrelesed, :][:, R_unrelesed]
        m12 = matrix[R_unrelesed, :][:, R_relesed]
        m21 = matrix[R_relesed, :][:, R_unrelesed]
        m22 = matrix[R_relesed, :][:, R_relesed]
        return m11, m12, m21, m22

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
    model.members = np.array(model.members)
    model.nodes_cord = np.array(model.nodes_cord)


    Ls = get_L_ARRAY(model)

    DOFs = []
    memberPartD = []
    T = []
    for i in range(len(model.members)):
        DOFs.append(buildDOFVector([model.members[i,0], model.members[i,1]]))
        memberPartD.append(member_PartD(model.members_Releases[i]))
        T.append(getMemberT(model,i,Ls))
    return DOFs, Ls, memberPartD, T

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

# def getGlobalFixedEndReactionVector(model: Frame3D_T):
#
#     # Get the partitioned global fixed end reaction vector
#     FER_val = np.zeros((len(model.nodes_cord) * 6, 1))
#
#
#     for phys_member in self.members.values():
#         for member in phys_member.sub_members.values():
#
#             member_FER = np.asarray(member.FER(combo.name), dtype=float).reshape(-1)
#             FER_val[member.dofs, 0] += member_FER
#     FER1, FER2 = partition(FER_val, D1_indices, D2_indices)
#
#     return FER1, FER2

def get_L_ARRAY(model: Frame3D_T):
    dx = model.nodes_cord[model.members[:, 1], 0] - model.nodes_cord[model.members[:, 0], 0]
    dy = model.nodes_cord[model.members[:, 1], 1] - model.nodes_cord[model.members[:, 0], 1]
    dz = model.nodes_cord[model.members[:, 1], 2] - model.nodes_cord[model.members[:, 0], 2]
    return np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

def get_k_local_relesed_ARRAY(model: Frame3D_T, k_local_array, k11, k12, k21, k22):
    k_ARRAY = []
    for i in range(len(k_local_array)):
        k_ARRAY.append(np.subtract(k11, np.matmul(np.matmul(k12, np.linalg.inv(k22)), k21)))
        j = 0
        for DOF in model.members_Releases[i]:
            if DOF:
                k_ARRAY[i] = np.insert(k_ARRAY[i], j, 0, axis=0)
                k_ARRAY[i] = np.insert(k_ARRAY[i], j, 0, axis=1)
            j += 1
    return k_ARRAY

def get_k_local_ARRAY(model: Frame3D_T, Ls):

    E = model.materials[model.members[:, 2], 0]
    G = model.materials[model.members[:, 2], 1]

    A = model.members_CrossSectionProps[:, 0]
    Iy = model.members_CrossSectionProps[:, 1]
    Iz = model.members_CrossSectionProps[:, 2]
    J = model.members_CrossSectionProps[:, 3]

    L = Ls

    # Create the uncondensed local stiffness matrix
    k_local_array = []
    for i in range(len(model.members)):
        k_local_array.append(np.array([[A[i] * E[i] / L[i], 0, 0, 0, 0, 0, -A[i] * E[i] / L[i], 0, 0, 0, 0, 0],
               [0, 12 * E[i] * Iz[i] / L[i] ** 3, 0, 0, 0, 6 * E[i] * Iz[i] / L[i] ** 2, 0, -12 * E[i] * Iz[i] / L[i] ** 3, 0, 0, 0,
                6 * E[i] * Iz[i] / L[i] ** 2],
               [0, 0, 12 * E[i] * Iy[i] / L[i] ** 3, 0, -6 * E[i] * Iy[i] / L[i] ** 2, 0, 0, 0, -12 * E[i] * Iy[i] / L[i] ** 3, 0,
                -6 * E[i] * Iy[i] / L[i] ** 2, 0],
               [0, 0, 0, G[i] * J[i] / L[i], 0, 0, 0, 0, 0, -G[i] * J[i] / L[i], 0, 0],
               [0, 0, -6 * E[i] * Iy[i] / L[i] ** 2, 0, 4 * E[i] * Iy[i] / L[i], 0, 0, 0, 6 * E[i] * Iy[i] / L[i] ** 2, 0, 2 * E[i] * Iy[i] / L[i], 0],
               [0, 6 * E[i] * Iz[i] / L[i] ** 2, 0, 0, 0, 4 * E[i] * Iz[i] / L[i], 0, -6 * E[i] * Iz[i] / L[i] ** 2, 0, 0, 0, 2 * E[i] * Iz[i] / L[i]],
               [-A[i] * E[i] / L[i], 0, 0, 0, 0, 0, A[i] * E[i] / L[i], 0, 0, 0, 0, 0],
               [0, -12 * E[i] * Iz[i] / L[i] ** 3, 0, 0, 0, -6 * E[i] * Iz[i] / L[i] ** 2, 0, 12 * E[i] * Iz[i] / L[i] ** 3, 0, 0, 0,
                -6 * E[i] * Iz[i] / L[i] ** 2],
               [0, 0, -12 * E[i] * Iy[i] / L[i] ** 3, 0, 6 * E[i] * Iy[i] / L[i] ** 2, 0, 0, 0, 12 * E[i] * Iy[i] / L[i] ** 3, 0,
                6 * E[i] * Iy[i] / L[i] ** 2, 0],
               [0, 0, 0, -G[i] * J[i] / L[i], 0, 0, 0, 0, 0, G[i] * J[i] / L[i], 0, 0],
               [0, 0, -6 * E[i] * Iy[i] / L[i] ** 2, 0, 2 * E[i] * Iy[i] / L[i], 0, 0, 0, 6 * E[i] * Iy[i] / L[i] ** 2, 0, 4 * E[i] * Iy[i] / L[i], 0],
               [0, 6 * E[i] * Iz[i] / L[i] ** 2, 0, 0, 0, 2 * E[i] * Iz[i] / L[i], 0, -6 * E[i] * Iz[i] / L[i] ** 2, 0, 0, 0, 4 * E[i] * Iz[i] / L[i]]]))

    return k_local_array

def memberPart_k_ARRAY(k, R_unrelesed, R_relesed):
    k1 = k[:,R_unrelesed, :]
    k2 = k[:,R_relesed, :]

    k11 = k1[:, :, R_unrelesed]
    k12 = k1[:, :, R_relesed]
    k21 = k2[:, :, R_unrelesed]
    k22 = k2[:, :, R_relesed]
    return k11, k12, k21, k22

def get_FER_ARRAY(fer, members_T, numM, numC):
    FER_array = []
    for i in range(numM):
        invT = np.linalg.inv(members_T[i])
        sub_FER_array = []
        for j in range(numC):
            sub_FER_array.append(np.matmul(invT, fer[i, j]))
        FER_array.append(sub_FER_array)
    return FER_array

def get_fer_ARRAY(model, k11, k12, k21, k22, fer1, fer2):
    ferCondensed_ARRAY = []
    for i in range(len(k12)):
        ferCondensed = np.subtract(fer1[i], np.matmul(np.matmul(k12[i], np.linalg.inv(k22[i])), fer2[i]))
        i = 0
        for DOF in model.members_Releases:
            if DOF:
                ferCondensed = np.insert(ferCondensed, i, 0, axis=0)
            i += 1
        ferCondensed_ARRAY.append(ferCondensed)
    return ferCondensed_ARRAY

def memberPart_fer_ARRAY(fer, R_unrelesed, R_relesed):
    fer1 = fer[:, R_unrelesed, :]
    fer2 = fer[:, R_relesed, :]
    return fer1, fer2

def get_fer_unc_ARRAY(model, pointLoads, distLoads):

    # pointLoads [memberINDEX, caseINDEX, [x, Px, Py, Pz, Mx, My, Mz]]
    # distLoads [memberINDEX, caseINDEX, [x1, x2, wz1, wz2, wy1, wy2]]

    fer_unc_ARRAY = np.zeros((len(model.members), 12, 1))
    fer_unc_ARRAY = np.add(fer_unc_ARRAY, get_FixedEndReactions_Pointload_ARRAY(model, pointLoads))
    fer_unc_ARRAY = np.add(fer_unc_ARRAY, get_FixedEndReactions_Distload_ARRAY(model, distLoads))

    return fer_unc_ARRAY


def get_FixedEndReactions_Pointload_ARRAY(model, pointLoads):
    # pointLoads [memberINDEX, caseINDEX, [x, Px, Py, Pz, Mx, My, Mz]] # TODO
    pass

def get_FixedEndReactions_Distload_ARRAY(model, distLoads):
    # distLoads [memberINDEX, caseINDEX, [x1, x2, wz1, wz2, wy1, wy2]] # TODO
    pass

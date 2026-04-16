import numpy as np
import math
import StructuralAnalysis.TEST.fixedEndReactionsCalculaters as ferCalc

from StructuralAnalysis.TEST.__main__ import Frame3D_T

def partD(model: Frame3D_T):
    D_unknown = []
    D_known = []
    for i in range(len(model.nodes_cord)):
        sup = model.nodes_Support[i]

        for direc in range(6):
            if not sup[direc]:
                D_unknown.append(i * 6 + direc)
            else:
                D_known.append(i * 6 + direc)
    return D_unknown, D_known

def prepMembers(model: Frame3D_T):
    model.members = np.array(model.members)
    model.nodes_cord = np.array(model.nodes_cord)
    Ls = get_L_ARRAY(model)
    DOFs = []
    memberPartD_unrelesed = []
    memberPartD_relesed = []
    T = []
    for i in range(len(model.members)):
        DOFs.append(buildDOFVector([model.members[i,0], model.members[i,1]]))
        p1, p2 = member_PartD(model.members_Releases[i])
        memberPartD_unrelesed.append(p1)
        memberPartD_relesed.append(p2)
        T.append(getMemberT(model,i,Ls))
    return DOFs, Ls, memberPartD_unrelesed, memberPartD_relesed, T

def buildDOFVector(listNodes_INDEX):
    dofs = np.empty(len(listNodes_INDEX) * 6, dtype=np.int64)
    local = np.arange(6, dtype=np.int64)
    for i, node in enumerate(listNodes_INDEX):
        start = i * 6
        dofs[start:start + 6] = node * 6 + local
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

def getGlobalFixedEndReactionVector_ARRAY(model: Frame3D_T, pointLoads, distLoads, R_unrelesed, R_relesed, k12, k22, members_T, D_unknown, D_known):
    numM = len(model.members)
    numC = len(model.casses)
    fer_unc_ARRAY = get_fer_unc_ARRAY(model, pointLoads, distLoads)
    fer1, fer2 = memberPart_fer_ARRAY(fer_unc_ARRAY, R_unrelesed, R_relesed)
    ferCondensed_ARRAY = get_fer_ARRAY(model, k12, k22, fer1, fer2, numM, numC)
    FER_array = get_FER_ARRAY(ferCondensed_ARRAY, members_T, numM, numC)
    FER_val = [np.zeros((len(model.nodes_cord) * 6, 1))] * numC
    for i in range(numM):
        for j in range(numC):
            FER = FER_array[i][j]
            member_FER = np.asarray(FER, dtype=float).reshape(-1)
            FER_val[j][model.members_DOF[i], 0] += member_FER
    FER1, FER2 = partFER(FER_val, D_unknown, D_known, numC)
    return FER1, FER2

def partFER(FER_val, D1_indices, D2_indices, numC):
    FER1 = []
    FER2 = []
    for j in range(numC):
        FER_val_sub = FER_val[j]
        FER1.append(FER_val_sub[D1_indices, :])
        FER2.append(FER_val_sub[D2_indices, :])
    return FER1, FER2

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
    return np.array(k_local_array)

def memberPart_k_ARRAY(k_array, R_unrelesed_array, R_relesed_array, numM):
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

def get_FER_ARRAY(fer, members_T, numM, numC):
    FER_array = []
    for i in range(numM):
        invT = np.linalg.inv(members_T[i])
        sub_FER_array = []
        for j in range(numC):
            sub_FER_array.append(np.matmul(invT, fer[i, j]))
        FER_array.append(sub_FER_array)
    return FER_array

def get_fer_ARRAY(model, k12, k22, fer1, fer2, numM, numC):
    ferCondensed_ARRAY = []
    for i in range(numM):
        k_12, k_22, fer_1, fer_2 = k12[i], k22[i], fer1[i], fer2[i]
        ferCondensed = np.subtract(fer_1, np.matmul(np.matmul(k_12, np.linalg.inv(k_22)), fer_2))[:,0] # TODO figure out why there is extra brakets
        j = 0
        product = []
        for a in range(numC):
            sub_ferCondensed = ferCondensed[a]
            for DOF in model.members_Releases[i]:
                if DOF:
                    sub_ferCondensed = np.insert(sub_ferCondensed, j, 0, axis=0)
                j += 1
            product.append(sub_ferCondensed)
        ferCondensed_ARRAY.append(product)
    return np.array(ferCondensed_ARRAY)

def memberPart_fer_ARRAY(fer, R_unrelesed, R_relesed):
    fer1_array = []
    fer2_array = []
    for i in range(len(fer)):
        fer1_a = []
        fer2_a = []
        for j in range(len(fer[i])):
            fer1_a.append(fer[i][j][R_unrelesed])
            fer2_a.append(fer[i][j][R_relesed])
        fer1_array.append(fer1_a)
        fer2_array.append(fer2_a)
    return np.array(fer1_array), np.array(fer2_array)

def get_fer_unc_ARRAY(model, pointLoads, distLoads):
    fer_unc_ARRAY = np.zeros((len(model.members), 12, 1))
    fer_unc_ARRAY = np.add(fer_unc_ARRAY, get_FixedEndReactions_Pointload_ARRAY(model, pointLoads))
    fer_unc_ARRAY = np.add(fer_unc_ARRAY, get_FixedEndReactions_Distload_ARRAY(model, distLoads))
    return fer_unc_ARRAY

def get_FixedEndReactions_Pointload_ARRAY(model, pointLoads):
    # pointLoads [memberINDEX, caseINDEX, location, data], data: [x, Px, Py, Pz, Mx, My, Mz]
    Reactions = []
    for i, loadCasses in enumerate(pointLoads):
        L = model.members_L[i]
        reactionsMember = []
        for j, loads in enumerate(loadCasses):
            reactions = np.zeros((12, 1))
            if loads is not None:
                for k, load in enumerate(loads):
                    if load[1] != 0:
                        reactions = np.add(reactions, ferCalc.pointLoadX(load[1], load[0], L))
                    if load[2] != 0:
                        reactions = np.add(reactions, ferCalc.pointLoadY(load[2], load[0], L))
                    if load[3] != 0:
                        reactions = np.add(reactions, ferCalc.pointLoadZ(load[2], load[0], L))
                    if load[4] != 0:
                        reactions = np.add(reactions, ferCalc.momentX(load[3], load[0], L))
                    if load[5] != 0:
                        reactions = np.add(reactions, ferCalc.momentY(load[4], load[0], L))
                    if load[6] != 0:
                        reactions = np.add(reactions, ferCalc.momentZ(load[5], load[0], L))
            reactionsMember.append(reactions)
        Reactions.append(reactionsMember)
    return Reactions

def get_FixedEndReactions_Distload_ARRAY(model, distLoads):
    # distLoads [memberINDEX, caseINDEX, location, data], data: [x1, x2, wx1, wx2, wy1, wy2, wz1, wz2]
    Reactions = []
    for i, loadCasses in enumerate(distLoads):
        L = model.members_L[i]
        reactionsMember = []
        for j, loads in enumerate(loadCasses):
            reactions = np.zeros((12, 1))
            for k, load in enumerate(loads):
                if not (load[2] == 0 and load[3] == 0):
                    reactions = np.add(reactions, ferCalc.distributedLoadX([load[2],load[3]], [load[0], load[1]], L))
                if not (load[4] == 0 and load[5] == 0):
                    reactions = np.add(reactions, ferCalc.distributedLoadY([load[4],load[5]], [load[0], load[1]], L))
                if not (load[6] == 0 and load[7] == 0):
                    reactions = np.add(reactions, ferCalc.distributedLoadZ([load[6],load[7]], [load[0], load[1]], L))
            reactionsMember.append(reactions)
        Reactions.append(reactionsMember)
    return Reactions

def assembleLoads(model: Frame3D_T):
    pointLoads = assemblePointLoads(model)
    distLoads = assembleDistLoads(model)
    return pointLoads, distLoads

def assemblePointLoads(model: Frame3D_T):
    # pointLoads [memberINDEX, caseINDEX, location, data], data: [x, Px, Py, Pz, Mx, My, Mz]
    newPointLoads = []
    for i, loadCasses in enumerate(model.members_PointLoads):
        newLoads  = [None] * len(model.casses)
        for loads in loadCasses:
            newLoad  = []
            for j in range(len(loads[1][0])):
                newLoad.append([loads[1][0][j]] + loads[1][1][j])
            newLoads[loads[0]] = newLoad
        newPointLoads.append(newLoads)
    return newPointLoads

def assembleDistLoads(model: Frame3D_T):
    # distLoads [memberINDEX, caseINDEX, location, data], data: [x1, x2, wx1, wx2, wy1, wy2, wz1, wz2]
    newDistLoads = []
    for i, loadCasses in enumerate(model.members_DistLoads):
        newLoads = [None] * len(model.casses)
        for loads in loadCasses:
            newLoad = []
            for j in range(len(loads[1])):
                newLoad.append(loads[1][j][0] + loads[1][j][1])
            newLoads[loads[0]] = newLoad
        for j in range(len(newLoads)):
            if newLoads[j] is None:
                # noinspection PyTypeChecker
                newLoads[j] = []
        newDistLoads.append(newLoads)
    return newDistLoads

def partitionedGlobalNodalForceVector(model, numN):
    p_array = []
    for loadCaseI in model.casses:
        p = np.zeros((numN * 6, 1))
        for i in range(numN):
            if np.size(model.nodes_loads[loadCaseI]) != 0:
                local = np.array(model.nodes_loads[loadCaseI][1])
            else:
                local = np.zeros(6, dtype=float)
            dofs = buildDOFVector([i])
            p[dofs, 0] += local
        p_array.append(p)
    return partFER(p_array, model.D_unknown, model.D_known, len(model.casses))

def k_member_make_global(k_local_array, T_array):
    k_member_global = []
    for i in range(len(T_array)):
        k_member_global.append(np.matmul(np.matmul(np.linalg.inv(T_array[i]), k_local_array[i]), T_array[i]))
    return k_member_global

def get_K_Global(model, k_global, numN, numM):
    K = np.zeros((numN * 6, numN * 6))
    for i in range(numM):
        K[np.ix_(model.members_DOF[i], model.members_DOF[i])] += np.asarray(k_global[i], dtype=float)
    return K

def partition_K_gloabl(K_global, R_unrelesed, R_relesed):
    K11 = (K_global[R_unrelesed, :][:, R_unrelesed])
    K12 = (K_global[R_unrelesed, :][:, R_relesed])
    K21 = (K_global[R_relesed, :][:, R_unrelesed])
    K22 = (K_global[R_relesed, :][:, R_relesed])
    return K11, K12, K21, K22

def get_D(K11, K12, P1_array, FER1_array, Index_Unsupported, Index_Supported, numN, numC):
    D1_array = []
    D2 = np.zeros((len(Index_Supported),1))
    for i in range(numC):
        P1 = np.array(P1_array[i])
        FER1 = np.array(FER1_array[i])
        if K11.shape == (0, 0):
            D1_array.append([])
        else:
            try:
                D1_array.append(np.linalg.solve(K11, np.subtract(np.subtract(P1, FER1), np.matmul(K12, D2))))
            except:
                raise Exception('The stiffness matrix is singular, which implies rigid body motion. The structure is unstable. Aborting analysis.')
    D_array = assemble_D_array(D1_array, D2, Index_Supported, Index_Unsupported, numN, numC)
    DX_array, DY_array, DZ_array, RX_array, RY_array, RZ_array = get_direcction_deflections(D_array, numN, numC)
    return D_array, DX_array, DY_array, DZ_array, RX_array, RY_array, RZ_array

def assemble_D_array(D1_array, D2, Index_Supported, Index_Unsupported, numN, numC):
    D_array = []
    for a in range(numC):
        D1 = np.array(D1_array[a])
        D = np.zeros((numN * 6, 1))
        for j in range(numN):
            for i in range(6):
                if j * 6 + i in Index_Supported:
                    D[(j * 6 + i, 0)] = D2[Index_Supported.index(j * 6 + i), 0]
                else:
                    D[(j * 6 + i, 0)] = D1[Index_Unsupported.index(j * 6 + i), 0]
        D_array.append(D)
    return D_array


def get_direcction_deflections(D_array, numN, numC):
    DX_array = []
    DY_array = []
    DZ_array = []
    RX_array = []
    RY_array = []
    RZ_array = []
    for i in range(numC):
        DX_case = []
        DY_case = []
        DZ_case = []
        RX_case = []
        RY_case = []
        RZ_case = []
        for j in range(numN):
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
    return DX_array, DY_array, DZ_array, RX_array, RY_array, RZ_array
import math
import numpy as np

def segment_Member(members, members_L, members_CrossSectionProps, materials, pointLoads, distLoads, f_array, fer_array, d_array, numM, numC):
    seg = []
    seg_InternalLoads = []
    seg_DistLoads = []
    seg_thata = []
    seg_delta = []
    for mINDEX in range(numM):
        seg_sub = []
        seg_sub_InternalLoads = []
        seg_sub_DistLoads = []
        seg_sub_thata = []
        seg_sub_delta = []
        for cINDEX in range(numC):
            # Create a list of discontinuity locations
            disconts = [0, float(members_L[mINDEX])]
            for load in pointLoads[mINDEX][cINDEX]:
                disconts.append(load[0])
            for load in distLoads[mINDEX][cINDEX]:
                disconts.append(load[0])  # Distributed load start locations
                disconts.append(load[1])  # Distributed load end locations
            disconts = sorted(set(disconts)) # Sort the list and eliminate duplicate values
            E = materials[members[mINDEX,2],0]
            A = members_CrossSectionProps[mINDEX,0]
            Iy = members_CrossSectionProps[mINDEX,1]
            Iz = members_CrossSectionProps[mINDEX,2]
            J = members_CrossSectionProps[mINDEX,3]
            L = members_L[mINDEX]
            seg_x1 = []
            seg_x2 = []
            seg_EIz = []
            seg_EIy = []
            seg_EA = []
            for index in range(len(disconts)-1):
                seg_x1.append(disconts[index])
                seg_x2.append(disconts[index+1])
                seg_EIz.append(E * Iz)
                seg_EIy.append(E * Iy)
                seg_EA.append(E * A)
            d = d_array[cINDEX, mINDEX]
            f = f_array[cINDEX, mINDEX]
            fer = fer_array[mINDEX, cINDEX]
            m1z = f[5,0]
            m2z = f[11,0]
            m1y = -f[4,0]
            m2y = -f[10,0]
            fem1z = fer[5,0]
            fem2z = fer[11,0]
            fem1y = -fer[4,0]
            fem2y = -fer[10,0]
            delta1y = d[1]
            delta2y = d[7]
            delta1z = d[2]
            delta2z = d[8]
            seg_deltaY1, seg_deltaZ1 = [delta1z], [delta1y]
            seg_thetaY1 = [-1 / 3 * ((m1y - fem1y) * L / (E * Iy) - (m2y - fem2y) * L / (2 * E * Iy) + 3 * (delta2z - delta1z) / L)]
            seg_thetaZ1 = [ 1 / 3 * ((m1z - fem1z) * L / (E * Iz) - (m2z - fem2z) * L / (2 * E * Iz) + 3 * (delta2y - delta1y) / L)]
            seg_deltaX1 = [d[0]]
            seg_L = []
            seg_distLoad_X, seg_distLoad_Y, seg_distLoad_Z = [], [], []
            seg_InternalLoad_X, seg_InternalLoad_Y, seg_InternalLoad_Z = [], [], [] # P1, X_T1 # V1, M1 # V1, M1
            for i in range(len(seg_x1)):
                x = seg_x1[i]
                seg_L.append(seg_x2[i] - seg_x1[i])

                seg_InternalLoad_X.append([f[0, 0], f[3, 0]])  # P1, X_T1
                seg_InternalLoad_Y.append([f[2, 0], f[4, 0] + f[2, 0] * x])  # V1, M1
                seg_InternalLoad_Z.append([f[1, 0], f[5, 0] - f[1, 0] * x])  # V1, M1

                seg_distLoad_Z.append([0,0]) # w1, w2
                seg_distLoad_Y.append([0,0]) # w1, w2
                seg_distLoad_X.append([0,0]) # w1, p2

                if i > 0:
                    seg_thetaZ1.append(seg_thetaZ1[i-1] - (-seg_InternalLoad_Z[i-1][0] * seg_L[i-1] **2 / 2 - seg_distLoad_Z[i-1][0] * seg_L[i-1] ** 3/6 + seg_L[i-1] * seg_InternalLoad_Z[i-1][1] + seg_L[i-1] **4*(seg_distLoad_Z[i-1][0] - seg_distLoad_Z[i-1][1])/(24 * seg_L[i-1])) / seg_EIz[i-1])
                    seg_deltaZ1.append((seg_deltaZ1[i-1] + seg_thetaZ1[i-1]* seg_L[i-1] + seg_InternalLoad_Z[i-1][0]*seg_L[i-1]**3/(6*seg_EIz[i-1]) + seg_distLoad_Z[i-1][0]*x**4/(24*seg_EIz[i-1]) +seg_L[i-1]**2*(-seg_InternalLoad_Z[i-1][1])/(2*seg_EIz[i-1]) + x**5*(-seg_distLoad_Z[i-1][0] + seg_distLoad_Z[i-1][1])/(120*seg_EIz[i-1]*seg_L[i-1])))
                    seg_thetaY1.append(seg_thetaY1[i-1] + (-seg_InternalLoad_Y[i-1][0]*seg_L[i-1]**2/2 - seg_distLoad_Y[i-1][0]*seg_L[i-1]**3/6 + seg_L[i-1]*(-seg_InternalLoad_Y[i-1][1]) + seg_L[i-1]**4*(seg_distLoad_Y[i-1][0] - seg_distLoad_Y[i-1][1])/(24*seg_L[i-1]))/seg_EIy[i-1])
                    seg_deltaY1.append((seg_deltaY1[i-1] - seg_thetaY1[i-1]*seg_L[i-1] + seg_InternalLoad_Y[i-1][0]*seg_L[i-1]**3/(6*seg_EIy[i-1]) + seg_distLoad_Y[i-1][0]*seg_L[i-1]**4/(24*seg_EIy[i-1]) -seg_L[i-1]**2*(-seg_distLoad_Y[i-1][1])/(2*seg_EIy[i-1]) - seg_L[i-1]**5*(seg_distLoad_Y[i-1][0] - seg_distLoad_Y[i-1][1])/(120*seg_EIy[i-1]*seg_L[i-1])))
                    seg_deltaX1.append(seg_deltaX1[i-1] - 1/seg_EA[i-1]*(seg_InternalLoad_X[i-1][0]*seg_L[i-1] + seg_distLoad_X[i-1][0]*seg_L[i-1]**2/2 + (seg_distLoad_X[i-1][1] - seg_distLoad_X[i-1][0])*seg_L[i-1]**3/(6*seg_L[i-1])))
                for pointLoad in pointLoads[mINDEX][cINDEX]:
                    if round(pointLoad[0], 10) <= round(x, 10):
                        seg_InternalLoad_X[i][0] += pointLoad[1]
                        seg_InternalLoad_Y[i][0] += pointLoad[2]
                        seg_InternalLoad_Y[i][1] -= pointLoad[2] * (x - pointLoad[0])
                        seg_InternalLoad_Z[i][0] += pointLoad[3]
                        seg_InternalLoad_Z[i][1] -= pointLoad[3] * (x - pointLoad[0])
                        seg_InternalLoad_X[i][1] += pointLoad[4]
                        seg_InternalLoad_Y[i][1] += pointLoad[5]
                        seg_InternalLoad_Z[i][1] += pointLoad[6]
                for distLoad in distLoads[mINDEX][cINDEX]:
                    if round(distLoad[0], 10) <= round(x, 10):
                        if round(distLoad[1], 10) > round(x, 10):
                            #X Direction
                            seg_distLoad_X[i][0] += (distLoad[3] - distLoad[2]) / (distLoad[1] - distLoad[0]) * (x - distLoad[0]) + distLoad[2]
                            seg_distLoad_X[i][1] += (distLoad[3] - distLoad[2]) / (distLoad[1] - distLoad[0]) * (seg_x2[i] - distLoad[0]) + distLoad[2]
                            seg_InternalLoad_X[i][0] += (distLoad[2] + (distLoad[2] + (distLoad[3] - distLoad[2]) / (distLoad[1] - distLoad[0]) * (x - distLoad[0]))) / 2 * (x - distLoad[0])
                            # Y Direction
                            seg_distLoad_Z[i][0] += (distLoad[5] - distLoad[4]) / (distLoad[1] - distLoad[0]) * (x - distLoad[0]) + distLoad[4]
                            seg_distLoad_Z[i][1] += (distLoad[5] - distLoad[4]) / (distLoad[1] - distLoad[0]) * (seg_x2[i] - distLoad[0]) + distLoad[4]
                            w2_int = (distLoad[4] + (distLoad[5] - distLoad[4]) / (distLoad[1] - distLoad[0]) * (x - distLoad[0]))
                            seg_InternalLoad_Z[i][0] += (distLoad[4] + w2_int) / 2 * (x - distLoad[0])
                            seg_InternalLoad_Z[i][1] -= (distLoad[0] - x) * (2 * distLoad[4] * distLoad[0] - 3 * distLoad[4] * x + distLoad[4] * x + w2_int * distLoad[0] - 3 * w2_int * x + 2 * w2_int * x) / 6
                            # Z Direction
                            seg_distLoad_Y[i][0] += (distLoad[7] - distLoad[6]) / (distLoad[1] - distLoad[0]) * (x - distLoad[0]) + distLoad[6]
                            seg_distLoad_Y[i][1] += (distLoad[7] - distLoad[6]) / (distLoad[1] - distLoad[0]) * (seg_x2[i] - distLoad[0]) + distLoad[6]
                            w2_int = distLoad[3][6] + (distLoad[7] - distLoad[6]) / (distLoad[1] - distLoad[0]) * (x - distLoad[0])
                            seg_InternalLoad_Y[i][0] += (distLoad[6] + w2_int) / 2 * (x - distLoad[0])
                            seg_InternalLoad_Y[i][1] += (distLoad[0] - x) * (2 * distLoad[6] * distLoad[0] - 3 * distLoad[6] * x + distLoad[6] * x + w2_int * distLoad[0] - 3 * w2_int * x + 2 * w2_int * x) / 6
                        else:
                            #X Direction
                            seg_InternalLoad_X[i][0] += (distLoad[2] + distLoad[3]) / 2 * (distLoad[1] - distLoad[0])
                            # Y direction
                            seg_InternalLoad_Z[i][0] += (distLoad[4] + distLoad[5]) / 2 * (distLoad[1] - distLoad[0])
                            seg_InternalLoad_Z[i][1] -= (distLoad[0] - distLoad[1]) * (2 * distLoad[4] * distLoad[0] - 3 * distLoad[4] * x + distLoad[4] * distLoad[1] + distLoad[5] * distLoad[0] - 3 * distLoad[5] * x + 2 * distLoad[5] * distLoad[1]) / 6
                            # z Direction
                            seg_InternalLoad_Y[i][0] += (distLoad[6] + distLoad[7]) / 2 * (distLoad[1] - distLoad[0])
                            seg_InternalLoad_Y[i][1] += (distLoad[0] - distLoad[1]) * (2 * distLoad[6] * distLoad[0] - 3 * distLoad[6] * x + distLoad[6] * distLoad[1] + distLoad[7] * distLoad[0] - 3 * distLoad[7] * x + 2 * distLoad[7] * distLoad[1]) / 6
            seg_sub.append([seg_x1, seg_x2, seg_EIz, seg_EIy, seg_EA])
            seg_sub_InternalLoads.append([seg_InternalLoad_X, seg_InternalLoad_Y, seg_InternalLoad_Z]) # [P1, T1], [Vy1, My1], [Vz1, Mz1]
            seg_sub_DistLoads.append([seg_distLoad_X, seg_distLoad_Y, seg_distLoad_Z]) # [wx1, wx2], [wy1, wy2], [wz1, wz2]
            seg_sub_thata.append([seg_thetaY1, seg_thetaZ1])
            seg_sub_delta.append([seg_deltaX1, seg_deltaY1, seg_deltaZ1])
        seg.append(seg_sub)
        seg_InternalLoads.append(seg_sub_InternalLoads)
        seg_DistLoads.append(seg_sub_DistLoads)
        seg_thata.append(seg_sub_thata)
        seg_delta.append(seg_sub_delta)
    return seg, seg_InternalLoads, seg_DistLoads, seg_thata, seg_delta

def extremFinder(model, mINDEX, seg, seg_InternalLoads, seg_DistLoads, absFunction, absFunctionDirection, POIFunction, Direc, sign = None, comboINDEXS = None):
    if comboINDEXS is None: comboINDEXS = model.casses
    global_, governing_combo = None, None
    seg, seg_InternalLoads, seg_DistLoads = seg[mINDEX], seg_InternalLoads[mINDEX], seg_DistLoads[mINDEX]
    for comboINDEX in comboINDEXS:
        abs_ = []
        L = np.array(seg[comboINDEX][1]) - np.array(seg[comboINDEX][0])
        for segINDEX in range(len(seg)):
            abs_.append(absFunction(POIFunction(seg_DistLoads[comboINDEX][Direc][segINDEX][0], seg_DistLoads[comboINDEX][Direc][segINDEX][1], seg_InternalLoads[comboINDEX][Direc][segINDEX][0], seg_InternalLoads[comboINDEX][Direc][segINDEX][1], L[comboINDEX], sign)))
        abs_ = absFunction(abs_)
        if absFunctionDirection == "max":
            if global_ is None or abs_ > global_: global_, governing_combo = abs_, comboINDEX
        elif absFunctionDirection == "min":
            if global_ is None or abs_ < global_: global_, governing_combo = abs_, comboINDEX
    return global_, governing_combo


""" --------------- SHEAR --------------- """


def shear(model, x, mINDEX, comboINDEX, seg, seg_InternalLoads, seg_DistLoads):
    seg = seg[mINDEX][comboINDEX]
    seg_InternalLoads = seg_InternalLoads[mINDEX][comboINDEX]
    seg_DistLoads = seg_DistLoads[mINDEX][comboINDEX]
    for i in range(len(seg[0])):
        if round(seg[0][i], 10) <= round(x, 10) < round(seg[1][i], 10):
            VY = seg_InternalLoads[1][i][0] + seg_DistLoads[1][i][0] * (x - seg[0][i]) + (x - seg[0][i]) ** 2 * (-seg_DistLoads[1][i][0] + seg_DistLoads[1][i][1]) / (2 * (seg[1][i]-seg[0][i]))
            VZ = seg_InternalLoads[2][i][0] + seg_DistLoads[2][i][0] * (x - seg[0][i]) + (x - seg[0][i]) ** 2 * (-seg_DistLoads[2][i][0] + seg_DistLoads[2][i][1]) / (2 * (seg[1][i] - seg[0][i]))
            return VY, VZ
    if math.isclose(x, model.members_L[mINDEX]):
        lastIndex = len(seg[0]) - 1
        VY = seg_InternalLoads[1][lastIndex][0] + seg_DistLoads[1][lastIndex][0] * (x - seg[0][lastIndex]) + (x - seg[0][lastIndex]) ** 2 * (-seg_DistLoads[1][lastIndex][0] + seg_DistLoads[1][lastIndex][1]) / (2 * (seg[1][lastIndex] - seg[0][lastIndex]))
        VZ = seg_InternalLoads[2][lastIndex][0] + seg_DistLoads[2][lastIndex][0] * (x - seg[0][lastIndex]) + (x - seg[0][lastIndex]) ** 2 * (-seg_DistLoads[2][lastIndex][0] + seg_DistLoads[2][lastIndex][1]) / (2 * (seg[1][lastIndex] - seg[0][lastIndex]))
        return VY, VZ
    return 0, 0

def shearLocOfInterist(w1, w2, L):
    if w1 - w2 == 0: x1 = 0
    else: x1 = w1 * L / (w1 - w2)
    if round(x1, 10) < 0 or round(x1, 10) > round(L, 10): x1 = 0
    return x1

def seg_V_POI(w1, w2, V1, M_1, L, sign):
    x = shearLocOfInterist(w1, w2, L)
    shear1 = V1 + w1 * x + x ** 2 * (-w1 + w2) / (2 * L)
    shear2 = V1
    shear3 = V1 + w1 * L + L ** 2 * (-w1 + w2) / (2 * L)
    return [shear1, shear2, shear3]

def max_shear(model, mINDEX, seg, seg_InternalLoads, seg_DistLoads, comboINDEXS = None):
    Y, yC = extremFinder(model, mINDEX, seg, seg_InternalLoads, seg_DistLoads, max, "max",seg_V_POI, 1, comboINDEXS=comboINDEXS)
    Z, zC = extremFinder(model, mINDEX, seg, seg_InternalLoads, seg_DistLoads, max, "max", seg_V_POI, 2, comboINDEXS=comboINDEXS)
    return Y, Z, yC, zC

def min_shear(model, mINDEX, seg, seg_InternalLoads, seg_DistLoads, comboINDEXS = None):
    Y, yC = extremFinder(model, mINDEX, seg, seg_InternalLoads, seg_DistLoads, min, "min", seg_V_POI, 1, comboINDEXS=comboINDEXS)
    Z, zC = extremFinder(model, mINDEX, seg, seg_InternalLoads, seg_DistLoads, min, "min", seg_V_POI, 2, comboINDEXS=comboINDEXS)
    return Y, Z, yC, zC


""" --------------- MOMENT --------------- """


def moment_calc(M1, V1, w1, w2, L, x, sign):
    return sign*M1 - V1*x - w1*x**2/2 - x**3*(-w1 + w2)/(6*L)

def moment(model, x, mINDEX, comboINDEX, seg, seg_InternalLoads, seg_DistLoads):
    seg = seg[mINDEX][comboINDEX]
    seg_InternalLoads = seg_InternalLoads[mINDEX][comboINDEX]
    seg_DistLoads = seg_DistLoads[mINDEX][comboINDEX]
    for i in range(len(seg[0])):
        if round(seg[0][i], 10) <= round(x, 10) < round(seg[1][i], 10):
            MZ = moment_calc(seg_InternalLoads[2][i][1], seg_InternalLoads[2][i][0], seg_DistLoads[2][i][0], seg_DistLoads[2][i][1], (seg[1][i]-seg[0][i]), (x - seg[0][i]), 1)
            MY = moment_calc(seg_InternalLoads[1][i][1], seg_InternalLoads[1][i][0], seg_DistLoads[1][i][0], seg_DistLoads[1][i][1], (seg[1][i]-seg[0][i]), (x - seg[0][i]), -1)
            return MZ, MY
    if math.isclose(x, model.members_L[mINDEX]):
        i = len(seg[0]) - 1
        MZ = moment_calc(seg_InternalLoads[2][i][1], seg_InternalLoads[2][i][0], seg_DistLoads[2][i][0],seg_DistLoads[2][i][1], (seg[1][i] - seg[0][i]), (x - seg[0][i]), 1)
        MY = moment_calc(seg_InternalLoads[1][i][1], seg_InternalLoads[1][i][0], seg_DistLoads[1][i][0],seg_DistLoads[1][i][1], (seg[1][i] - seg[0][i]), (x - seg[0][i]), -1)
        return MY, MZ
    return 0, 0

def momentLocOfInterist(w1, w2, V1, L):
    a = -(w2 - w1) / (2 * L)
    b = -w1
    c = -V1
    if a == 0:
        if b != 0: x1 = -c / b
        else: x1 = 0
        x2 = 0
    elif b ** 2 - 4 * a * c < 0:
        x1, x2 = 0, 0
    else:
        x1 = (-b + (b ** 2 - 4 * a * c) ** 0.5) / (2 * a)
        x2 = (-b - (b ** 2 - 4 * a * c) ** 0.5) / (2 * a)
    if round(x1, 10) < 0 or round(x1, 10) > round(L, 10): x1 = 0
    if round(x2, 10) < 0 or round(x2, 10) > round(L, 10): x2 = 0
    return x1, x2

def seg_M_POI(w1, w2, V1, M_1, L, sign):
    x1, x2 = momentLocOfInterist(w1, w2, V1, L)
    M1 = moment_calc(M_1, V1, w1, w2, L, x1, sign)
    M2 = moment_calc(M_1, V1, w1, w2, L, x2, sign)
    M3 = moment_calc(M_1, V1, w1, w2, L, 0, sign)
    M4 = moment_calc(M_1, V1, w1, w2, L, L, sign)
    return [M1, M2, M3, M4]

def max_moment(model, mINDEX, seg, seg_InternalLoads, seg_DistLoads, comboINDEXS = None):
    Y, yC = extremFinder(model, mINDEX, seg, seg_InternalLoads, seg_DistLoads, max, "max",seg_M_POI, 1, sign=-1, comboINDEXS = comboINDEXS)
    Z, zC = extremFinder(model, mINDEX, seg, seg_InternalLoads, seg_DistLoads, max, "max",seg_M_POI, 2, sign=1, comboINDEXS=comboINDEXS)
    return Y, Z, yC, zC

def min_moment(model, mINDEX, seg, seg_InternalLoads, seg_DistLoads, comboINDEXS = None):
    Y, yC = extremFinder(model, mINDEX, seg, seg_InternalLoads, seg_DistLoads, min, "min",seg_M_POI, 1, sign=-1, comboINDEXS=comboINDEXS)
    Z, zC = extremFinder(model, mINDEX, seg, seg_InternalLoads, seg_DistLoads, min, "min",seg_M_POI, 2, sign=1, comboINDEXS=comboINDEXS)
    return Y, Z, yC, zC


""" --------------- TORQUE --------------- """


def torque(model, x, mINDEX, comboINDEX, seg, seg_InternalLoads):
    seg = seg[mINDEX][comboINDEX]
    seg_InternalLoads = seg_InternalLoads[mINDEX][comboINDEX]
    for i in range(len(seg[0])):
        if round(seg[0][i], 10) <= round(x, 10) < round(seg[1][i], 10):
            return seg_InternalLoads[0][i][1]
    if math.isclose(x, model.members_L[mINDEX]):
        i = len(seg[0]) - 1
        return seg_InternalLoads[0][i][1]
    return 0

def seg_T_POI(w1, w2, P1, M1, L, sign):
    return [M1]

def max_tourque(model, mINDEX, seg, seg_InternalLoads, seg_DistLoads, comboINDEXS = None):
    return extremFinder(model, mINDEX, seg, seg_InternalLoads, seg_DistLoads, max, "max",seg_T_POI, 0, comboINDEXS=comboINDEXS)

def min_tourque(model, mINDEX, seg, seg_InternalLoads, seg_DistLoads, comboINDEXS = None):
    return extremFinder(model, mINDEX, seg, seg_InternalLoads, seg_DistLoads, min, "min",seg_T_POI, 0, comboINDEXS=comboINDEXS)


""" --------------- AXIAL --------------- """


def axial_calc(p1, p2, P1, L, x):
    return P1 + (p2 - p1)/(2*L)*x**2 + p1*x

def axial(model, x, mINDEX, comboINDEX, seg, seg_InternalLoads, seg_DistLoads):
    seg = seg[mINDEX][comboINDEX]
    seg_InternalLoads = seg_InternalLoads[mINDEX][comboINDEX]
    seg_DistLoads = seg_DistLoads[mINDEX][comboINDEX]
    for i in range(len(seg[0])):
        if round(seg[0][i], 10) <= round(x, 10) < round(seg[1][i], 10):
            return axial_calc(seg_DistLoads[0][1][0], seg_DistLoads[0][1][1], seg_InternalLoads[0][i][0], (seg[1][i] - seg[0][i]), (x-seg[0][i]))
    if math.isclose(x, model.members_L[mINDEX]):
        i = len(seg[0]) - 1
        return axial_calc(seg_DistLoads[0][1][0], seg_DistLoads[0][1][1], seg_InternalLoads[0][i][0], (seg[1][i] - seg[0][i]), (x-seg[0][i]))
    return 0

def axialLocOfInterist(p1, p2, L):
    if p1 - p2 != 0:
        x1 = L * p1 / (p1 - p2)
    else:
        x1 = 0
    if round(x1, 10) < 0 or round(x1, 10) > round(L, 10):
        x1 = 0
    return x1

def seg_P_POI(w1, w2, P1, M1, L, sign):
    x1 = axialLocOfInterist(w1, w2, L)
    P_1 = axial_calc(w1, w2, P1, L, x1)
    P_2 = axial_calc(w1, w2, P1, L, 0)
    P_3 = axial_calc(w1, w2, P1, L, L)
    return [P_1, P_2, P_3]

def max_axial(model, mINDEX, seg, seg_InternalLoads, seg_DistLoads, comboINDEXS = None):
    return extremFinder(model, mINDEX, seg, seg_InternalLoads, seg_DistLoads, max, "max",seg_P_POI, 0, comboINDEXS)

def min_axial(model, mINDEX, seg, seg_InternalLoads, seg_DistLoads, comboINDEXS = None):
    return extremFinder(model, mINDEX, seg, seg_InternalLoads, seg_DistLoads, min, "min",seg_P_POI, 0, comboINDEXS)


""" --------------- DEFLECTION --------------- """


def axial_deflection_calc(delta_x1, EA, P1, w1, w2, L, x):
    return delta_x1 - 1/EA*(P1 * x + w1 * x ** 2 / 2 + (w2 - w1) * x ** 3 / (6 * L))

def deflection_calc(delta1, theta1, V1, EI, w1, w2, M1, L, sign, x):
    return delta1 + sign* theta1*x + V1*x**3/(6 * EI) + w1*x**4/(24 * EI) + sign*x**2*(-M1)/(2 * EI) + sign*x**5*(sign* (-w1) + sign*w2)/(120 * EI * L)

def deflection(model, x, mINDEX, comboINDEX, seg, seg_InternalLoads, seg_DistLoads, seg_thata, seg_delta):
    seg = seg[mINDEX][comboINDEX]
    seg_InternalLoads = seg_InternalLoads[mINDEX][comboINDEX]
    seg_DistLoads = seg_DistLoads[mINDEX][comboINDEX]
    seg_thata = seg_thata[mINDEX][comboINDEX]
    seg_delta = seg_delta[mINDEX][comboINDEX]
    for i in range(len(seg[0])):
        if round(seg[0][i], 10) <= round(x, 10) < round(seg[1][i], 10):
            DX = axial_deflection_calc(seg_delta[0][i], seg[4][i], seg_InternalLoads[0][i][0], seg_DistLoads[0][i][0], seg_DistLoads[0][i][1], (seg[1][i] - seg[0][i]), (x-seg[0][i]))
            DY = deflection_calc(seg_delta[1][i], seg_thata[0][i], seg_InternalLoads[1][i][0], seg[3][i], seg_DistLoads[1][i][0], seg_DistLoads[1][i][1], seg_InternalLoads[1][i][1], (seg[1][i] - seg[0][i]), -1, (x-seg[0][i]))
            DZ = deflection_calc(seg_delta[2][i], seg_thata[1][i], seg_InternalLoads[2][i][0], seg[2][i], seg_DistLoads[2][i][0], seg_DistLoads[2][i][1], seg_InternalLoads[2][i][1], (seg[1][i] - seg[0][i]), 1, (x-seg[0][i]))
            return [DX, DY, DZ]
    if math.isclose(x, model.members_L[mINDEX]):
        i = len(seg[0]) - 1
        DX = axial_deflection_calc(seg_delta[0][i], seg[4][i], seg_InternalLoads[0][i][0], seg_DistLoads[0][i][0], seg_DistLoads[0][i][1], (seg[1][i] - seg[0][i]), (x-seg[0][i]))
        DY = deflection_calc(seg_delta[1][i], seg_thata[0][i], seg_InternalLoads[1][i][0], seg[3][i], seg_DistLoads[1][i][0], seg_DistLoads[1][i][1], seg_InternalLoads[1][i][1], (seg[1][i] - seg[0][i]), -1, (x-seg[0][i]))
        DZ = deflection_calc(seg_delta[2][i], seg_thata[1][i], seg_InternalLoads[2][i][0], seg[2][i], seg_DistLoads[2][i][0], seg_DistLoads[2][i][1], seg_InternalLoads[2][i][1], (seg[1][i] - seg[0][i]), 1, (x-seg[0][i]))
        return [DX, DY, DZ]
    return [0,0,0]

def difflectionExtremFinder(model, mINDEX, seg, seg_InternalLoads, seg_DistLoads, seg_thata, seg_delta, absFunctionDirection, comboINDEXS, members_L):
    if comboINDEXS is None: comboINDEXS = model.casses
    global_, governing_combo = [None, None, None], [None, None, None]
    L = members_L[mINDEX]
    for comboINDEX in comboINDEXS:
        dmax = deflection(model, 0, mINDEX, comboINDEX, seg, seg_InternalLoads, seg_DistLoads, seg_thata, seg_delta)
        for i in range(100):
            d = deflection(model, L * i / 99, mINDEX, comboINDEX, seg, seg_InternalLoads, seg_DistLoads, seg_thata, seg_delta)
            for j in range(3):
                if absFunctionDirection == "max":
                    if d[j] > dmax[j]:
                        dmax[j] = d[j]
                elif absFunctionDirection == "min":
                    if d[j] < dmax[j]:
                        dmax[j] = d[j]
        for j in range(3):
            if global_[j] is None:
                global_[j] = dmax[j]
                governing_combo[j] = comboINDEX
            else:
                if absFunctionDirection == "max":
                    # noinspection PyTypeChecker
                    if dmax[j] > global_[j]:
                        global_[j] = dmax[j]
                        governing_combo[j] = comboINDEX
                elif absFunctionDirection == "min":
                    # noinspection PyTypeChecker
                    if dmax[j] < global_[j]:
                        global_[j] = dmax[j]
                        governing_combo[j] = comboINDEX
    return global_, governing_combo

def max_difflection(model, mINDEX, seg, seg_InternalLoads, seg_DistLoads, seg_thata, seg_delta, comboINDEXS = None):
    return difflectionExtremFinder(model, mINDEX, seg, seg_InternalLoads, seg_DistLoads, seg_thata, seg_delta,"max", comboINDEXS, model.members_L)

def min_difflection(model, mINDEX, seg, seg_InternalLoads, seg_DistLoads, seg_thata, seg_delta, comboINDEXS = None):
    return difflectionExtremFinder(model, mINDEX, seg, seg_InternalLoads, seg_DistLoads, seg_thata, seg_delta,"min", comboINDEXS, model.members_L)
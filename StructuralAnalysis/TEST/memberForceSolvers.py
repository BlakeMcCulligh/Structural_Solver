import math
import numpy as np

from StructuralAnalysis.TEST.__main__ import Frame3D_T

def segment_Member(model: Frame3D_T, pointLoads, distLoads, f_array, fer_array, d_array, numM, numC):
    seg = []
    seg_InternalLoads = []
    seg_DistLoads = []
    for mINDEX in range(numM):
        seg_sub = []
        seg_sub_InternalLoads = []
        seg_sub_DistLoads = []
        for cINDEX in range(numC):
            # Create a list of discontinuity locations
            disconts = [0, model.members_L[mINDEX]]
            for load in pointLoads[mINDEX][cINDEX]:
                disconts.append(load[0])
            for load in distLoads[mINDEX][cINDEX]:
                disconts.append(load[0])  # Distributed load start locations
                disconts.append(load[1])  # Distributed load end locations
            disconts = sorted(set(disconts)) # Sort the list and eliminate duplicate values
            E = model.materials[model.members[mINDEX,2],0]
            A = model.members_CrossSectionProps[mINDEX,0]
            Iy = model.members_CrossSectionProps[mINDEX,1]
            Iz = model.members_CrossSectionProps[mINDEX,2]
            J = model.members_CrossSectionProps[mINDEX,3]
            L = model.members_L[mINDEX]
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
            for i in range(len(seg_x1)-1):
                x = seg_x1[i]
                seg_L.append(seg_x1[i+1] - seg_x1[i])
                seg_distLoad_Z.append([0,0]) # w1, w2
                seg_distLoad_Y.append([0,0]) # w1, w2
                seg_distLoad_X.append([0,0]) # p1, p2
                if i > 0:
                    seg_thetaZ1.append(seg_thetaZ1[i-1] - (-seg_InternalLoad_Z[i-1][0] * seg_L[i-1] **2 / 2 - seg_distLoad_Z[i-1][0] * seg_L[i-1] ** 3/6 + seg_L[i-1] * seg_InternalLoad_Z[i-1][1] + seg_L[i-1] **4*(seg_distLoad_Z[i-1][0] - seg_distLoad_Z[i-1][1])/(24 * seg_L[i-1])) / seg_EIz[i-1])
                    seg_deltaZ1.append((seg_deltaZ1[i-1] + seg_thetaZ1[i-1]* seg_L[i-1] + seg_InternalLoad_Z[i-1][0]*seg_L[i-1]**3/(6*seg_EIz[i-1]) + seg_distLoad_Z[i-1][0]*x**4/(24*seg_EIz[i-1]) +seg_L[i-1]**2*(-seg_InternalLoad_Z[i-1][1])/(2*seg_EIz[i-1]) + x**5*(-seg_distLoad_Z[i-1][0] + seg_distLoad_Z[i-1][1])/(120*seg_EIz[i-1]*seg_L[i-1])))
                    seg_thetaY1.append(seg_thetaY1[i-1] + (-seg_InternalLoad_Y[i-1][0]*seg_L[i-1]**2/2 - seg_distLoad_Y[i-1][0]*seg_L[i-1]**3/6 + seg_L[i-1]*(-seg_InternalLoad_Y[i-1][1]) + seg_L[i-1]**4*(seg_distLoad_Y[i-1][0] - seg_distLoad_Y[i-1][1])/(24*seg_L[i-1]))/seg_EIy[i-1])
                    seg_deltaY1.append((seg_deltaY1[i-1] - seg_thetaY1[i-1]*seg_L[i-1] + seg_InternalLoad_Y[i-1][0]*seg_L[i-1]**3/(6*seg_EIy[i-1]) + seg_distLoad_Y[i-1][0]*seg_L[i-1]**4/(24*seg_EIy[i-1]) -seg_L[i-1]**2*(-seg_distLoad_Y[i-1][1])/(2*seg_EIy[i-1]) - seg_L[i-1]**5*(seg_distLoad_Y[i-1][0] - seg_distLoad_Y[i-1][1])/(120*seg_EIy[i-1]*seg_L[i-1])))
                    seg_deltaX1.append(seg_deltaX1[i-1] - 1/seg_EA[i-1]*(seg_InternalLoad_X[i-1][0]*seg_L[i-1] + seg_distLoad_X[i-1][0]*seg_L[i-1]**2/2 + (seg_distLoad_X[i-1][1] - seg_distLoad_X[i-1][0])*seg_L[i-1]**3/(6*seg_L[i-1])))
                seg_InternalLoad_X.append([f[0, 0],f[3, 0]]) # P1, X_T1
                seg_InternalLoad_Y.append([f[2, 0],f[4, 0] + f[2, 0] * x]) # V1, M1
                seg_InternalLoad_Z.append([f[1, 0],f[5, 0] - f[1, 0] * x]) # V1, M1
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
        seg.append(seg_sub)
        seg_InternalLoads.append(seg_sub_InternalLoads)
        seg_DistLoads.append(seg_sub_DistLoads)
    return seg, seg_InternalLoads, seg_DistLoads

def shear(model, x, mINDEX, comboINDEX, seg, seg_InternalLoads, seg_DistLoads):
    seg = seg[mINDEX][comboINDEX]
    seg_InternalLoads = seg_InternalLoads[mINDEX][comboINDEX]
    seg_DistLoads = seg_DistLoads[mINDEX][comboINDEX]
    for i in range(len(seg)):
        if round(seg[i][0], 10) <= round(x, 10) < round(seg[i][1], 10):
            VY = seg_InternalLoads[1][i][0] + seg_DistLoads[1][i][0] * (x - seg[i][0]) + (x - seg[i][0]) ** 2 * (-seg_DistLoads[1][i][0] + seg_DistLoads[1][i][1]) / (2 * (seg[i][1]-seg[i][0]))
            VZ = seg_InternalLoads[2][i][0] + seg_DistLoads[2][i][0] * (x - seg[i][0]) + (x - seg[i][0]) ** 2 * (-seg_DistLoads[2][i][0] + seg_DistLoads[2][i][1]) / (2 * (seg[i][1] - seg[i][0]))
            return VY, VZ
    if math.isclose(x, model.members_L[mINDEX]):
        lastIndex = len(seg) - 1
        VY = seg_InternalLoads[1][lastIndex][0] + seg_DistLoads[1][lastIndex][0] * (x - seg[lastIndex][0]) + (x - seg[lastIndex][0]) ** 2 * (-seg_DistLoads[1][lastIndex][0] + seg_DistLoads[1][lastIndex][1]) / (2 * (seg[lastIndex][1]-seg[lastIndex][0]))
        VZ = seg_InternalLoads[2][lastIndex][0] + seg_DistLoads[2][lastIndex][0] * (x - seg[lastIndex][0]) + (x - seg[lastIndex][0]) ** 2 * (-seg_DistLoads[2][lastIndex][0] + seg_DistLoads[2][lastIndex][1]) / (2 * (seg[lastIndex][1] - seg[lastIndex][0]))
        return VY, VZ
    return 0, 0

def seg_V_POI(i, d_L, i_L, L, d_i):
    w1, w2, V1 = d_L[d_i][i][0], d_L[d_i][i][1], i_L[d_i][i][0]
    x = shearLocOfInterist(w1, w2, L)
    shear1 = V1 + w1 * x + x ** 2 * (-w1 + w2) / (2 * L)
    shear2 = V1
    shear3 = V1 + w1 * L + L ** 2 * (-w1 + w2) / (2 * L)
    return [shear1, shear2, shear3]

def shearLocOfInterist(w1, w2, L):
    if w1 - w2 == 0: x1 = 0
    else: x1 = w1 * L / (w1 - w2)
    if round(x1, 10) < 0 or round(x1, 10) > round(L, 10): x1 = 0
    return x1

def max_shear(model, mINDEX, seg, seg_InternalLoads, seg_DistLoads, comboINDEXS = None):
    if comboINDEXS is None: comboINDEXS = model.casses
    Vmax_global_Y, Vmax_global_Z, governing_combo_Y, governing_combo_Z = None, None, None, None
    seg, seg_InternalLoads, seg_DistLoads = seg[mINDEX], seg_InternalLoads[mINDEX], seg_DistLoads[mINDEX]
    for comboINDEX in comboINDEXS:
        Vmax_Y, Vmax_Z = [], []
        L = np.array(seg[comboINDEX][1]) - np.array(seg[comboINDEX][0])
        for segINDEX in range(len(seg)):
            Vmax_Y.append(max(seg_V_POI(segINDEX, seg_DistLoads[comboINDEX], seg_InternalLoads[comboINDEX], L[segINDEX], 1)))
            Vmax_Z.append(max(seg_V_POI(segINDEX, seg_DistLoads[comboINDEX], seg_InternalLoads[comboINDEX], L[segINDEX], 2)))
        Vmax_Y, Vmax_Z = min(Vmax_Y), min(Vmax_Z)
        if Vmax_global_Y is None or Vmax_Y > Vmax_global_Y: Vmax_global_Y, governing_combo_Y = Vmax_Y, comboINDEX
        if Vmax_global_Z is None or Vmax_Z > Vmax_global_Z: Vmax_global_Z, governing_combo_Z = Vmax_Z, comboINDEX
    return Vmax_global_Y, governing_combo_Y, Vmax_global_Z, governing_combo_Z

def min_shear(model, mINDEX, seg, seg_InternalLoads, seg_DistLoads, comboINDEXS=None):
    if comboINDEXS is None: comboINDEXS = model.casses
    Vmin_global_Y, Vmin_global_Z, governing_combo_Y, governing_combo_Z = None, None, None, None
    seg, seg_InternalLoads, seg_DistLoads = seg[mINDEX], seg_InternalLoads[mINDEX], seg_DistLoads[mINDEX]
    for comboINDEX in comboINDEXS:
        Vmin_Y, Vmin_Z = [], []
        L = np.array(seg[comboINDEX][1]) - np.array(seg[comboINDEX][0])
        for segINDEX in range(len(seg)):
            Vmin_Y.append(min(seg_V_POI(segINDEX, seg_DistLoads[comboINDEX], seg_InternalLoads[comboINDEX], L[segINDEX], 1)))
            Vmin_Z.append(min(seg_V_POI(segINDEX, seg_DistLoads[comboINDEX], seg_InternalLoads[comboINDEX], L[segINDEX], 2)))
        Vmin_Y, Vmin_Z = min(Vmin_Y), min(Vmin_Z)
        if Vmin_global_Y is None or Vmin_Y < Vmin_global_Y: Vmin_global_Y, governing_combo_Y = Vmin_Y, comboINDEX
        if Vmin_global_Z is None or Vmin_Z < Vmin_global_Z: Vmin_global_Z, governing_combo_Z = Vmin_Z, comboINDEX
    return Vmin_global_Y, governing_combo_Y, Vmin_global_Z, governing_combo_Z
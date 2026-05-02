import numpy as np

def pointLoadX_batch(F, x, L):
    F = np.asarray(F)
    x = np.asarray(x)
    b = L - x
    n = len(F)
    R = np.zeros((n, 12))
    R[:, 0] = -F * b / L
    R[:, 6] = -F * x / L
    return R

def pointLoadY_batch(F, x, L):
    F = np.asarray(F)
    x = np.asarray(x)
    b = L - x
    n = len(F)
    R = np.zeros((n, 12))
    R[:, 1]  = -F * b**2 * (L + 2*x) / L**3
    R[:, 5]  = -F * x * b**2 / L**2
    R[:, 7]  = -F * x**2 * (L + 2*b) / L**3
    R[:,11]  =  F * x**2 * b / L**2
    return R

def pointLoadZ_batch(F, x, L):
    F = np.asarray(F)
    x = np.asarray(x)
    b = L - x
    n = len(F)
    R = np.zeros((n, 12))
    R[:, 2]  = -F * b**2 * (L + 2*x) / L**3
    R[:, 4]  =  F * x * b**2 / L**2
    R[:, 8]  = -F * x**2 * (L + 2*b) / L**3
    R[:,10]  = -F * x**2 * b / L**2
    return R

def momentX_batch(M, x, L):
    M = np.asarray(M)
    x = np.asarray(x)
    b = L - x
    n = len(M)
    R = np.zeros((n, 12))
    R[:, 3] = -M * b / L
    R[:, 9] = -M * x / L
    return R

def momentY_batch(M, x, L):
    M = np.asarray(M)
    x = np.asarray(x)
    b = L - x
    n = len(M)
    R = np.zeros((n, 12))
    R[:, 2]  = -6 * M * x * b / L**3
    R[:, 4]  =  M * b * (2*x - b) / L**2
    R[:, 8]  =  6 * M * x * b / L**3
    R[:,10]  =  M * x * (2*b - x) / L**2
    return R

def momentZ_batch(M, x, L):
    M = np.asarray(M)
    x = np.asarray(x)
    b = L - x
    n = len(M)
    R = np.zeros((n, 12))
    R[:, 1]  =  6 * M * x * b / L**3
    R[:, 5]  =  M * b * (2*x - b) / L**2
    R[:, 7]  = -6 * M * x * b / L**3
    R[:,11]  =  M * x * (2*b - x) / L**2
    return R

def distributedLoadX(w: list, loc: list, L: float):
    reactions = np.zeros((12, 1))
    reactions[0, 0] = 1 / (6 * L) * (loc[0] - loc[1]) * (
                3 * L * w[0] + 3 * L * w[1] - 2 * w[0] * loc[0] - w[0] * loc[1] - w[1] * loc[0] - 2 * w[1] * loc[1])
    reactions[6, 0] = 1 / (6 * L) * (loc[0] - loc[1]) * (
                2 * w[0] * loc[0] + w[0] * loc[1] + w[1] * loc[0] + 2 * w[1] * loc[1])
    return reactions

def distributedLoadY(w: list, loc: list, L: float):
    reactions = np.zeros((12, 1))
    reactions[1, 0] = (loc[0] - loc[1]) * (
            10 * L ** 3 * w[0] + 10 * L ** 3 * w[1] - 15 * L * w[0] * loc[0] ** 2 - 10 * L * w[0] * loc[0] * loc[
        1] - 5 * L * w[0] * loc[1] ** 2 - 5 * L * w[1] * loc[0] ** 2 - 10 * L * w[1] * loc[0] * loc[1] - 15 * L * w[1] *
            loc[1] ** 2 + 8 * w[0] * loc[0] ** 3 + 6 * w[0] * loc[0] ** 2 * loc[1] + 4 * w[0] * loc[0] * loc[
                1] ** 2 + 2 * w[0] * loc[1] ** 3 + 2 * w[1] * loc[0] ** 3 + 4 * w[1] * loc[0] ** 2 * loc[1] + 6 * w[1] *
            loc[0] * loc[1] ** 2 + 8 * w[1] * loc[1] ** 3) / (20 * L ** 3)
    reactions[5, 0] = (loc[0] - loc[1]) * (
            20 * L ** 2 * w[0] * loc[0] + 10 * L ** 2 * w[0] * loc[1] + 10 * L ** 2 * w[1] * loc[0] + 20 * L ** 2 * w[
        1] * loc[1] - 30 * L * w[0] * loc[0] ** 2 - 20 * L * w[0] * loc[0] * loc[1] - 10 * L * w[0] * loc[
                1] ** 2 - 10 * L * w[1] * loc[0] ** 2 - 20 * L * w[1] * loc[0] * loc[1] - 30 * L * w[1] * loc[
                1] ** 2 + 12 * w[0] * loc[0] ** 3 + 9 * w[0] * loc[0] ** 2 * loc[1] + 6 * w[0] * loc[0] * loc[
                1] ** 2 + 3 * w[0] * loc[1] ** 3 + 3 * w[1] * loc[0] ** 3 + 6 * w[1] * loc[0] ** 2 * loc[1] + 9 * w[1] *
            loc[0] * loc[1] ** 2 + 12 * w[1] * loc[1] ** 3) / (60 * L ** 2)
    reactions[7, 0] = -(loc[0] - loc[1]) * (
            -15 * L * w[0] * loc[0] ** 2 - 10 * L * w[0] * loc[0] * loc[1] - 5 * L * w[0] * loc[1] ** 2 - 5 * L * w[1] *
            loc[0] ** 2 - 10 * L * w[1] * loc[0] * loc[1] - 15 * L * w[1] * loc[1] ** 2 + 8 * w[0] * loc[0] ** 3 + 6 *
            w[0] * loc[0] ** 2 * loc[1] + 4 * w[0] * loc[0] * loc[1] ** 2 + 2 * w[0] * loc[1] ** 3 + 2 * w[1] * loc[
                0] ** 3 + 4 * w[1] * loc[0] ** 2 * loc[1] + 6 * w[1] * loc[0] * loc[1] ** 2 + 8 * w[1] * loc[
                1] ** 3) / (20 * L ** 3)
    reactions[11, 0] = (loc[0] - loc[1]) * (
            -15 * L * w[0] * loc[0] ** 2 - 10 * L * w[0] * loc[0] * loc[1] - 5 * L * w[0] * loc[1] ** 2 - 5 * L * w[1] *
            loc[0] ** 2 - 10 * L * w[1] * loc[0] * loc[1] - 15 * L * w[1] * loc[1] ** 2 + 12 * w[0] * loc[0] ** 3 + 9 *
            w[0] * loc[0] ** 2 * loc[1] + 6 * w[0] * loc[0] * loc[1] ** 2 + 3 * w[0] * loc[1] ** 3 + 3 * w[1] * loc[
                0] ** 3 + 6 * w[1] * loc[0] ** 2 * loc[1] + 9 * w[1] * loc[0] * loc[1] ** 2 + 12 * w[1] * loc[
                1] ** 3) / (60 * L ** 2)
    return reactions

def distributedLoadZ(w: list, loc: list, L: float):
    reactions = np.zeros((12, 1))
    reactions[2, 0] = (loc[0] - loc[1]) * (
            10 * L ** 3 * w[0] + 10 * L ** 3 * w[1] - 15 * L * w[0] * loc[0] ** 2 - 10 * L * w[0] * loc[0] * loc[
        1] - 5 * L * w[0] * loc[1] ** 2 - 5 * L * w[1] * loc[0] ** 2 - 10 * L * w[1] * loc[0] * loc[1] - 15 * L * w[1] *
            loc[1] ** 2 + 8 * w[0] * loc[0] ** 3 + 6 * w[0] * loc[0] ** 2 * loc[1] + 4 * w[0] * loc[0] * loc[
                1] ** 2 + 2 * w[0] * loc[1] ** 3 + 2 * w[1] * loc[0] ** 3 + 4 * w[1] * loc[0] ** 2 * loc[1] + 6 * w[1] *
            loc[0] * loc[1] ** 2 + 8 * w[1] * loc[1] ** 3) / (20 * L ** 3)
    reactions[4, 0] = -(loc[0] - loc[1]) * (
            20 * L ** 2 * w[0] * loc[0] + 10 * L ** 2 * w[0] * loc[1] + 10 * L ** 2 * w[1] * loc[0] + 20 * L ** 2 * w[
        1] * loc[1] - 30 * L * w[0] * loc[0] ** 2 - 20 * L * w[0] * loc[0] * loc[1] - 10 * L * w[0] * loc[
                1] ** 2 - 10 * L * w[1] * loc[0] ** 2 - 20 * L * w[1] * loc[0] * loc[1] - 30 * L * w[1] * loc[
                1] ** 2 + 12 * w[0] * loc[0] ** 3 + 9 * w[0] * loc[0] ** 2 * loc[1] + 6 * w[0] * loc[0] * loc[
                1] ** 2 + 3 * w[0] * loc[1] ** 3 + 3 * w[1] * loc[0] ** 3 + 6 * w[1] * loc[0] ** 2 * loc[1] + 9 * w[1] *
            loc[0] * loc[1] ** 2 + 12 * w[1] * loc[1] ** 3) / (60 * L ** 2)
    reactions[8, 0] = (-(loc[0] - loc[1]) * (
            -15 * L * w[0] * loc[0] ** 2 - 10 * L * w[0] * loc[0] * loc[1] - 5 * L * w[0] * loc[1] ** 2 - 5 * L * w[1] *
            loc[0] ** 2 - 10 * L * w[1] * loc[0] * loc[1] - 15 * L * w[1] * loc[1] ** 2 + 8 * w[0] * loc[0] ** 3 + 6 *
            w[0] * loc[0] ** 2 * loc[1] + 4 * w[0] * loc[0] * loc[1] ** 2 + 2 * w[0] * loc[1] ** 3 + 2 * w[1] *
            loc[0] ** 3 + 4 * w[1] * loc[0] ** 2 * loc[1] + 6 * w[1] * loc[0] * loc[1] ** 2 + 8 * w[1] * loc[1] ** 3)
                       / (20 * L ** 3))
    reactions[10, 0] = -(loc[0] - loc[1]) * (
            -15 * L * w[0] * loc[0] ** 2 - 10 * L * w[0] * loc[0] * loc[1] - 5 * L * w[0] * loc[1] ** 2 - 5 * L * w[1] *
            loc[0] ** 2 - 10 * L * w[1] * loc[0] * loc[1] - 15 * L * w[1] * loc[1] ** 2 + 12 * w[0] * loc[0] ** 3 + 9 *
            w[0] * loc[0] ** 2 * loc[1] + 6 * w[0] * loc[0] * loc[1] ** 2 + 3 * w[0] * loc[1] ** 3 + 3 * w[1] * loc[
                0] ** 3 + 6 * w[1] * loc[0] ** 2 * loc[1] + 9 * w[1] * loc[0] * loc[1] ** 2 + 12 * w[1] * loc[
                1] ** 3) / (60 * L ** 2)
    return reactions
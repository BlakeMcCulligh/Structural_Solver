import numpy as np

def pointLoadX(magnatude: float, location: float, L: float):
    reactions = np.zeros((12, 1))
    reactions[0, 0] = -magnatude * (L - location) / L
    reactions[6, 0] = -magnatude * location / L
    return reactions

def pointLoadY(magnatude: float, location: float, L: float):
    b = L - location
    reactions = np.zeros((12, 1))
    reactions[1, 0] = -magnatude* b ** 2 * (L + 2 * location) / L ** 3
    reactions[5, 0] = -magnatude* location * b ** 2 / L ** 2
    reactions[7, 0] = -magnatude* location ** 2 * (L + 2 * b) / L ** 3
    reactions[11, 0] = magnatude* location ** 2 * b / L ** 2
    return reactions

def pointLoadZ(magnatude: float, location: float, L: float):
    b = L - location
    reactions = np.zeros((12, 1))
    reactions[2, 0] = -magnatude* b ** 2 * (L + 2 * location) / L ** 3
    reactions[4, 0] = magnatude* location * b ** 2 / L ** 2
    reactions[8, 0] = -magnatude* location ** 2 * (L + 2 * b) / L ** 3
    reactions[10, 0] = -magnatude* location ** 2 * b / L ** 2
    return reactions

def momentX(magnatude: float, location: float, L: float):
    reactions = np.zeros((12, 1))
    reactions[3, 0] = -magnatude * (L - location) / L
    reactions[9, 0] = -magnatude * location / L
    return reactions

def momentY(magnatude: float, location: float, L: float):
    b = L - location
    reactions = np.zeros((12, 1))
    reactions[2, 0] = -6 * magnatude * location * b / L ** 3
    reactions[4, 0] = magnatude * b * (2 * location - b) / L ** 2
    reactions[8, 0] = 6 * magnatude * location * b / L ** 3
    reactions[10, 0] = magnatude * location * (2 * b - location) / L ** 2
    return reactions

def momentZ(magnatude: float, location: float, L: float):
    b = L - location
    reactions = np.zeros((12, 1))
    reactions[1, 0] = 6 * magnatude * location * b / L ** 3
    reactions[5, 0] = magnatude * b * (2 * location - b) / L ** 2
    reactions[7, 0] = -6 * magnatude * location * b / L ** 3
    reactions[11, 0] = magnatude * location * (2 * b - location) / L ** 2
    return reactions

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
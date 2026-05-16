"""
Holds functions used to calculate the end reactions of a member if it were fixed.
"""

import numpy as np

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

def point_load_x(F: np.ndarray, x: np.ndarray, L: float) -> np.ndarray:
    """
    Calculates the reactions of a member from x direction point loads.

    :param F: ndarray. Force applied to the member.
    :param x: ndarray. Location of the force applied to the member.
    :param L: float. Length of the member.
    :return: ndarray, Reactions of the member.
    """

    b = L - x
    n = len(F)

    R = np.zeros((n, 12))

    R[:, 0] = -F * b / L
    R[:, 6] = -F * x / L

    return R

def point_load_y(F: np.ndarray, x: np.ndarray, L: float) -> np.ndarray:
    """
    Calculates the reactions of a member from y direction point loads.

    :param F: ndarray. Force applied to the member.
    :param x: ndarray. Location of the force applied to the member.
    :param L: float. Length of the member.
    :return: ndarray, Reactions of the member.
    """

    b = L - x
    n = len(F)

    R = np.zeros((n, 12))

    R[:, 1]  = -F * b**2 * (L + 2*x) / L**3
    R[:, 5]  = -F * x * b**2 / L**2
    R[:, 7]  = -F * x**2 * (L + 2*b) / L**3
    R[:,11]  =  F * x**2 * b / L**2

    return R

def point_load_z(F: np.ndarray, x: np.ndarray, L: float) -> np.ndarray:
    """
    Calculates the reactions of a member from z direction point loads.

    :param F: ndarray. Force applied to the member.
    :param x: ndarray. Location of the force applied to the member.
    :param L: float. Length of the member.
    :return: ndarray, Reactions of the member.
    """

    b = L - x
    n = len(F)

    R = np.zeros((n, 12))

    R[:, 2]  = -F * b**2 * (L + 2*x) / L**3
    R[:, 4]  =  F * x * b**2 / L**2
    R[:, 8]  = -F * x**2 * (L + 2*b) / L**3
    R[:,10]  = -F * x**2 * b / L**2

    return R

def moment_x(M: np.ndarray, x: np.ndarray, L: float) -> np.ndarray:
    """
    Calculates the reactions of a member from x direction moments.

    :param M: ndarray. Moment applied to the member.
    :param x: ndarray. Location of the force applied to the member.
    :param L: float. Length of the member.
    :return: ndarray, Reactions of the member.
    """

    b = L - x
    n = len(M)

    R = np.zeros((n, 12))

    R[:, 3] = -M * b / L
    R[:, 9] = -M * x / L

    return R

def moment_y(M: np.ndarray, x: np.ndarray, L: float) -> np.ndarray:
    """
    Calculates the reactions of a member from y direction moments.

    :param M: ndarray. Moment applied to the member.
    :param x: ndarray. Location of the force applied to the member.
    :param L: float. Length of the member.
    :return: ndarray, Reactions of the member.
    """

    b = L - x
    n = len(M)

    R = np.zeros((n, 12))

    R[:, 2]  = -6 * M * x * b / L**3
    R[:, 4]  =  M * b * (2*x - b) / L**2
    R[:, 8]  =  6 * M * x * b / L**3
    R[:,10]  =  M * x * (2*b - x) / L**2

    return R

def moment_z(M: np.ndarray, x: np.ndarray, L: float) -> np.ndarray:
    """
    Calculates the reactions of a member from z direction moments.

    :param M: ndarray. Moment applied to the member.
    :param x: ndarray. Location of the force applied to the member.
    :param L: float. Length of the member.
    :return: ndarray, Reactions of the member.
    """

    b = L - x
    n = len(M)

    R = np.zeros((n, 12))

    R[:, 1]  =  6 * M * x * b / L**3
    R[:, 5]  =  M * b * (2*x - b) / L**2
    R[:, 7]  = -6 * M * x * b / L**3
    R[:,11]  =  M * x * (2*b - x) / L**2

    return R

def distributed_load_x(w: list, loc: list, L: float) -> np.ndarray:
    """
    Calculates the reactions of a member from an x direction distributed load.

    :param w: list. Distributed load represented by the magnitude at both ends of it.
    :param loc: list. Start and end points of the distributed load.
    :param L: float. Length of the member.
    :return: ndarray, Reactions of the member.
    """

    reactions = np.zeros((12, 1))

    reactions[0, 0] = 1 / (6 * L) * (loc[0] - loc[1]) * (
                3 * L * w[0] + 3 * L * w[1] - 2 * w[0] * loc[0] - w[0] * loc[1] - w[1] * loc[0] - 2 * w[1] * loc[1])

    reactions[6, 0] = 1 / (6 * L) * (loc[0] - loc[1]) * (
                2 * w[0] * loc[0] + w[0] * loc[1] + w[1] * loc[0] + 2 * w[1] * loc[1])

    return reactions

def distributed_load_y(w: list, loc: list, L: float) -> np.ndarray:
    """
    Calculates the reactions of a member from a y direction distributed load.

    :param w: list. Distributed load represented by the magnitude at both ends of it.
    :param loc: list. Start and end points of the distributed load.
    :param L: float. Length of the member.
    :return: ndarray, Reactions of the member.
    """

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

def distributed_load_z(w: list, loc: list, L: float) -> np.ndarray:
    """
    Calculates the reactions of a member from an x direction distributed load.

    :param w: list. Distributed load represented by the magnitude at both ends of it.
    :param loc: list. Start and end points of the distributed load.
    :param L: float. Length of the member.
    :return: ndarray, Reactions of the member.
    """

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
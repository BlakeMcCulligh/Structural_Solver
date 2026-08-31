"""
Holds functions to calculate properties of square-HSS cross-sections.
"""

__author__ = "Blake McCulligh"
__copyright__ = "Copyright 2026 Blake McCulligh"
__credits__ = ["Blake McCulligh"]

__license__ = "MIT"
__version__ = "0.1.0b1"
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = "Beta"

def get_A(d: float, t: float) -> float:
    """
    Calculates the area of a square-HSS cross-section.

    :param d: float. Depth
    :param t: float. Thickness
    :return: float. Cross-Section Area
    """

    return 4*t*(d - t)

def get_I(d: float, t: float) -> float:
    """
    Calculates the moment of inertia of a square-HSS cross-section.

    :param d: float. Depth
    :param t: float. Thickness
    :return: float. Moment of Inertia
    """

    return 1/12*(d**4-(d-2*t)**4)

def get_J(d: float, t: float) -> float:
    """
    Calculates the torsional constant of a square-HSS cross-section.

    :param d: float. Diameter
    :param t: float. Thickness
    :return: float. Torsional Constant
    """

    return t*d**3
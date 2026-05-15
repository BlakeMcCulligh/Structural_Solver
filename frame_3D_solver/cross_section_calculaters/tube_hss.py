"""
Holds functions to calculate properties of tube-HSS cross-sections.
"""

from math import pi

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

def get_A(d: float, t: float) -> float:
    """
    Calculates the area of a tube-HSS cross-section.

    :param d: float. Diameter
    :param t: float. Thickness
    :return: float. Cross-Section Area
    """

    return pi*(d**2-(d-2*t)**2)/4

def get_I(d: float, t: float) -> float:
    """
    Calculates the moment of inertia of a tube-HSS cross-section.

    :param d: float. Diameter
    :param t: float. Thickness
    :return: float. Moment of Inertia
    """

    return pi*(d**4-(d-2*t)**4)/64

def get_J(d: float, t: float) -> float:
    """
    Calculates the torsional constant of a tube-HSS cross-section.

    :param d: float. Diameter
    :param t: float. Thickness
    :return: float. Torsional Constant
    """

    return pi*(d**4-(d-2*t)**4)/32
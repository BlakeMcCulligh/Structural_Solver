"""
Holds functions to calculate properties of rectangular-HSS cross-sections.
"""

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

def get_A(d: float, b: float, t: float) -> float:
    """
    Calculates the area of a rectangular-HSS cross-section.

    :param d: float. Depth
    :param b: float. Width
    :param t: float. Thickness
    :return: float. Cross-Section Area
    """

    return b*d-(b-2*t)*(d-2*t)

def get_Ix(d: float, b: float, t: float) -> float:
    """
    Calculates the moment of inertia around the strong axis of a rectangular-HSS cross-section.

    :param d: float. Depth
    :param b: float. Width
    :param t: float. Thickness
    :return: float. Moment of Inertia X
    """

    return 1/12*(b*d**3-(b-2*t)*(d-2*t)**3)

def get_Iy(d: float, b: float, t: float) -> float:
    """
    Calculates the moment of inertia around the weak axis of a rectangular-HSS cross-section.

    :param d: float. Depth
    :param b: float. Width
    :param t: float. Thickness
    :return: float. Moment of Inertia Y
    """

    return 1/12*(d*b**3-(d-2*t)*(b-2*t)**3)

def get_J(d: float, b: float, t: float) -> float:
    """
    Calculates the torsional constant of a rectangular-HSS cross-section.

    :param d: float. Depth
    :param b: float. Width
    :param t: float. Thickness
    :return: float. Torsional Constant
    """

    return 2*t*d**2*b**2/(d+b)
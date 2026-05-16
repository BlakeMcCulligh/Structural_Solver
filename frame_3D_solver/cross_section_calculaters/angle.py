"""
Holds functions to calculate properties of angle cross-sections.
"""

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

def _get_x(d: float, b: float, t: float) -> float:
    """
    Calculates the horizontal distance form the outside of the corner to the center of the angle.

    :param d: float. Depth
    :param b: float. Width
    :param t: float. Thickness
    :return: float. Horizontal distance form the outside of the corner to the center of the angle.
    """

    return (b**2+(d-t)*t)/(2*(b+d-t))

def _get_y(d: float, b: float, t: float) -> float:
    """
    Calculates the vertical distance form the outside of the corner to the center of the angle.

    :param d: float. Depth
    :param b: float. Width
    :param t: float. Thickness
    :return: float. Vertical distance form the outside of the corner to the center of the angle.
    """

    return (d**2+(b-t)*t)/(2*(b+d-t))

def get_A(d: float, b: float, t: float) -> float:
    """
    Calculates the area of an angle cross-section.

    :param d: float. Depth
    :param b: float. Width
    :param t: float. Thickness
    :return: float. Cross-Section Area
    """

    return t*(b+d-t)

def get_Ix(d: float, b: float, t: float) -> float:
    """
    Calculates the moment of inertia around the strong axis of an angle cross-section.

    :param d: float. Depth
    :param b: float. Width
    :param t: float. Thickness
    :return: float. Moment of Inertia X
    """

    y = _get_y(d, b, t)
    return (t*(d-y)**3+b*y**3-(b-t)*(y-t)**3)/3

def get_Iy(d: float, b: float, t: float) -> float:
    """
    Calculates the moment of inertia around the weak axis of an angle cross-section.

    :param d: float. Depth
    :param b: float. Width
    :param t: float. Thickness
    :return: float. Moment of Inertia Y
    """

    x = _get_x(d, b, t)
    return (t*(b-x)**3+d*x**3-(d-t)*(x-t)**3)/3

def get_J(d: float, b: float, t: float) -> float:
    """
    Calculates the torsional constant of an angle cross-section.
    Assumes t is small compared to d and b

    :param d: float. Depth
    :param b: float. Width
    :param t: float. Thickness
    :return: float. Torsional Constant
    """

    return t**3*(b+d)/3
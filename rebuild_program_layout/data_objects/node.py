"""
Handels everything to do with a node.
"""

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

class Node:
    """
    Object storing a node.
    """

    def __init__(self, x: float = 0, y: float = 0, z: float = 0, cords: list[float] | None = None):
        """
        Constructor: sets the location of the node.

        :param x: x coordinate of the node.
        :type x: float
        :param y: y coordinate of the node.
        :type y: float
        :param z: z coordinate of the node.
        :type z: float
        :param cords:  Alternative to input a list of the coordanits.
        :type cords: list[float]
        """

        if cords is None:
            self.x: float = x
            self.y: float = y
            self.z: float = z
        else:
            self.x: float = cords[0]
            self.y: float = cords[1]
            self.z: float = cords[2]

        self.support: list[bool] = [False, False, False, False, False, False]

        self.load: list[list[float | int]] = []

    def get_cords(self) -> list[float]:
        """
        Return the coordinates of the node.

        :return: list of the coordinates of the node.
        :rtype: list[float]
        """

        return [self.x, self.y, self.z]

    def set_support(self, support: list[float]) -> None:
        """
        Set the support of the node. False for unsupported, True for supported.

        :param support: list of support degress of freedom. [DX, DY, DZ, RX, RY, RZ]
        :type support: list[float]
        """

        self.support = support

    def add_load(self, case: int, load: list[float | int]) -> None:
        """
        Adds a load to the node.

        :param case: Load case the load is to be in.
        :type case: int
        :param load: list of magnatudes of the load in each direction. [DX, DY, DZ, RX, RY, RZ]
        :type load: list[float]
        """

        self.load.append(load + [case])
"""
Handels everything to do with a material.
"""

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

class Material:
    """
    Object storing a material.
    """

    def __init__(self, E: float, G: float, nu: float, rho: float, fy: float | None = None):

        self.E: float = E
        self.G: float = G
        self.nu: float = nu
        self.rho: float = rho
        self.fy: float | None = fy
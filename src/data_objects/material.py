"""
Handels everything to do with a material.
"""

__author__ = "Blake McCulligh"
__copyright__ = "Copyright 2026 Blake McCulligh"
__credits__ = ["Blake McCulligh"]

__license__ = "MIT"
__version__ = "0.1.0b1"
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = "Beta"

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
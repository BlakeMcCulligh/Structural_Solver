"""
Handels everything to do with a member.
"""

from rebuild_program_layout.data_objects.node import Node
from rebuild_program_layout.frame_3D import Frame3D

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

class Member:
    """
    Object storing a member.
    """

    def __init__(self, node_1_index: int, node_2_index: int, material_index: int, set_cross_section_props: bool,
                 cross_section_props: list[float], frame: Frame3D):
        """
        Constructor. sets the data for a member.

        :param node_1_index: Index of the first node of the member.
        :type node_1_index: int
        :param node_2_index: Index of the second node of the member.
        :type node_2_index: int
        :param material_index: Index of the material of the member.
        :type material_index: int
        :param set_cross_section_props: If the cross-section properties should be set for an optimization or
                                        if they should be optimized.
        :type set_cross_section_props: bool
        :param cross_section_props: Cross-section properties of the member. [A, Iy, Iz, J]
        :type cross_section_props: list[float]
        :param frame: 3D frame the member is in.
        :type frame: Frame3D
        """

        self.frame = frame

        self.node_1_index = node_1_index
        self.node_2_index = node_2_index
        self.material_index = material_index
        self.set_cross_section_props = set_cross_section_props
        self.A = cross_section_props[0]
        self.Iy = cross_section_props[1]
        self.Iz = cross_section_props[2]
        self.J = cross_section_props[3]

        self.node_1: Node = self.frame.nodes[self.node_1_index]
        self.node_2: Node = self.frame.nodes[self.node_2_index]

        self.materials = self.frame.materials[self.material_index]




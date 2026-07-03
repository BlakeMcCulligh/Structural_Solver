"""
Handels all operations to do with the 3D frame.
"""
import pandas as pd

from rebuild_program_layout.data_objects.member import Member
from rebuild_program_layout.data_objects.material import Material
from rebuild_program_layout.data_objects.node import Node
from rebuild_program_layout.optimization import start_optimization

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

class Frame3D:
    def __init__(self, controller):

        self._controller = controller

        self.nodes: list[Node] = []
        self.materials: list[Material] = []
        self.members: list[Member] = []

    def linear_analysis(self):
        pass

    def GlobalOptimization(self, group_assignments, group_types, lower_bounds, upper_bounds,
                                                      cost_function, weight_run, reaction_run, internal_forces_run):
        # tempory testing
        bounds = [(-5, 5), (-5, 5)]
        start_optimization(self._controller.executor, self._controller.root, bounds)

        pass

    def save(self, file_path):
        pass

    def export_results(self, file_path):
        pass

    def open_frame(self, file_path):
        pass

    def open_results(self, file_path):
        pass

    def import_nodes(self, file_path: str) -> None:
        """
        Imports coordinates for nodes from the specified Excel file.
        Headers must be "X, Y, Z" and there must only be one sheet.

        :param file_path: file path to the Excel file to be imported.
        :type file_path: str
        """

        nodes_df = pd.read_excel(file_path)
        nodes = [nodes_df["X"].tolist(), nodes_df["Y"].tolist(), nodes_df["Z"].tolist()]
        for i in nodes[0]: self.nodes.append(Node(x=nodes[0][i], y=nodes[1][i], z=nodes[2][i]))

    def import_members(self, file_path: str) -> None:
        """
        Imports member end node indices, material index, if cross-sections properties are to be set,
        and the cross-section properties for the specified Excel file.
        Headers must be "i Node, j Node, Material, Set Cross-Section Properties, A, Iy, Iz, J"
        and there must only be one sheet.

        :param file_path: file path to the Excel file to be imported.
        :type file_path: str
        """

        m_df = pd.read_excel(file_path)
        m = [m_df["i Node"].tolist(), m_df["j Node"].tolist(), m_df["Material"].tolist(),
                   m_df["Set Cross-Section Properties"].tolist(), m_df["A"].tolist(), m_df["Iy"].tolist(),
                   m_df["Iz"].tolist(), m_df["J"].tolist()]
        for i in m[0]: self.members.append(Member(m[0][i], m[1][i], m[2][i], m[3][i],
                                       [m[4][i], m[5][i], m[6][i], m[7][i]], self))

    def import_materials(self, file_path: str) -> None:
        """
        Imports material properties from the specified Excel file.
        Headers must be "E, G, nu, rho, fy" and there must only be one sheet.

        :param file_path: file path to the Excel file to be imported.
        :type file_path: str
        """

        m_df = pd.read_excel(file_path)
        m = [m_df["E"].tolist(), m_df["G"].tolist(), m_df["nu"].tolist(), m_df["rho"].tolist(),
                     m_df["fy"].tolist()]
        for i in m[0]: self.materials.append(Material(m[0][i], m[1][i], m[2][i], m[3][i], m[4][i]))

    def import_supports(self, file_path):
        pass

    def import_releases(self, file_path):
        pass

    def import_node_loads(self, file_path):
        pass

    def import_member_point_loads(self, file_path):
        pass

    def import_member_dist_loads(self, file_path):
        pass

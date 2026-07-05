"""
Handels all operations to do with the 3D frame.
"""
import numpy as np
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
    def __init__(self, controller, display_frame):

        self._controller = controller
        self._display_frame = display_frame

        self.nodes: list[Node] = []
        self.materials: list[Material] = []
        self.members: list[Member] = []

    def updateDisplays(self):
        """
        Updates all places that display info about the truss.
        """

        # Updating Adding Tables if the window is open
        if self._display_frame.TableWindow is not None:
            tables = self._display_frame.TableWindow.tables

            # nodes
            tables[0].delete(*tables[0].get_children())
            tables[3].delete(*tables[3].get_children())
            tables[5].delete(*tables[5].get_children())
            j = 0
            k = 0
            for i in range(len(self.nodes)):
                tables[0].insert('', 'end', values=[str(i), self.nodes[i].x, self.nodes[i].y, self.nodes[i].z])
                if any(self.nodes[i].support):
                    tables[3].insert('', 'end', values=[str(j), str(i)] + self.nodes[i].support)
                    j += 1
                for a in range(len(self.nodes[i].load)):
                    tables[5].insert('', 'end', values=[str(k), str(i)] + self.nodes[i].load[a])
                    k += 1

            # materials
            tables[1].delete(*tables[1].get_children())
            for i in range(len(self.materials)):
                tables[1].insert('', 'end', values=[str(i), self.materials[i].E, self.materials[i].G,
                                                    self.materials[i].nu, self.materials[i].rho, self.materials[i].fy])

            # members
            tables[2].delete(*tables[2].get_children())
            tables[4].delete(*tables[4].get_children())
            tables[6].delete(*tables[6].get_children())
            tables[7].delete(*tables[7].get_children())
            j = 0
            k = 0
            c = 0
            for i in range(len(self.members)):
                tables[2].insert('', 'end', values=[str(i), self.members[i].node_1_index, self.members[i].node_2_index,
                                                    self.members[i].material_index,
                                                    self.members[i].set_cross_section_props, self.members[i].A,
                                                    self.members[i].Iy, self.members[i].Iz, self.members[i].J])
                if any(self.members[i].releces):
                    tables[4].insert('', 'end', values=[str(j), str(i)] + self.members[i].releces)
                    j += 1
                for a in range(len(self.members[i].point_loads)):
                    tables[5].insert('', 'end', values=[str(k), str(i)] + self.members[i].point_loads[a])
                    k += 1
                for b in range(len(self.members[i].dist_loads)):
                    tables[7].insert('', 'end', values=[str(c), str(i)] + self.members[i].dist_loads[b])
                    c += 1

        # TODO

    def linear_analysis(self):
        pass #TODO

    def GlobalOptimization(self, group_assignments, group_types, lower_bounds, upper_bounds,
                                                      cost_function, weight_run, reaction_run, internal_forces_run):
        # tempory testing
        bounds = [(-5, 5), (-5, 5)]
        start_optimization(self._controller.executor, self._controller.root, bounds)

        pass #TODO

    def save(self, file_path):
        pass #TODO

    def export_results(self, file_path):
        pass #TODO

    def open_frame(self, file_path: str) -> None:
        """
        Opens a .structframe file from the provided file path.

        :param file_path: Path to the .structframe file
        :type file_path: str
        """

        with open(file_path, "r") as f:
            lines = f.readlines()
            f.close()

            # clearing current file
            self.nodes = []
            self.materials = []
            self.members = []

            # Nodes
            num_nodes = int(lines[0])
            for i in range(num_nodes):
                node_line = lines[i + 1].replace(" ", "")
                n = []
                for x in node_line.split(','):
                    try:
                        n.append(float(x))
                    except ValueError:
                        pass
                self.nodes.append(Node(cords=n))

            # Materials
            index = num_nodes + 1
            num_materials = int(lines[index])
            for i in range(num_materials):
                material_line = lines[i + index + 1].replace(" ", "")
                n = []
                for x in material_line.split(','):
                    try:
                        n.append(float(x))
                    except ValueError:
                        pass
                self.materials.append(Material(n[0],n[1],n[2],n[3],n[4]))

            # Members
            index += num_materials + 1
            num_members = int(lines[index])
            for i in range(num_members):
                member_line = lines[i + index + 1].replace(" ", "")
                n = []
                for x in member_line.split(','):
                    try:
                        if x == "True":
                            n.append(True)
                        elif x == "False":
                            n.append(False)
                        else:
                            n.append(float(x))
                    except ValueError:
                        pass
                self.members.append(Member(int(n[0]),int(n[1]),int(n[2]),n[3],[n[4],n[5],n[6],n[7]],
                                           self))

            # Supports
            index += num_members + 1
            num_supports = int(lines[index])
            for i in range(num_supports):
                support_line = lines[i + index + 1].replace(" ", "")
                n = []
                for x in support_line.split(','):
                    try:
                        if x == "True":
                            n.append(True)
                        elif x == "False":
                            n.append(False)
                        else:
                            n.append(float(x))
                    except ValueError:
                        pass
                self.nodes[int(n[0])].set_support([n[1],n[2],n[3],n[4],n[5],n[6]])

            # Releases
            index += num_supports + 1
            num_releases = int(lines[index])
            for i in range(num_releases):
                release_line = lines[i + index + 1].replace(" ", "")
                n = []
                for x in release_line.split(','):
                    try:
                        if x == "True":
                            n.append(True)
                        elif x == "False":
                            n.append(False)
                        else:
                            n.append(float(x))
                    except ValueError:
                        pass
                self.members[int(n[0])].set_releces([n[1],n[2],n[3],n[4],n[5],n[6],n[7],n[8],n[9],n[10],n[11],n[12]])

            # Node Load
            index += num_releases + 1
            num_node_loads = int(lines[index])
            for i in range(num_node_loads):
                node_load_line = lines[i + index + 1].replace(" ", "")
                n = []
                for x in node_load_line.split(','):
                    try:
                        if x == "True":
                            n.append(True)
                        elif x == "False":
                            n.append(False)
                        else:
                            n.append(float(x))
                    except ValueError:
                        pass
                self.nodes[int(n[0])].add_load(int(n[7]),[n[1],n[2],n[3],n[4],n[5],n[6]])

            # Member Point Load
            index += num_node_loads + 1
            num_member_point_loads = int(lines[index])
            for i in range(num_member_point_loads):
                member_point_line = lines[i + index + 1].replace(" ", "")
                n = []
                for x in member_point_line.split(','):
                    try:
                        if x == "True":
                            n.append(True)
                        elif x == "False":
                            n.append(False)
                        else:
                            n.append(float(x))
                    except ValueError:
                        pass
                self.members[int(n[0])].add_point_load(int(n[8]),n[1],[n[2],n[3],n[4],n[5],n[6],n[7]])

            # Member Distributed Load
            index += num_member_point_loads + 1
            num_member_dist_loads = int(lines[index])
            for i in range(num_member_dist_loads):
                member_dist_line = lines[i + index + 1].replace(" ", "")
                n = []
                for x in member_dist_line.split(','):
                    try:
                        if x == "True":
                            n.append(True)
                        elif x == "False":
                            n.append(False)
                        else:
                            n.append(float(x))
                    except ValueError:
                        pass
                self.members[int(n[0])].add_dist_load(int(n[9]),[n[1],n[2]],[n[3],n[4],n[5],n[6],n[7],n[8]])

        self.updateDisplays()

    def open_results(self, file_path):
        pass #TODO

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

        self.updateDisplays()

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

        self.updateDisplays()

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

        self.updateDisplays()

    def import_supports(self, file_path: str) -> None:
        """
        Imports supports node indices, and what degrees of freedom are to be supported for the specified Excel file.
        Headers must be "Node, DX, DY, DZ, RX, RY, RZ" and there must only be one sheet.

        :param file_path: file path to the Excel file to be imported.
        :type file_path: str
        """

        s_df = pd.read_excel(file_path)
        supports = [s_df["Node"].tolist(), s_df["DX"].tolist(), s_df["DY"].tolist(), s_df["DZ"].tolist(),
                    s_df["RX"].tolist(), s_df["RY"].tolist(), s_df["RZ"].tolist()]
        supports = np.array(supports).T
        for i in range(len(supports)):
            self.nodes[supports[i,0]].set_support(supports[i,1:].tolist())

        self.updateDisplays()

    def import_releases(self, file_path: str) -> None:
        """
        Import released member indices, and what ends and degrees of freedom are to be released for the specified
        Excel file.
        Headers must be "Member, i DX, i DY, i DZ, i RX, i RY, i RZ, j DX, j DY, j DZ, j RX, j RY, j RZ" and there
        must only be one sheet.

        :param file_path: file path to the Excel file to be imported.
        :type file_path: str
        """

        r_df = pd.read_excel(file_path)
        releases = [r_df["Member"].tolist(), r_df["i DX"].tolist(), r_df["i DY"].tolist(), r_df["i DZ"].tolist(),
                    r_df["i RX"].tolist(), r_df["i RY"].tolist(), r_df["i RZ"].tolist(), r_df["j DX"].tolist(),
                    r_df["j DY"].tolist(), r_df["j DZ"].tolist(), r_df["j RX"].tolist(), r_df["j RY"].tolist(),
                    r_df["j RZ"].tolist()]
        releases = np.array(releases).T
        for i in range(len(releases)):
            self.members[releases[i,0]].set_releces(releases[i,1:].tolist())

        self.updateDisplays()

    def import_node_loads(self, file_path: str) -> None:
        """
        Imports node load, node indices, what directions and magnitudes the load is to be in,
        and the index of the load case.
        Headers must be "Node, PX, PY, PZ, MX, MY, MZ, Case" and there must only be one sheet.

        :param file_path: file path to the Excel file to be imported.
        :type file_path: str
        """

        n_df = pd.read_excel(file_path)
        node_loads = [n_df["Node"].tolist(), n_df["Case"].tolist(), n_df["PX"].tolist(), n_df["PY"].tolist(),
                      n_df["PZ"].tolist(),n_df["MX"].tolist(), n_df["MY"].tolist(), n_df["MZ"].tolist()]
        node_loads = np.array(node_loads).T
        for i in range(len(node_loads)):
            self.nodes[node_loads[i,0]].add_load(node_loads[i,1].tolist(),node_loads[i,2:].tolist())

        self.updateDisplays()

    def import_member_point_loads(self, file_path: str) -> None:
        """
        Imports member point loads, member indices, location, what directions and magnitudes the load is to be in
        and the index of the load case.
        Headers must be "Member, X, PX, PY, PZ, MX, MY, MZ, Case" and there must only be one sheet.
        """

        m_df = pd.read_excel(file_path)
        member_point_loads = [m_df["Member"].tolist(), m_df["Case"].tolist(), m_df["X"].tolist(), m_df["PX"].tolist(),
                              m_df["PY"].tolist(),m_df["PZ"].tolist(), m_df["MX"].tolist(), m_df["MY"].tolist(),
                              m_df["MZ"].tolist()]
        member_point_loads = np.array(member_point_loads).T
        for i in range(len(member_point_loads)):
            self.members[member_point_loads[i,0]].add_point_load(member_point_loads[i,1].tolist(),
                                                                 member_point_loads[i,2].tolist(),
                                                                 member_point_loads[i,3:].tolist())

        self.updateDisplays()

    def import_member_dist_loads(self, file_path: str) -> None:
        """
        Imports member distributed loads: member indices, start and end locations, start and end force magnitude and
        directions, and the index of the load case
        Headers must be "Member, X1, X2, WX1, WX2, WY1, WY2, WZ1, WZ2, Case" and there must only be one sheet.
        """

        m_df = pd.read_excel(file_path)
        member_dist_loads = [m_df["Member"].tolist(), m_df["Case"].tolist(), m_df["X1"].tolist(), m_df["X2"].tolist(),
                             m_df["WX1"].tolist(), m_df["WX2"].tolist(), m_df["WY1"].tolist(), m_df["WY2"].tolist(),
                             m_df["WZ1"].tolist(), m_df["WZ2"].tolist()]
        member_dist_loads = np.array(member_dist_loads).T
        for i in range(len(member_dist_loads)):
            self.members[member_dist_loads[i,0]].add_dist_load(member_dist_loads[i,1].tolist(),
                                                               member_dist_loads[i,2:3].tolist(),
                                                               member_dist_loads[i,4:].tolist())

        self.updateDisplays()
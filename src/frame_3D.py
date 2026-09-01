"""
Handels all operations to do with the 3D frame.
"""
import numpy as np
import pandas as pd
import scipy.optimize as opt

from data_objects.member import Member
from data_objects.material import Material
from data_objects.node import Node
from src.excel_analysis import results_frame_3d_opt
from optimization import start_3D_frame_global_optimization, start_3D_frame_local_optimization
import frame_3D_solver.main
import frame_3D_solver.helper_functions as hf

__author__ = "Blake McCulligh"
__copyright__ = "Copyright 2026 Blake McCulligh"
__credits__ = ["Blake McCulligh"]

__license__ = "MIT"
__version__ = "0.1.0b1"
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = "Beta"

class Frame3D:
    def __init__(self, controller, display_frame):

        self._controller = controller
        self._display_frame = display_frame

        self.nodes: list[Node] = []
        self.materials: list[Material] = []
        self.members: list[Member] = []

        self.solver_frame = None
        self.import_data = None

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
                    tables[6].insert('', 'end', values=[str(k), str(i)] + self.members[i].point_loads[a])
                    k += 1
                for b in range(len(self.members[i].dist_loads)):
                    tables[7].insert('', 'end', values=[str(c), str(i)] + self.members[i].dist_loads[b])
                    c += 1

        # TODO update 3D display system

    def linear_analysis(self):
        """
        Runs linear analysis on the frame.
        """

        self.load_data_solver()

        self.solver_frame.PreAnalysisLinear()
        D, DX, DY, DZ, RX, RY, RZ, weight, reactions, internal_forces = self.solver_frame.AnalysisLinear(
                                                                                             get_weight=True,
                                                                                             get_reactions=True,
                                                                                             get_internal_forces=True)
        #TODO handeling results

    def global_optimization(self, group_assignments: list[int], group_types: list[str], lower_bounds: list[float],
                            upper_bounds: list[float], cost_function: str, weight_run: bool, reaction_run: bool,
                            internal_forces_run: bool, import_data: list | None = None) -> None:
        """
        Runs global optimization on the cross-sections of unset members frame.

        :param group_assignments: list of indices of member groups for non set members to be assigned to.
                                  Must be length of non set members.
        :type group_assignments: list[int]
        :param group_types: List of cross-section types for each member group to be assigned. Must be length of number
                            of member groups.
        :type group_types: list[str]
        :param lower_bounds: Lower bound on the optimization variables. If list must be length of number of variables.
        :type lower_bounds: list[float]
        :param upper_bounds: Upper bound on the optimization variables. If list must be length of number of variables.
        :type upper_bounds: list[float]
        :param cost_function: Cost function for optimization stored in a string.
        :type cost_function: str
        :param weight_run: Weather the weight of all the members is needed for the cost function.
        :type weight_run: bool
        :param reaction_run: Weather the reactions are needed for the cost function.
        :type reaction_run: bool
        :param internal_forces_run: Weather the internal forces are needed for the cost function.
        :type internal_forces_run: bool
        :param import_data: List of the data from an Excel file analysis
        :type import_data: list
        """

        self.load_data_solver()
        self.import_data = import_data
        start_3D_frame_global_optimization(self._controller.executor, self, self.solver_frame,
                                           group_assignments, group_types, lower_bounds, upper_bounds, cost_function,
                                           weight_run, reaction_run, internal_forces_run)

    def local_optimization(self, group_assignments: list[int], group_types: list[str], lower_bounds: list[float],
                           upper_bounds: list[float], inital: list[float], cost_function: str, weight_run: bool,
                           reaction_run: bool, internal_forces_run: bool, import_data: list | None = None):
        """
        Runs local optimization on the cross-sections of unset members frame.

        :param group_assignments: list of indices of member groups for non set members to be assigned to.
                                    Must be length of non set members.
        :type group_assignments: list[int]
        :param group_types: List of cross-section types for each member group to be assigned. Must be length of number
                            of member groups.
        :type group_types: list[str]
        :param lower_bounds: Lower bound on the optimization variables. If list must be length of number of variables.
        :type lower_bounds: list[float]
        :param upper_bounds: Upper bound on the optimization variables. If list must be length of number of variables.
        :type upper_bounds: list[float]
        :param inital: Initial guess for the optimization variables.
        :type inital: list[float]
        :param cost_function: Cost function for optimization stored in a string.
        :type cost_function: str
        :param weight_run: Weather the weight of all the members is needed for the cost function.
        :type weight_run: bool
        :param reaction_run: Weather the reactions are needed for the cost function.
        :type reaction_run: bool
        :param internal_forces_run: Weather the internal forces are needed for the cost function.
        :type internal_forces_run: bool
        :param import_data: List of the data from an Excel file analysis
        :type import_data: list
        """

        self.load_data_solver()
        self.import_data = import_data
        start_3D_frame_local_optimization(self._controller.executor, self, self.solver_frame,
                                          group_assignments, group_types, lower_bounds, upper_bounds, inital,
                                          cost_function, weight_run, reaction_run, internal_forces_run)

    def handel_optimization_results(self, result_opt: opt.OptimizeResult, group_assignments: list[int],
                                    group_types: list[str], optimization_type: str) -> None:
        """
        Updates the frame to the optimized results and exports the optimization results to an Excel file.

        :param result_opt: scipy's OptimizeResult object
        :type result_opt: scipy.optimize.OptimizeResult
        :param group_assignments: list of indices of member groups for non set members to be assigned to.
                                  Must be length of non set members.
        :type group_assignments: list[int]
        :param group_types: List of cross-section types for each member group to be assigned. Must be length of number
                            of member groups.
        :type group_types: list[str]
        :param optimization_type: Type of optimization that was performed.
        :type optimization_type: str
        """

        x = result_opt.x
        cost = result_opt.fun

        if x is not None:
            cross_section_props = hf.get_cross_section_props(x, group_assignments, group_types)
            j = 0
            member_indices = []
            for i in range(len(self.members)):
                if not self.members[i].set_cross_section_props:
                    member_indices.append(i)
                    cs = cross_section_props[j]
                    self.members[i].set_cross_section_props = True
                    self.members[i].A = cs[0]
                    self.members[i].Iy = cs[1]
                    self.members[i].Iz = cs[2]
                    self.members[i].J = cs[3]
                    j += 1

            self.linear_analysis()

            groupLength = []
            startLocation = []

            for i in range(len(group_types)):
                group_type = group_types[i]
                startLocation.append(sum(groupLength) - 1)
                if group_type == ["Angle"]:
                    groupLength.append(3)
                elif group_type == ["RectHSS"]:
                    groupLength.append(3)
                elif group_type == ["SquareHSS"]:
                    groupLength.append(2)
                else:
                    groupLength.append(2)

            results = []
            current_x_index = 0
            for i in range(len(cross_section_props)):
                member_index = [member_indices[i]]
                group_type = [group_types[group_assignments[i]]]

                dim = []
                if group_type == ["Angle"]:
                    dim.append(x[startLocation[group_assignments[i]]])
                    dim.append(x[startLocation[group_assignments[i]] + 1])
                    dim.append(x[startLocation[group_assignments[i]] + 2])
                    current_x_index += 3
                elif group_type == ["RectHSS"]:
                    dim.append(x[startLocation[group_assignments[i]]])
                    dim.append(x[startLocation[group_assignments[i]] + 1])
                    dim.append(x[startLocation[group_assignments[i]] + 2])
                    current_x_index += 3
                elif group_type == ["SquareHSS"]:
                    dim.append(x[startLocation[group_assignments[i]]])
                    dim.append("N/A")
                    dim.append(x[startLocation[group_assignments[i]] + 1])
                    current_x_index += 2
                elif group_type == ["TubeHSS"]:
                    dim.append(x[startLocation[group_assignments[i]]])
                    dim.append("N/A")
                    dim.append(x[startLocation[group_assignments[i]] + 1])
                    current_x_index += 2

                cs = cross_section_props[i]

                results.append(member_index + group_type + dim + cs)

            if self._display_frame is None:
                results_frame_3d_opt(cost, results, self.import_data, optimization_type)

            else:

                self.updateDisplays()

                file_path = self._display_frame.get_new_file_path(".xlsx", [("Excel Files", "*.xlsx")])

                if file_path is not None:

                    df1 = pd.DataFrame(results,
                                       columns=["Member Index", "Cross-Section _type", "d", "b", "t", "A", "Iy", "Iz","J"])
                    df2 = pd.DataFrame({'Cost': [cost], 'Optimization Type': [optimization_type]})

                    with pd.ExcelWriter(file_path) as writer:
                        df1.to_excel(writer, sheet_name="Results", index=False)
                        df2.to_excel(writer, sheet_name="Cost", index=False)

                else:
                    print("Failed to export optimization Results")

    def save(self, file_path):
        pass #TODO saving

    def export_results(self, file_path):
        pass #TODO exporting resutls

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
                                                               member_dist_loads[i,2:4].tolist(),
                                                               member_dist_loads[i,4:].tolist())

        self.updateDisplays()

    def load_data_solver(self):
        """
        Creates a 3D frame solver object from the data stored in this frame.
        """

        self.solver_frame = frame_3D_solver.main.Frame3D_Solver()

        for n in self.nodes:
            self.solver_frame.AddNode(n.x, n.y, n.z)

        for m in self.materials:
            self.solver_frame.AddMaterial(m.E, m.G, m.nu, m.rho, m.fy)

        for m in self.members:
            self.solver_frame.AddMember(m.node_1_index, m.node_2_index, m.material_index, m.set_cross_section_props,
                            m.A, m.Iy, m.Iz, m.J)

        for i, n in enumerate(self.nodes):
            self.solver_frame.AddSupport(i, n.support[0], n.support[1], n.support[2], n.support[3], n.support[4],
                                         n.support[5])

        for i, m in enumerate(self.members):
            self.solver_frame.AddReleases(i, m.releces[0], m.releces[1], m.releces[2], m.releces[3], m.releces[4],
                                          m.releces[5], m.releces[6], m.releces[7], m.releces[8], m.releces[9],
                                          m.releces[10], m.releces[11])

        for i, n in enumerate(self.nodes):
            for l in n.load:
                self.solver_frame.AddNodeLoad(i, l[0], l[1], l[2], l[3], l[4], l[5], l[6])

        for i, m in enumerate(self.members):
            for l in m.point_loads:
                self.solver_frame.AddMemberPointLoad(i, l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7])

        for i, m in enumerate(self.members):
            for l in m.dist_loads:
                self.solver_frame.AddMemberDistLoad(i, l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8])
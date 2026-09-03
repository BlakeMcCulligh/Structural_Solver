"""
Handels all operations to do with the 3D frame.
"""

import scipy.optimize as opt

from data_objects.member import Member
from data_objects.material import Material
from data_objects.node import Node
from src.excel_analysis import results_frame_2d_opt
from optimization import start_2D_frame_global_optimization, start_2D_frame_local_optimization
import frame_2D_solver.main
from frame_3D_solver.helper_functions import get_cross_section_props

__author__ = "Blake McCulligh"
__copyright__ = "Copyright 2026 Blake McCulligh"
__credits__ = ["Blake McCulligh"]

__license__ = "MIT"
__version__ = "0.1.0b1"
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = "Beta"

class Frame2D:
    def __init__(self, controller, display_frame):

        self._controller = controller
        self._display_frame = display_frame

        self.nodes: list[Node] = []
        self.materials: list[Material] = []
        self.members: list[Member] = []

        self.solver_frame = None
        self.import_data = None

    def linear_analysis(self):
        """
        Runs linear analysis on the frame.
        """

        self.load_data_solver()

        self.solver_frame.PreAnalysisLinear()
        D, DX, DY, R, weight, reactions = self.solver_frame.AnalysisLinear(get_weight=True,get_reactions=True)
        #TODO handeling results

    def global_optimization(self, group_assignments: list[int], group_types: list[str], lower_bounds: list[float],
                            upper_bounds: list[float], cost_function: str, weight_run: bool, reaction_run: bool,
                            import_data: list | None = None) -> None:
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
        :param import_data: List of the data from an Excel file analysis
        :type import_data: list
        """

        self.load_data_solver()
        self.import_data = import_data
        start_2D_frame_global_optimization(self._controller.executor, self, self.solver_frame,
                                           group_assignments, group_types, lower_bounds, upper_bounds, cost_function,
                                           weight_run, reaction_run)

    def local_optimization(self, group_assignments: list[int], group_types: list[str], lower_bounds: list[float],
                           upper_bounds: list[float], inital: list[float], cost_function: str, weight_run: bool,
                           reaction_run: bool, import_data: list | None = None):
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
        :param import_data: List of the data from an Excel file analysis
        :type import_data: list
        """

        self.load_data_solver()
        self.import_data = import_data
        start_2D_frame_local_optimization(self._controller.executor, self, self.solver_frame,
                                          group_assignments, group_types, lower_bounds, upper_bounds, inital,
                                          cost_function, weight_run, reaction_run)

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
            cross_section_props = get_cross_section_props(x, group_assignments, group_types)
            j = 0
            member_indices = []
            for i in range(len(self.members)):
                if not self.members[i].set_cross_section_props:
                    member_indices.append(i)
                    cs = cross_section_props[j]
                    self.members[i].set_cross_section_props = True
                    self.members[i].A = cs[0]
                    self.members[i].Iy = 0
                    self.members[i].Iz = cs[1]
                    self.members[i].J = 0
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

            results_frame_2d_opt(cost, results, self.import_data, optimization_type)

    def load_data_solver(self):
        """
        Creates a 3D frame solver object from the data stored in this frame.
        """

        self.solver_frame = frame_2D_solver.main.Frame2D_Solver()

        for n in self.nodes:
            self.solver_frame.AddNode(n.x, n.y)

        for m in self.materials:
            self.solver_frame.AddMaterial(m.E, m.G, m.nu, m.rho, m.fy)

        for m in self.members:
            self.solver_frame.AddMember(m.node_1_index, m.node_2_index, m.material_index, m.set_cross_section_props,
                            m.A, m.Iz)

        for i, n in enumerate(self.nodes):
            self.solver_frame.AddSupport(i, n.support[0], n.support[1], n.support[2])

        for i, m in enumerate(self.members):
            self.solver_frame.AddReleases(i, m.releces[0], m.releces[1], m.releces[2], m.releces[3], m.releces[4],
                                          m.releces[5])

        for i, n in enumerate(self.nodes):
            for l in n.load:
                self.solver_frame.AddNodeLoad(i, l[0], l[1], l[2], l[3])

        for i, m in enumerate(self.members):
            for l in m.point_loads:
                self.solver_frame.AddMemberPointLoad(i, l[0], l[1], l[2], l[3], l[4])

        for i, m in enumerate(self.members):
            for l in m.dist_loads:
                self.solver_frame.AddMemberDistLoad(i, l[0], l[1], l[2], l[3], l[4], l[5], l[6])
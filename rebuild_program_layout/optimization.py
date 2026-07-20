"""
Handels the optimization of the structure.
"""

from __future__ import annotations
from typing import TYPE_CHECKING
from concurrent.futures.thread import ThreadPoolExecutor

if TYPE_CHECKING:
    from rebuild_program_layout.frame_3D import Frame3D
    from frame_3D_solver.main import Frame3D_Solver

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

def start_global_optimization(executor: ThreadPoolExecutor, frame3D: Frame3D, solver_frame: Frame3D_Solver,
                              group_assignments: list[int], group_types: list[str], lower_bounds: list[float],
                              upper_bounds: list[float], cost_function: str, weight_run: bool, reaction_run: bool,
                              internal_forces_run: bool):
    """
    Runs the global optimization of the non set cross-sections of the frame.

    :param executor: Runs functions in a diferent thread then the rest of the program.
    :type executor: concurrent.futures.thread.ThreadPoolExecutor
    :param frame3D: 3D frame that is to be optimized.
    :type frame3D: frame_3D.Frame3D
    :param solver_frame: Solver object for the 3D frame that is to be optimized.
    :type solver_frame: frame_3D_solver.main.Frame3D_Solver
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
    """

    global_opt = executor.submit(lambda: solver_frame.optimize(group_assignments, group_types, lower_bounds,
                                                               upper_bounds, cost_function, weight_run, reaction_run,
                                                               internal_forces_run))

    global_opt.done()
    frame3D.handel_global_optimization_results(global_opt.result(), group_assignments, group_types)

    error = global_opt.exception()
    if error:
        print(f"The thread failed with error type: {type(error).__name__}")
        print(f"Error message: {error}")

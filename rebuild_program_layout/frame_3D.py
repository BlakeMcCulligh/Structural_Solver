"""
Handels all operations to do with the 3D frame.
"""

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

    def import_nodes(self, file_path):
        pass

    def import_members(self, file_path):
        pass

    def import_materials(self, file_path):
        pass

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

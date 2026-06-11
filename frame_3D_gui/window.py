"""
Holds the class that handles everything to do with 3D frames gui
including their windows, inputting, exporting, and saving.
"""

import tkinter as tk
from tkinter import ttk,  filedialog
from copy import copy
import pandas as pd
import numpy as np
from scipy.optimize import OptimizeResult

from frame_3D_gui import export
from frame_3D_gui.data import Data
from frame_3D_gui.export import export_results
from frame_3D_gui.opening import open_frame, open_results
from frame_3D_gui.optimize_pop_up import OptimizationPopUp
from frame_3D_gui.results import Results
from frame_3D_gui.save import save_frame, save_results
from frame_3D_gui.text_validate import validate_float, validate_index, validate_bool
from frame_3D_solver.main import Frame3D
# noinspection PyPep8Naming
import drawing_3D.engine_3D as TDE
import frame_3D_solver.helper_functions as hf
from frame_3D_gui.display import Display

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

class MainWindow(tk.Frame):
    """
    Main window Object for 3D Frame Analyses and anything else to do with 3D Frames.
    """

    def __init__(self, root):
        """
        Constructor for the main window of the 3D frame.

        :param root: Root of the window.
        """

        tk.Frame.__init__(self, root)

        self.FilePath: str | None = None # File path for saving.

        self.Root = root

        self._center_window()

        self.box_set = None
        self.popup_settings = None

        self._create_top_menu()

        self.Data: Data = Data()

        self.val_float = (self.Root.register(validate_float), '%P')
        self.val_index = (self.Root.register(validate_index), '%P')
        self.val_bool = (self.Root.register(validate_bool), '#P')

        # setting analysis objects to none.
        self.Frame: Frame3D | None = None
        self.Results: Results | None = None
        self.OptimizationResults: OptimizeResult | None = None

        # 3D rendering canvas
        self.graph = tk.Canvas(root, bg="white")
        self.graph.pack(fill="both", expand=True)

        # 3D camera Data
        self.Camera: np.ndarray = np.array([0.0, -2.0, 0.0])
        self.Up: np.ndarray = np.array([0.0, 0.0, 1.0])
        self.LookDir: np.ndarray = np.array([0.0, 1.0, 0.0])
        self.CamYRot: float = 0
        self.CamXRot: float = -np.pi/2
        self.Target: np.ndarray = np.array([0.0, -1.0, 0.0])
        self.FOV: float = 90
        self.Z_FAR: float = 1000
        self.Z_NEAR: float = 0.1
        self.move_speed: float = 0.01

        self.LIGHT_DIR: np.ndarray = np.array([0, 0, -1]) # 3D rendering lighting direction unit vector

        # arrays of geometry being displayed in the 3D rendering
        self.DisplayData: Display = Display(self)

        # binding 3D rendering movement inputs
        self.Root.bind("<MouseWheel>", self._zoom)
        self._scroll_cords_last = []
        self._shift_state_last = False
        self.Root.bind("<Button-2>", self._scroll_down)
        self.Root.bind("<B2-Motion>", self._scroll_update)
        self.Root.bind("<ButtonRelease-2>", self._scroll_up)

        # Variables for the input table window
        self._table_window = None # Root of the input Tables Window
        self._table_tabs = None # ttk notebook of the tabs
        self.Tables = []
        self._boxes = [] # input boxes
        self._buttons = []

        # roots for the frame on each tab
        self._node_tab = None
        self._mat_tab = None
        self._member_tab = None
        self._support_tab = None
        self._release_tab = None
        self._node_load_tab = None
        self._member_point_load_tab = None
        self._member_dist_load_tab = None

        self._create_table_window()

        self.Root.protocol("WM_DELETE_WINDOW", self._exit)

    def _center_window(self):
        """
        Centers the main window on the screen and sets its dimensions.
        """

        WIDTH = 1000
        HEIGHT = 800

        screen_width = self.Root.winfo_screenwidth()
        screen_height = self.Root.winfo_screenheight()

        if WIDTH > screen_width or HEIGHT > screen_height:
            WIDTH = screen_width * 0.8
            HEIGHT = screen_height * 0.8
        x = (screen_width // 2) - (WIDTH // 2)
        y = (screen_height // 2) - (HEIGHT // 2) - 50

        self.Root.geometry(f"{WIDTH}x{HEIGHT}+{x}+{y}")

    def _create_top_menu(self):
        """
        Creates the top menu of the window.
        """

        menubar = tk.Menu(self.Root)
        self.Root.config(menu=menubar)

        file_menu = tk.Menu(menubar, tearoff=False)
        file_menu.add_command(label='Save', command=self._save)
        file_menu.add_command(label='Save As', command=self._save_as)
        file_menu.add_command(label='Open Frame', command=self._open)
        file_menu.add_command(label='Export Results', command=self._export_results)
        menubar.add_cascade(label="File", menu=file_menu)

        import_menu = tk.Menu(menubar, tearoff=False)
        import_menu.add_command(label='Nodes', command=self._import_nodes)
        import_menu.add_command(label='Members', command=self._import_members)
        import_menu.add_command(label='Materials', command=self._import_materials)
        import_menu.add_command(label='Supports', command=self._import_supports)
        import_menu.add_command(label='Releases', command=self._import_releases)
        import_menu.add_command(label='Node Loads', command=self._import_node_loads)
        import_menu.add_command(label='Members Point Loads', command=self._import_member_point_loads)
        import_menu.add_command(label='Members Distributed Loads', command=self._import_member_dist_loads)
        menubar.add_cascade(label="Import", menu=import_menu)

        add_menu = tk.Menu(menubar, tearoff=False)
        add_menu.add_command(label='Open Add Tabels', command = self._open_table_window)
        menubar.add_cascade(label="Add", menu=add_menu)

        analysis_menu = tk.Menu(menubar, tearoff=False)
        analysis_menu.add_command(label='Linear Analysis', command=self._linear_analysis)
        analysis_menu.add_command(label='Global Optimization', command=self._optimization_window)
        menubar.add_cascade(label="Analysis", menu=analysis_menu)

        display_menu = tk.Menu(menubar, tearoff=False)
        display_menu.add_command(label='Move Speed', command=self._chamge_move_speed)
        display_menu.add_command(label='Support Scale', command=self._change_support_scale)
        menubar.add_cascade(label="Display Settings", menu=display_menu)

    def _chamge_move_speed(self):
        """
        Changes the movement moltiplier to the input value.
        """

        self.popup_settings = tk.Toplevel(self.Root)

        WIDTH = 400
        HEIGHT = 100
        self._center_popup_window(self.popup_settings, WIDTH, HEIGHT)
        self.popup_settings.title("Move Speed")  # Set the title
        self.popup_settings.resizable(False, False)

        title = tk.Label(self.popup_settings, text="Move Speed: ", font=('Helvetica', 12))
        title.place(x=10, y=10)
        self.box_set = tk.Entry(self.popup_settings, validate='key', validatecommand=self.val_float)
        self.box_set.place(x=150, y=10)
        but = (tk.Button(self.popup_settings, text="okay", command=self._set_move_speed))
        but.place(x=150, y=50)

    def _set_move_speed(self):
        """
        Sets the movement moltiplier to the value in the text box and closes the popup window.
        """

        new_val = self.box_set.get()
        self.move_speed = float(new_val)
        self.popup_settings.destroy()

    def _change_support_scale(self):
        """
        Changes the support scale moltiplier to the input value.
        """

        self.popup_settings = tk.Toplevel(self.Root)

        WIDTH = 400
        HEIGHT = 100
        self._center_popup_window(self.popup_settings, WIDTH, HEIGHT)
        self.popup_settings.title("Support Scale")  # Set the title
        self.popup_settings.resizable(False, False)

        title = tk.Label(self.popup_settings, text="Support Scale: ", font=('Helvetica', 12))
        title.place(x=10, y=10)
        self.box_set = tk.Entry(self.popup_settings, validate='key', validatecommand=self.val_float)
        self.box_set.place(x=150, y=10)
        but = (tk.Button(self.popup_settings, text="okay", command=self._set_support_scale))
        but.place(x=150, y=50)

    def _set_support_scale(self):
        """
        Sets the support scale to the value in the text box and closes the popup window.
        """

        new_val = self.box_set.get()
        self.DisplayData.scale[0] = float(new_val)
        self.popup_settings.destroy()
        self.DisplayData.ConvertToPrint()

    def _center_popup_window(self, popup, width, height):
        """
        Centers the pop-up window on the screen and sets its size.

        :param popup: Popup root.
        :param width: Width to set the pop-Up window to.
        :param height: Height to set the pop-Up window to.
        """

        main_window_x = self.Root.winfo_x()
        main_window_y = self.Root.winfo_y()
        main_window_width = self.Root.winfo_width()
        main_window_height = self.Root.winfo_height()

        x = (main_window_width // 2) + main_window_x - (width // 2)
        y = (main_window_height // 2) + main_window_y - (height // 2)

        popup.geometry(f"{width}x{height}+{x}+{y}")

    def _exit(self):
        """
        Closes program
        """
        self.Root.quit()
        self.Root.destroy()

    def _save_as(self):
        """
        Saves files under specified name and location.
        """

        self.FilePath = save_frame(self.Data)

    def _save(self):
        """
        Saves files in current directory. If no current directory is specified, save under specified name and location.
        """

        if self.FilePath is None:
            self._save_as()
        else:
            save_frame(self.Data, self.FilePath)

    def _export_results(self):
        """
        Exports the analysis results to an Excel file.
        """

        export_results(self.Results)

    def _open(self):
        """
        Opens files under specified name and location.
        """

        open_frame(self)
        open_results(self, self.FilePath)

    def _import_nodes(self):
        """
        Imports coordinates for nodes from the specified Excel file.
        Headers must be "X, Y, Z" and there must only be one sheet.
        """

        file_types = [("Excel Files", "*.xlsx")]
        file_path = _select_file_gui(file_types) # getting file path from user
        nodes_df = pd.read_excel(file_path)
        nodes = [nodes_df["X"].tolist(), nodes_df["Y"].tolist(), nodes_df["Z"].tolist()]
        self.Data.AddNodes(self, nodes, True, True)

    def _import_materials(self):
        """
        Imports material properties from the specified Excel file.
        Headers must be "E, G, nu, rho, fy" and there must only be one sheet.
        """

        file_types = [("Excel Files", "*.xlsx")]
        file_path = _select_file_gui(file_types) # getting file path from user
        m_df = pd.read_excel(file_path)
        materials = [m_df["E"].tolist(),m_df["G"].tolist(),m_df["nu"].tolist(),m_df["rho"].tolist(),m_df["fy"].tolist()]
        self.Data.AddMaterials(self, materials, True, True)

    def _import_members(self):
        """
        Imports member end node indices, material index, if cross-sections properties are to be set,
        and the cross-section properties for the specified Excel file.
        Headers must be "i Node, j Node, Material, Set Cross-Section Properties, A, Iy, Iz, J"
        and there must only be one sheet.
        """

        file_types = [("Excel Files", "*.xlsx")]
        file_path = _select_file_gui(file_types) # getting file path from user
        m_df = pd.read_excel(file_path)
        members = [m_df["i Node"].tolist(), m_df["j Node"].tolist(), m_df["Material"].tolist(),
                   m_df["Set Cross-Section Properties"].tolist(), m_df["A"].tolist(), m_df["Iy"].tolist(),
                   m_df["Iz"].tolist(), m_df["J"].tolist()]
        self.Data.AddMembers(self, members, True, True)

    def _import_supports(self):
        """
        Imports supports node indices, and what degrees of freedom are to be supported for the specified Excel file.
        Headers must be "Node, DX, DY, DZ, RX, RY, RZ" and there must only be one sheet.
        """

        file_types = [("Excel Files", "*.xlsx")]
        file_path = _select_file_gui(file_types) # getting file path from user
        s_df = pd.read_excel(file_path)
        supports = [s_df["Node"].tolist(),s_df["DX"].tolist(),s_df["DY"].tolist(),s_df["DZ"].tolist(),
                    s_df["RX"].tolist(),s_df["RY"].tolist(),s_df["RZ"].tolist()]
        self.Data.AddSupports(self, supports, True, True)

    def _import_releases(self):
        """
        Import released member indices, and what ends and degrees of freedom are to be released for the specified
        Excel file.
        Headers must be "Member, i DX, i DY, i DZ, i RX, i RY, i RZ, j DX, j DY, j DZ, j RX, j RY, j RZ" and there
        must only be one sheet.
        """

        file_types = [("Excel Files", "*.xlsx")]
        file_path = _select_file_gui(file_types) # getting file path from user
        r_df = pd.read_excel(file_path)
        releases = [r_df["Member"].tolist(), r_df["i DX"].tolist(), r_df["i DY"].tolist(), r_df["i DZ"].tolist(),
                    r_df["i RX"].tolist(), r_df["i RY"].tolist(), r_df["i RZ"].tolist(),r_df["j DX"].tolist(),
                    r_df["j DY"].tolist(), r_df["j DZ"].tolist(), r_df["j RX"].tolist(), r_df["j RY"].tolist(),
                    r_df["j RZ"].tolist()]
        self.Data.AddReleases(self, releases, True, True)

    def _import_node_loads(self):
        """
        Imports node load, node indices, what directions and magnitudes the load is to be in,
         and the index of the load case.
        Headers must be "Node, PX, PY, PZ, MX, MY, MZ, Case" and there must only be one sheet.
        """

        file_types = [("Excel Files", "*.xlsx")]
        file_path = _select_file_gui(file_types) # getting file path from user
        n_df = pd.read_excel(file_path)
        node_loads = [n_df["Node"].tolist(),n_df["PX"].tolist(),n_df["PY"].tolist(),n_df["PZ"].tolist(),
                     n_df["MX"].tolist(),n_df["MY"].tolist(),n_df["MZ"].tolist(),n_df["Case"].tolist()]
        self.Data.AddNodeLoads(self, node_loads, True, True)

    def _import_member_point_loads(self):
        """
        Imports member point loads, member indices, location, what directions and magnitudes the load is to be in
        and the index of the load case.
        Headers must be "Member, X, PX, PY, PZ, MX, MY, MZ, Case" and there must only be one sheet.
        """

        file_types = [("Excel Files", "*.xlsx")]
        file_path = _select_file_gui(file_types) # getting file path from user
        m_df = pd.read_excel(file_path)
        member_point_loads = [m_df["Member"].tolist(),m_df["X"].tolist(),m_df["PX"].tolist(),m_df["PY"].tolist(),
                            m_df["PZ"].tolist(),m_df["MX"].tolist(),m_df["MY"].tolist(),m_df["MZ"].tolist(),
                            m_df["Case"].tolist()]
        self.Data.AddMemberPointLoads(self, member_point_loads, True, True)

    def _import_member_dist_loads(self):
        """
        Imports member distributed loads: member indices, start and end locations, start and end force magnitude and
        directions, and the index of the load case
        Headers must be "Member, X1, X2, WX1, WX2, WY1, WY2, WZ1, WZ2, Case" and there must only be one sheet.
        """

        file_types = [("Excel Files", "*.xlsx")]
        file_path = _select_file_gui(file_types) # getting file path from user
        m_df = pd.read_excel(file_path)
        member_dist_loads = [m_df["Member"].tolist(),m_df["X1"].tolist(),m_df["X2"].tolist(),m_df["WX1"].tolist(),
                           m_df["WX2"].tolist(),m_df["WY1"].tolist(),m_df["WY2"].tolist(),m_df["WZ1"].tolist(),
                           m_df["WZ2"].tolist(),m_df["Case"].tolist()]
        self.Data.AddMemberDistLoads(self, member_dist_loads, True, True)

    def _linear_analysis(self):
        """
         Runs linear analysis on the defined frame, and saves the results.
        """

        frame = _add_data_to_frame(self, self.Data)

        frame.PreAnalysisLinear()
        D, DX, DY, DZ, RX, RY, RZ, weight, reactions, internal_forces = frame.AnalysisLinear(get_weight= True,
                                                                                             get_reactions= True,
                                                                                             get_internal_forces= True)
        self.Frame = frame
        self.Results = Results()

        if D is not None:
            self.Results.AddNodalDeflections(DX, DY, DZ, RX, RY, RZ)
            self.Results.add_weight(weight)
            self.Results.AddReactions(reactions)
            self.Results.AddInternalForces(internal_forces)
            self._save()
            save_results(self.Results, self.FilePath)
        elif self.FilePath is not None:
            save_results(self.Results, self.FilePath)

    def _optimization_window(self):
        """
        Opens the window to get input data for an optimization and the start the optimization.
        """

        OptimizationPopUp(self, self.Root, self.Data)

    def GlobalOptimization(self, GroupAssignments, GroupTypes, LowerBound, UpperBound, CostFunction, WeightRun,
                           ReactionRun, InternalForcesRun):
        """
        Runs global cross-section optimization on the frames members that do not have specified cross-sections. Saves
        the optimized frame under a specified file name and path, and its linear analysis results.

        :param GroupAssignments: Member groups assignments for members without specified cross-sections.
        :param GroupTypes: What kind of cross-section each member group is.
        :param LowerBound: Lower bounds on optimization variables (dimensions on the cross-sections).
        :param UpperBound: Upper bounds on optimization variables (dimensions on the cross-sections).
        :param CostFunction: String that carries the equation to get the cost of the frame.
                             The optimizer minimizes the value this equation returns.
        :param WeightRun: If the weight of the members are needed for the cost equation.
        :param ReactionRun: If the reactions of the frame are needed for the cost equation.
        :param InternalForcesRun: If the internal forces of the frame are needed for the cost equation.
        """

        save_frame(self.Data)

        frame = _add_data_to_frame(self, self.Data)

        frame.PreAnalysisLinear()
        results = frame.optimize(GroupAssignments, GroupTypes, LowerBound, UpperBound, CostFunction, WeightRun,
                                 ReactionRun, InternalForcesRun)
        self.Frame = frame
        self.OptimizationResults = results

        self.FilePath = None
        member_indices =  self._update_frame_to_optimization_results(GroupAssignments, GroupTypes)
        if member_indices is not None:
            self._linear_analysis() # runs an analysis of the optimized frame to get the full analysis Results
            self._save_optimization_results(member_indices, GroupAssignments, GroupTypes)

    def _update_frame_to_optimization_results(self, group_assignments, group_types):
        """
        Updates the frame's cross-sections to the results found the optimizer.

        :param group_assignments: Member groups assignments for members without specified cross-sections.
        :param group_types: What kind of cross-section each member group is.
        :return: member_indices: Indices of the optimized member cross-sections.
        """

        x = self.OptimizationResults.x
        if x is not None:
            opt_cross_section_props = hf.get_cross_section_props(x, group_assignments, group_types)
            member_indices = []
            j = 0
            for i in range(len(self.Data.Members[0])):
                if not self.Data.Members[3][i]:
                    member_indices.append(i)
                    cs = opt_cross_section_props[j]
                    self.Data.Members[3][i] = True
                    self.Data.Members[4][i] = cs[0]
                    self.Data.Members[5][i] = cs[1]
                    self.Data.Members[6][i] = cs[2]
                    self.Data.Members[7][i] = cs[3]
                    j += 1
        else:
            member_indices = None

        return member_indices

    def _save_optimization_results(self, member_indices, group_assignments, group_types):
        """
        Exports the optimization results to an Excel file.

        :param member_indices: Indices of the optimized member cross-sections.
        :param group_assignments: Member groups assignments for members without specified cross-sections.
        :param group_types: What kind of cross-section each member group is.
        """

        x = self.OptimizationResults.x
        cost = [self.OptimizationResults.fun]
        opt_cross_section_props = hf.get_cross_section_props(x, group_assignments, group_types)

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
        for i in range(len(opt_cross_section_props)):
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

            cs = opt_cross_section_props[i]

            results.append(member_index + group_type + dim + cs)

        export.export_optimization_results(results, cost)

    """ ----------------------------------------------------------------------------------------------"""
    """ ---------------------------------------- INPUT TABLES ----------------------------------------"""
    """ ----------------------------------------------------------------------------------------------"""

    def _open_table_window(self) -> None:
        """
        Opens the window with all the input tables.
        """

        # Reseting window variables
        self._table_window = None
        self._table_tabs = None
        self.Tables = []
        self._boxes = []
        self._buttons = []
        self._node_tab = None
        self._mat_tab = None
        self._member_tab = None
        self._support_tab = None
        self._release_tab = None
        self._node_load_tab = None
        self._member_point_load_tab = None
        self._member_dist_load_tab = None

        self._create_table_window()

        # adding Data to tables
        self.Data.AddNodeToTable(self)
        self.Data.AddMaterialsToTable(self)
        self.Data.AddMemberToTable(self)
        self.Data.AddSupportToTable(self)
        self.Data.AddReleasesToTable(self)
        self.Data.AddNodeLoadsToTables(self)
        self.Data.AddMemberPointLoadToTable(self)
        self.Data.AddMemberDistLoadsToTables(self)

    def _create_table_window(self):
        """
        Creates window with all the input tables.
        """

        # TODO add delete buttons

        self._table_window = tk.Toplevel(self.Root)
        self._table_window.title("Sheets")
        self._table_window.geometry("800x700")
        self._table_window.resizable(False, True)

        self._table_tabs = ttk.Notebook(self._table_window)
        self._table_tabs.pack(expand=True, fill="both")

        self._node_tab = ttk.Frame(self._table_tabs)

        # Nodes
        self._table_tabs.add(self._node_tab, text="Nodes")
        self.Tables.append(self._create_table(self._node_tab, ("Index", "X", "Y", "Z")))
        self._boxes.append([])
        self._boxes[0].append(tk.Entry(self._node_tab, validate='key', validatecommand=self.val_float))
        self._boxes[0].append(tk.Entry(self._node_tab, validate='key', validatecommand=self.val_float))
        self._boxes[0].append(tk.Entry(self._node_tab, validate='key', validatecommand=self.val_float))
        label = tk.Label(self._node_tab, text="X:", font=('Helvetica', 12))
        label.place(x=325, y=297)
        self._boxes[0][0].place(x=350, y=300) # X
        label = tk.Label(self._node_tab, text="Y:", font=('Helvetica', 12))
        label.place(x=325, y=347)
        self._boxes[0][1].place(x=350, y=350) # Y
        label = tk.Label(self._node_tab, text="Z:", font=('Helvetica', 12))
        label.place(x=325, y=397)
        self._boxes[0][2].place(x=350, y=400) # Z
        self._buttons.append([])
        self._buttons[0].append(tk.Button(self._node_tab, text="Edit", command=self._edit_values))
        self._buttons[0].append(tk.Button(self._node_tab, text="Add", command=self._add_values))
        self._buttons[0][0].place(x=550, y=300)
        self._buttons[0][1].place(x=550, y=350)

        # Materials
        self._mat_tab = ttk.Frame(self._table_tabs)
        self._table_tabs.add(self._mat_tab, text="Materials")
        self.Tables.append(self._create_table(self._mat_tab, ("Index", "E", "G", "nu", "rho", "fy")))
        self._boxes.append([])
        self._boxes[1].append(tk.Entry(self._mat_tab, validate='key', validatecommand=self.val_float))
        self._boxes[1].append(tk.Entry(self._mat_tab, validate='key', validatecommand=self.val_float))
        self._boxes[1].append(tk.Entry(self._mat_tab, validate='key', validatecommand=self.val_float))
        self._boxes[1].append(tk.Entry(self._mat_tab, validate='key', validatecommand=self.val_float))
        self._boxes[1].append(tk.Entry(self._mat_tab, validate='key', validatecommand=self.val_float))
        label = tk.Label(self._mat_tab, text="E:", font=('Helvetica', 12))
        label.place(x=125, y=297)
        self._boxes[1][0].place(x=150, y=300)
        label = tk.Label(self._mat_tab, text="G:", font=('Helvetica', 12))
        label.place(x=325, y=297)
        self._boxes[1][1].place(x=350, y=300)
        label = tk.Label(self._mat_tab, text="Poisons Ratio:", font=('Helvetica', 12))
        label.place(x=40, y=347)
        self._boxes[1][2].place(x=150, y=350)
        label = tk.Label(self._mat_tab, text="Density:", font=('Helvetica', 12))
        label.place(x=70, y=397)
        self._boxes[1][3].place(x=150, y=400)
        label = tk.Label(self._mat_tab, text="Yield Strength:", font=('Helvetica', 12))
        label.place(x=30, y=447)
        self._boxes[1][4].place(x=150, y=450)
        self._buttons.append([])
        self._buttons[1].append(tk.Button(self._mat_tab, text="Edit", command=self._edit_values))
        self._buttons[1].append(tk.Button(self._mat_tab, text="Add", command=self._add_values))
        self._buttons[1][0].place(x=550, y=300)
        self._buttons[1][1].place(x=550, y=350)

        # Members
        self._member_tab = ttk.Frame(self._table_tabs)
        self._table_tabs.add(self._member_tab, text="Members")
        self.Tables.append(self._create_table(self._member_tab, ("Index", "i Node", "j Node", "Material Id",
                                                             "Set C.S.", "A", "Iy", "Iz", "J")))
        self._boxes.append([])
        self._boxes[2].append(tk.Entry(self._member_tab, validate='key', validatecommand=self.val_index))
        self._boxes[2].append(tk.Entry(self._member_tab, validate='key', validatecommand=self.val_index))
        self._boxes[2].append(tk.Entry(self._member_tab, validate='key', validatecommand=self.val_index))
        self._boxes[2].append(ttk.Combobox(self._member_tab, values=("True", "False"),state="readonly"))
        self._boxes[2].append(tk.Entry(self._member_tab, validate='key', validatecommand=self.val_float))
        self._boxes[2].append(tk.Entry(self._member_tab, validate='key', validatecommand=self.val_float))
        self._boxes[2].append(tk.Entry(self._member_tab, validate='key', validatecommand=self.val_float))
        self._boxes[2].append(tk.Entry(self._member_tab, validate='key', validatecommand=self.val_float))
        label = tk.Label(self._member_tab, text="i Node:", font=('Helvetica', 12))
        label.place(x=90, y=297)
        self._boxes[2][0].place(x=150, y=300)
        label = tk.Label(self._member_tab, text="j Node:", font=('Helvetica', 12))
        label.place(x=290, y=297)
        self._boxes[2][1].place(x=350, y=300)
        label = tk.Label(self._member_tab, text="Material ID:", font=('Helvetica', 12))
        label.place(x=60, y=347)
        self._boxes[2][2].place(x=150, y=350)
        label = tk.Label(self._member_tab, text="Set C.S.:", font=('Helvetica', 12))
        label.place(x=275, y=347)
        self._boxes[2][3].place(x=350, y=350, width=125)
        label = tk.Label(self._member_tab, text="A:", font=('Helvetica', 12))
        label.place(x=125, y=397)
        self._boxes[2][4].place(x=150, y=400)
        label = tk.Label(self._member_tab, text="Iy:", font=('Helvetica', 12))
        label.place(x=325, y=397)
        self._boxes[2][5].place(x=350, y=400)
        label = tk.Label(self._member_tab, text="Iz:", font=('Helvetica', 12))
        label.place(x=125, y=447)
        self._boxes[2][6].place(x=150, y=450)
        label = tk.Label(self._member_tab, text="J:", font=('Helvetica', 12))
        label.place(x=325, y=447)
        self._boxes[2][7].place(x=350, y=450)
        self._buttons.append([])
        self._buttons[2].append(tk.Button(self._member_tab, text="Edit", command=self._edit_values))
        self._buttons[2].append(tk.Button(self._member_tab, text="Add", command=self._add_values))
        self._buttons[2][0].place(x=550, y=300)
        self._buttons[2][1].place(x=550, y=350)

        # Supports
        self._support_tab = ttk.Frame(self._table_tabs)
        self._table_tabs.add(self._support_tab, text="Supports")
        self.Tables.append(self._create_table(self._support_tab, ("Index", "i Node", "D X", "D Y", "D Z.", "R X",
                                                             "R Y", "R Z")))
        self._boxes.append([])
        self._boxes[3].append(tk.Entry(self._support_tab, validate='key', validatecommand=self.val_index))
        self._boxes[3].append(ttk.Combobox(self._support_tab, values=("True", "False"), state="readonly"))
        self._boxes[3].append(ttk.Combobox(self._support_tab, values=("True", "False"), state="readonly"))
        self._boxes[3].append(ttk.Combobox(self._support_tab, values=("True", "False"), state="readonly"))
        self._boxes[3].append(ttk.Combobox(self._support_tab, values=("True", "False"), state="readonly"))
        self._boxes[3].append(ttk.Combobox(self._support_tab, values=("True", "False"), state="readonly"))
        self._boxes[3].append(ttk.Combobox(self._support_tab, values=("True", "False"), state="readonly"))
        label = tk.Label(self._support_tab, text="i Node:", font=('Helvetica', 12))
        label.place(x=90, y=297)
        self._boxes[3][0].place(x=150, y=300)
        label = tk.Label(self._support_tab, text="D X:", font=('Helvetica', 12))
        label.place(x=312, y=297)
        self._boxes[3][1].place(x=350, y=300, width=125)
        label = tk.Label(self._support_tab, text="D Y:", font=('Helvetica', 12))
        label.place(x=115, y=347)
        self._boxes[3][2].place(x=150, y=350, width=125)
        label = tk.Label(self._support_tab, text="D Z:", font=('Helvetica', 12))
        label.place(x=315, y=347)
        self._boxes[3][3].place(x=350, y=350, width=125)
        label = tk.Label(self._support_tab, text="R X:", font=('Helvetica', 12))
        label.place(x=112, y=397)
        self._boxes[3][4].place(x=150, y=400, width=125)
        label = tk.Label(self._support_tab, text="R Y:", font=('Helvetica', 12))
        label.place(x=315, y=397)
        self._boxes[3][5].place(x=350, y=400, width=125)
        label = tk.Label(self._support_tab, text="R Z:", font=('Helvetica', 12))
        label.place(x=115, y=447)
        self._boxes[3][6].place(x=150, y=450, width=125)
        self._buttons.append([])
        self._buttons[3].append(tk.Button(self._support_tab, text="Edit", command=self._edit_values))
        self._buttons[3].append(tk.Button(self._support_tab, text="Add", command=self._add_values))
        self._buttons[3][0].place(x=550, y=300)
        self._buttons[3][1].place(x=550, y=350)

        # Releases
        self._release_tab = ttk.Frame(self._table_tabs)
        self._table_tabs.add(self._release_tab, text="Releases")
        self.Tables.append(self._create_table(self._release_tab, ("Index", "i Member", "i DX", "i DY", "i DZ",
                                                             "i RX", "i RY", "i RZ", "j DX", "j DY", "j DZ", "j RX",
                                                             "j RY", "j RZ")))
        self._boxes.append([])
        self._boxes[4].append(tk.Entry(self._release_tab, validate='key', validatecommand=self.val_index))
        for i in range(12): self._boxes[4].append(ttk.Combobox(self._release_tab, values=("True", "False"),
                                                               state="readonly"))
        label = tk.Label(self._release_tab, text="i Member:", font=('Helvetica', 12))
        label.place(x=70, y=297)
        self._boxes[4][0].place(x=150, y=300)
        label = tk.Label(self._release_tab, text="i DX:", font=('Helvetica', 12))
        label.place(x=310, y=297)
        self._boxes[4][1].place(x=350, y=300, width=125)
        label = tk.Label(self._release_tab, text="i DY:", font=('Helvetica', 12))
        label.place(x=110, y=347)
        self._boxes[4][2].place(x=150, y=350, width=125)
        label = tk.Label(self._release_tab, text="i DZ:", font=('Helvetica', 12))
        label.place(x=310, y=347)
        self._boxes[4][3].place(x=350, y=350, width=125)
        label = tk.Label(self._release_tab, text="i RX:", font=('Helvetica', 12))
        label.place(x=110, y=397)
        self._boxes[4][4].place(x=150, y=400, width=125)
        label = tk.Label(self._release_tab, text="i RY:", font=('Helvetica', 12))
        label.place(x=310, y=397)
        self._boxes[4][5].place(x=350, y=400, width=125)
        label = tk.Label(self._release_tab, text="i RZ:", font=('Helvetica', 12))
        label.place(x=110, y=447)
        self._boxes[4][6].place(x=150, y=450, width=125)
        label = tk.Label(self._release_tab, text="j DX:", font=('Helvetica', 12))
        label.place(x=310, y=447)
        self._boxes[4][7].place(x=350, y=450, width=125)
        label = tk.Label(self._release_tab, text="j DY:", font=('Helvetica', 12))
        label.place(x=110, y=497)
        self._boxes[4][8].place(x=150, y=500, width=125)
        label = tk.Label(self._release_tab, text="j DZ:", font=('Helvetica', 12))
        label.place(x=310, y=497)
        self._boxes[4][9].place(x=350, y=500, width=125)
        label = tk.Label(self._release_tab, text="j RX:", font=('Helvetica', 12))
        label.place(x=110, y=547)
        self._boxes[4][10].place(x=150, y=550, width=125)
        label = tk.Label(self._release_tab, text="j RY:", font=('Helvetica', 12))
        label.place(x=310, y=547)
        self._boxes[4][11].place(x=350, y=550, width=125)
        label = tk.Label(self._release_tab, text="j RZ:", font=('Helvetica', 12))
        label.place(x=110, y=597)
        self._boxes[4][12].place(x=150, y=600, width=125)
        self._buttons.append([])
        self._buttons[4].append(tk.Button(self._release_tab, text="Edit", command=self._edit_values))
        self._buttons[4].append(tk.Button(self._release_tab, text="Add", command=self._add_values))
        self._buttons[4][0].place(x=550, y=300)
        self._buttons[4][1].place(x=550, y=350)

        # Nodal Loads
        self._node_load_tab = ttk.Frame(self._table_tabs)
        self._table_tabs.add(self._node_load_tab, text="Node Loads")
        self.Tables.append(self._create_table(self._node_load_tab, ("Index", "i Node", "P X", "P Y", "P Z.",
                                                                    "M X","M Y", "M Z","Cases")))
        self._boxes.append([])
        self._boxes[5].append(tk.Entry(self._node_load_tab, validate='key', validatecommand=self.val_index))
        self._boxes[5].append(tk.Entry(self._node_load_tab, validate='key', validatecommand=self.val_float))
        self._boxes[5].append(tk.Entry(self._node_load_tab, validate='key', validatecommand=self.val_float))
        self._boxes[5].append(tk.Entry(self._node_load_tab, validate='key', validatecommand=self.val_float))
        self._boxes[5].append(tk.Entry(self._node_load_tab, validate='key', validatecommand=self.val_float))
        self._boxes[5].append(tk.Entry(self._node_load_tab, validate='key', validatecommand=self.val_float))
        self._boxes[5].append(tk.Entry(self._node_load_tab, validate='key', validatecommand=self.val_float))
        self._boxes[5].append(tk.Entry(self._node_load_tab, validate='key', validatecommand=self.val_index))
        label = tk.Label(self._node_load_tab, text="i Node:", font=('Helvetica', 12))
        label.place(x=90, y=297)
        self._boxes[5][0].place(x=150, y=300)
        label = tk.Label(self._node_load_tab, text="P X:", font=('Helvetica', 12))
        label.place(x=310, y=297)
        self._boxes[5][1].place(x=350, y=300)
        label = tk.Label(self._node_load_tab, text="P Y:", font=('Helvetica', 12))
        label.place(x=110, y=347)
        self._boxes[5][2].place(x=150, y=350)
        label = tk.Label(self._node_load_tab, text="P Z:", font=('Helvetica', 12))
        label.place(x=310, y=347)
        self._boxes[5][3].place(x=350, y=350)
        label = tk.Label(self._node_load_tab, text="M X:", font=('Helvetica', 12))
        label.place(x=110, y=397)
        self._boxes[5][4].place(x=150, y=400)
        label = tk.Label(self._node_load_tab, text="M Y:", font=('Helvetica', 12))
        label.place(x=310, y=397)
        self._boxes[5][5].place(x=350, y=400)
        label = tk.Label(self._node_load_tab, text="M Z:", font=('Helvetica', 12))
        label.place(x=110, y=447)
        self._boxes[5][6].place(x=150, y=450)
        label = tk.Label(self._node_load_tab, text="Cases:", font=('Helvetica', 12))
        label.place(x=290, y=447)
        self._boxes[5][7].place(x=350, y=450)
        self._buttons.append([])
        self._buttons[5].append(tk.Button(self._node_load_tab, text="Edit", command=self._edit_values))
        self._buttons[5].append(tk.Button(self._node_load_tab, text="Add", command=self._add_values))
        self._buttons[5][0].place(x=550, y=300)
        self._buttons[5][1].place(x=550, y=350)

        # Member Point Loads
        self._member_point_load_tab = ttk.Frame(self._table_tabs)
        self._table_tabs.add(self._member_point_load_tab, text="Member Point Loads")
        self.Tables.append(self._create_table(self._member_point_load_tab, ("Index", "i Member", "x", "P X", "P Y",
                                                                      "P Z.", "M X", "M Y", "M Z", "Cases")))
        self._boxes.append([])
        self._boxes[6].append(tk.Entry(self._member_point_load_tab, validate='key', validatecommand=self.val_index))
        self._boxes[6].append(tk.Entry(self._member_point_load_tab, validate='key', validatecommand=self.val_float))
        self._boxes[6].append(tk.Entry(self._member_point_load_tab, validate='key', validatecommand=self.val_float))
        self._boxes[6].append(tk.Entry(self._member_point_load_tab, validate='key', validatecommand=self.val_float))
        self._boxes[6].append(tk.Entry(self._member_point_load_tab, validate='key', validatecommand=self.val_float))
        self._boxes[6].append(tk.Entry(self._member_point_load_tab, validate='key', validatecommand=self.val_float))
        self._boxes[6].append(tk.Entry(self._member_point_load_tab, validate='key', validatecommand=self.val_float))
        self._boxes[6].append(tk.Entry(self._member_point_load_tab, validate='key', validatecommand=self.val_float))
        self._boxes[6].append(tk.Entry(self._member_point_load_tab, validate='key', validatecommand=self.val_index))
        label = tk.Label(self._member_point_load_tab, text="i Member:", font=('Helvetica', 12))
        label.place(x=70, y=297)
        self._boxes[6][0].place(x=150, y=300)
        label = tk.Label(self._member_point_load_tab, text="x:", font=('Helvetica', 12))
        label.place(x=325, y=297)
        self._boxes[6][1].place(x=350, y=300)
        label = tk.Label(self._member_point_load_tab, text="P X:", font=('Helvetica', 12))
        label.place(x=110, y=347)
        self._boxes[6][2].place(x=150, y=350)
        label = tk.Label(self._member_point_load_tab, text="P Y:", font=('Helvetica', 12))
        label.place(x=310, y=347)
        self._boxes[6][3].place(x=350, y=350)
        label = tk.Label(self._member_point_load_tab, text="P Z:", font=('Helvetica', 12))
        label.place(x=110, y=397)
        self._boxes[6][4].place(x=150, y=400)
        label = tk.Label(self._member_point_load_tab, text="M X:", font=('Helvetica', 12))
        label.place(x=310, y=397)
        self._boxes[6][5].place(x=350, y=400)
        label = tk.Label(self._member_point_load_tab, text="M Y:", font=('Helvetica', 12))
        label.place(x=110, y=447)
        self._boxes[6][6].place(x=150, y=450)
        label = tk.Label(self._member_point_load_tab, text="M Z:", font=('Helvetica', 12))
        label.place(x=310, y=447)
        self._boxes[6][7].place(x=350, y=450)
        label = tk.Label(self._member_point_load_tab, text="Cases:", font=('Helvetica', 12))
        label.place(x=70, y=497)
        self._boxes[6][8].place(x=150, y=500)
        self._buttons.append([])
        self._buttons[6].append(tk.Button(self._member_point_load_tab, text="Edit", command=self._edit_values))
        self._buttons[6].append(tk.Button(self._member_point_load_tab, text="Add", command=self._add_values))
        self._buttons[6][0].place(x=550, y=300)
        self._buttons[6][1].place(x=550, y=350)

        # Member Distributed Loads
        self._member_dist_load_tab = ttk.Frame(self._table_tabs)
        self._table_tabs.add(self._member_dist_load_tab, text="Member Distributed Loads")
        self.Tables.append(self._create_table(self._member_dist_load_tab, ("Index", "i Member", "x 1", "x 2", "wx 1",
                                                                    "wx 2", "wy 1", "wy 2", "wz 1", "wz 2", "Case")))
        self._boxes.append([])
        self._boxes[7].append(tk.Entry(self._member_dist_load_tab, validate='key', validatecommand=self.val_index))
        self._boxes[7].append(tk.Entry(self._member_dist_load_tab, validate='key', validatecommand=self.val_float))
        self._boxes[7].append(tk.Entry(self._member_dist_load_tab, validate='key', validatecommand=self.val_float))
        self._boxes[7].append(tk.Entry(self._member_dist_load_tab, validate='key', validatecommand=self.val_float))
        self._boxes[7].append(tk.Entry(self._member_dist_load_tab, validate='key', validatecommand=self.val_float))
        self._boxes[7].append(tk.Entry(self._member_dist_load_tab, validate='key', validatecommand=self.val_float))
        self._boxes[7].append(tk.Entry(self._member_dist_load_tab, validate='key', validatecommand=self.val_float))
        self._boxes[7].append(tk.Entry(self._member_dist_load_tab, validate='key', validatecommand=self.val_float))
        self._boxes[7].append(tk.Entry(self._member_dist_load_tab, validate='key', validatecommand=self.val_float))
        self._boxes[7].append(tk.Entry(self._member_dist_load_tab, validate='key', validatecommand=self.val_index))
        label = tk.Label(self._member_dist_load_tab, text="i Member:", font=('Helvetica', 12))
        label.place(x=70, y=297)
        self._boxes[7][0].place(x=150, y=300)
        label = tk.Label(self._member_dist_load_tab, text="x 1:", font=('Helvetica', 12))
        label.place(x=310, y=297)
        self._boxes[7][1].place(x=350, y=300)
        label = tk.Label(self._member_dist_load_tab, text="x 2:", font=('Helvetica', 12))
        label.place(x=110, y=347)
        self._boxes[7][2].place(x=150, y=350)
        label = tk.Label(self._member_dist_load_tab, text="wx 1:", font=('Helvetica', 12))
        label.place(x=300, y=347)
        self._boxes[7][3].place(x=350, y=350)
        label = tk.Label(self._member_dist_load_tab, text="wx 2:", font=('Helvetica', 12))
        label.place(x=100, y=397)
        self._boxes[7][4].place(x=150, y=400)
        label = tk.Label(self._member_dist_load_tab, text="wy 1:", font=('Helvetica', 12))
        label.place(x=300, y=397)
        self._boxes[7][5].place(x=350, y=400)
        label = tk.Label(self._member_dist_load_tab, text="wy 2:", font=('Helvetica', 12))
        label.place(x=100, y=447)
        self._boxes[7][6].place(x=150, y=450)
        label = tk.Label(self._member_dist_load_tab, text="wz 1:", font=('Helvetica', 12))
        label.place(x=300, y=447)
        self._boxes[7][7].place(x=350, y=450)
        label = tk.Label(self._member_dist_load_tab, text="wz 2:", font=('Helvetica', 12))
        label.place(x=100, y=497)
        self._boxes[7][8].place(x=150, y=500)
        label = tk.Label(self._member_dist_load_tab, text="Cases:", font=('Helvetica', 12))
        label.place(x=290, y=497)
        self._boxes[7][9].place(x=350, y=500)
        self._buttons.append([])
        self._buttons[7].append(tk.Button(self._member_dist_load_tab, text="Edit", command=self._edit_values))
        self._buttons[7].append(tk.Button(self._member_dist_load_tab, text="Add", command=self._add_values))
        self._buttons[7][0].place(x=550, y=300)
        self._buttons[7][1].place(x=550, y=350)

    def _create_table(self, tab, headings):
        """
        Creates an input table.

        :param tab: Root of tab to put the table on.
        :param headings: Headings to put on the table.
        :return: table
        """

        table = ttk.Treeview(tab, columns=headings, show='headings', height=12)

        # adding scroll bar to table
        scrollbar = ttk.Scrollbar(tab, orient="vertical", command=table.yview)
        scrollbar.pack(side="right", fill="y")
        table.configure(yscrollcommand=scrollbar.set)

        # adding headings to table
        for i in range(len(headings)):
            table.heading(headings[i], text=headings[i])
            table.column(headings[i], width=50)

        # Updates the input boxes to be the current values in the selected table row
        table.bind("<<TreeviewSelect>>", self._get_selected_row)

        table.pack()

        return table

    # noinspection PyUnusedLocal
    def _get_selected_row(self, event):
        """
        Updates the input boxes to be the current values in the selected table row.
        """

        num_col = [3,5,8,7,13,8,9,10]
        t_id = self._table_tabs.index("current")

        row_id = self.Tables[t_id].focus()
        row_info = self.Tables[t_id].item(row_id).get('values')
        for i in range(num_col[t_id]):
            self._boxes[t_id][i].delete(0, 'end')
            self._boxes[t_id][i].insert(0, row_info[i + 1])

    def _edit_values(self):
        """
        Changes the values in the selected row of the table to be the values in the input boxes.
        """

        num_col = [3, 5, 8, 7, 13, 8, 9, 10]
        t_id = self._table_tabs.index("current")

        row_id = self.Tables[t_id].focus()
        row_info = self.Tables[t_id].item(row_id).get('values')
        new_val = [row_info[0]]
        for i in range(num_col[t_id]): new_val.append(self._boxes[t_id][i].get()) # getting values from text boxes

        # Converting strings to floats and bools
        data_new = []
        for i in range(num_col[t_id]):
            if new_val[i + 1] == "True":
                data_new.append(True)
            elif new_val[i + 1] == "False":
                data_new.append(False)
            else:
                data_new.append(float(new_val[i + 1]))

        # Updating Data
        if t_id == 0:   self.Data.EditNode(self, data_new, new_val[0], True, True, row_id)
        elif t_id == 1: self.Data.EditMaterials(self, data_new, new_val[0], True, True, row_id)
        elif t_id == 2: self.Data.EditMembers(self, data_new, new_val[0], True, True, row_id)
        elif t_id == 3: self.Data.EditSupports(self, data_new, new_val[0], True, True, row_id)
        elif t_id == 4: self.Data.EditReleases(self, data_new, new_val[0], True, True, row_id)
        elif t_id == 5: self.Data.EditNodeLoads(self, data_new, new_val[0], True, True, row_id)
        elif t_id == 6: self.Data.EditMemberPointLoad(self, data_new, new_val[0], True, True, row_id)
        elif t_id == 7: self.Data.EditMemberDistLoad(self, data_new, new_val[0], True, True, row_id)

    def _add_values(self):
        """
        Adds a new row to the table with the values in the input boxes.
        """

        num_col = [3, 5, 8, 7, 13, 8, 9, 10]
        t_id = self._table_tabs.index("current")
        num_row = len(self.Tables[t_id].get_children())
        new_val = [num_row]
        for i in range(num_col[t_id]): new_val.append(self._boxes[t_id][i].get()) # getting values from text boxes

        # Converting strings to floats and bools
        data_new = []
        for i in range(num_col[t_id]):
            if new_val[i + 1] == "True":
                data_new.append(True)
            elif new_val[i + 1] == "False":
                data_new.append(False)
            else:
                data_new.append(float(new_val[i + 1]))

        # Updating Data
        if t_id == 0: self.Data.AddNodes(self, data_new, True, True)
        elif t_id == 1: self.Data.AddMaterials(self, data_new, True, True)
        elif t_id == 2: self.Data.AddMembers(self, data_new, True, True)
        elif t_id == 3: self.Data.AddSupports(self, data_new, True, True)
        elif t_id == 4: self.Data.AddReleases(self, data_new, True, True)
        elif t_id == 5: self.Data.AddNodeLoads(self, data_new, True, True)
        elif t_id == 6: self.Data.AddMemberPointLoads(self, data_new, True, True)
        elif t_id == 7: self.Data.AddMemberDistLoads(self, data_new, True, True)

    """ ----------------------------------------------------------------------------------------------"""
    """ ---------------------------------------- 3D Rendering ----------------------------------------"""
    """ ----------------------------------------------------------------------------------------------"""
    # TODO add display system for turning off and on parts of the model (Nodes, Members, Supports, Releases, loads (by case) x 6)
        # todo pop-up window with check boxes to get way to display
        # todo converting all data into text, nodes, lines, and surfaces
        # todo updating print lists
    # TODO add display system for turning off and on results (deflections x 6, reactions x 6, internal loads x 6)
        # todo pop-up window with check boxes to get way to display
        # todo solving for list of deflections along each member and getting key points
        # todo solving for list of internal loads along each member and getting key points
        # todo converting all data into text, nodes, lines, and surfaces
        # todo updating print lists

    # def AddPrintSurface(self, Surface):
    #     """
    #     Adds a surface to be printed in the 3D rendering.
    #
    #     param Surface: Surface's corner coordinates. [X,Y,Z].
    #     """
    #
    #     tri = TDE.triangalize_surface(self, Surface, False)
    #     for i in range(len(tri)): self.PrintSurfaceTri = np.vstack((self.PrintSurfaceTri, tri[i]))
    #     self.UpdateCanves()
    #
    # def AddPrintSolid(self, Solid, FlipNormal: list | None = None):
    #     """
    #     Adds a solid to be printed in the 3D rendering.
    #
    #     param Solid: Solid's corner coordinates. [X,Y,Z].
    #     :param FlipNormal: List of weather each of the solid's face's normals are to be flipped or not.
    #                        Is not needed if the solid is already transformed into triangles.
    #     """
    #
    #     num_surfaces = len(Solid)
    #     for i in range(num_surfaces):
    #         if isinstance(FlipNormal, list): # needs transformed into triangles
    #             tri = TDE.triangalize_surface(self, Solid[i], FlipNormal[i])
    #             for j in range(len(tri)): self.PrintSolidTri = np.vstack((self.PrintSolidTri, [tri[j]]))
    #         else: # already transformed into triangles
    #             self.PrintSolidTri = np.vstack((self.PrintSolidTri, [Solid[i]]))

    def _scroll_down(self, event):
        """
        Ran when the scroll wheel is pressed down.
        Checks if the shift button is held down and if it is start a rotation movement,
        and if not start a pan movement.

        :param event: Event of the scroll wheel being down.
        """

        if bool(event.state & 0x0001):
            self._shift_state_last = True
            self._rot_start(event)
        else:
            self._shift_state_last = False
            self._pan_start(event)

    def _scroll_update(self, event):
        """
        Ran if the cursor moves while the scroll wheel is held down.
        Updates the rotation or panning movement if the shift button is still in the same state, if not end the
        current movement and start the other one.

        :param event: Event of the scroll wheel being down.
        """

        if bool(event.state & 0x0001):
            if not self._shift_state_last:
                self._shift_state_last = True
                self._pan_end(event)
                self._rot_start(event)
            else: self._rot_update(event)
        else:
            if self._shift_state_last:
                self._shift_state_last = False
                self._rot_end(event)
                self._pan_start(event)
            else: self._pan_update(event)

    def _scroll_up(self, event):
        """
        Ran when the scroll wheel is released.
        Ends Rotation or panning movement.

        :param event: Event of the scroll wheel being released.
        """

        if bool(event.state & 0x0001):
            if not self._shift_state_last: self._pan_end(event)
            else: self._rot_end(event)
        else:
            if self._shift_state_last: self._rot_end(event)
            else: self._pan_end(event)

    def _zoom(self, event):
        """
        Ran when the scroll wheel is turned.
        Runs the zooming movement.

        :param event: Event of the scroll wheel being turned.
        """

        if event.delta > 0:
            self.Camera = self.Camera + self.LookDir * self.move_speed * 10
        else:
            self.Camera = self.Camera - self.LookDir * self.move_speed * 10
        self.UpdateCanves()

    def _pan_start(self, event):
        """
        Saves the coordinates of the cursor when a panning movement is started.

        :param event: Event of the scroll wheel being down.
        """

        self._scroll_cords_last = [event.x, event.y]

    def _pan_update(self, event):
        """
        Updates the coordinates of the cursor when a pan movement is in progress and updates the rendering.

        :param event: Event of the scroll wheel being down.
        """

        cords = [event.x, event.y]
        self._pan_move(cords)

    def _pan_end(self, event):
        """
        Updates the coordinates of the cursor when a pan movement is ended and updates the rendering.

        :param event: Event of the scroll wheel being released.
        """

        cords = [event.x, event.y]
        self._pan_move(cords)

    def _pan_move(self, cords):
        """
        Updates the render when a panning movement is taking place.

        :param cords: Current coordinates of the cursor.
        """

        change_cords = [cords[0] - self._scroll_cords_last[0], cords[1] - self._scroll_cords_last[1]]
        self.Camera = self.Camera - self.Up * change_cords[1] * self.move_speed
        self.Camera = self.Camera + np.cross(self.LookDir, self.Up) * change_cords[0] * self.move_speed
        self.UpdateCanves()
        self._scroll_cords_last = cords

    def _rot_start(self, event):
        """
        Saves the coordinates of the cursor when a rotation movement is started.

        :param event: Event of the scroll wheel being down.
        """

        self._scroll_cords_last = [event.x, event.y]

    def _rot_update(self, event):
        """
        Updates the coordinates of the cursor when a rotation movement is in progress and updates the rendering.

        :param event: Event of the scroll wheel being down.
        """

        cords = [event.x, event.y]
        self._rot_move(cords)

    def _rot_end(self, event):
        """
        Updates the coordinates of the cursor when a rotation movement is ended and updates the rendering.

        :param event: Event of the scroll wheel being released.
        """

        cords = [event.x, event.y]
        self._rot_move(cords)

    def _rot_move(self, cords):
        """
        Updates the render when a rotation movement is taking place.

        :param cords: Current coordinates of the cursor.
        """

        change_cords = [cords[0] - self._scroll_cords_last[0], cords[1] - self._scroll_cords_last[1]]
        self.CamYRot = self.CamYRot + change_cords[0] * 0.005
        self.CamXRot = self.CamXRot + change_cords[1] * 0.005
        self.UpdateCanves()
        self._scroll_cords_last = cords

    def UpdateCanves(self):
        """
        Updates the 3D rendering.
        """

        # creates a copy of all the geometry to be shown
        # node = copy(self.PrintNodes)
        # line = copy(self.PrintLines)
        # surf_tri = copy(self.PrintSurfaceTri)
        # solid_tri = copy(self.PrintSolidTri)
        node = copy(self.DisplayData.PrintNodes)
        line = copy(self.DisplayData.PrintLines)
        surf_tri = copy(self.DisplayData.PrintSurfaceTri)
        solid_tri = copy(self.DisplayData.PrintSolidTri)
        text_nodes = copy(self.DisplayData.PrintText[1])
        text_strings = copy(self.DisplayData.PrintText[0])

        # updating all the data about the location and direction of the camera.
        self._update_cam_data()

        # gets all the normals of the tris used for solids being displayed
        solid_tri_normals = TDE.get_normals(solid_tri)

        # removes solid tris that are facing away from the camera.
        solid_tri, solid_tri_normals = TDE.remove_tri_facing_away(self, solid_tri, solid_tri_normals)

        # handles shading and illumination of solids
        tri_color = TDE.illumination(self, solid_tri_normals, len(surf_tri))

        # combing the surface tri and the solid tri into one list
        for i in range(len(surf_tri)): solid_tri = np.append(solid_tri, [surf_tri[i]], axis=0)
        tri = solid_tri

        # transforming all the coordinates to the local coordinate system based on the camera.
        node, line, tri, text_nodes = TDE.transform_to_local(self, node, line, tri, text_nodes)

        # trims all nodes, lines, and tris that are eather to close to the camera or are behind it.
        node, line, tri, tri_color, text_nodes, text_strings = (
            TDE.clip_close(node, line, tri, tri_color, text_nodes, text_strings))

        self.graph.delete("all") # clearing canvas of old rendering

        if len(node) > 0 or len(line) > 0 or len(tri) > 0: # skip if nothing is in the display space.

            # projecting the nodes, lines, and tris and scales to window size.
            w = self.graph.winfo_width()
            h = self.graph.winfo_height()
            node, line, tri, text_nodes = TDE.project(self, node, line, tri, text_nodes)

            # trimming all nodes, lines, and tris that are outside the frame of the camera.
            if len(tri) > 0: tri = tri[:, :, :-1]
            if len(line) > 0: line = line[:, :, :-1]
            node, line, tri, tri_color, text_nodes, text_strings = TDE.clip_edges(node, line, tri, tri_color, text_nodes, text_strings, h, w)

            # printing tris, lines, and nodes
            if len(tri) > 0: self._print_tri(np.array(tri), np.array(tri_color))
            if len(line) > 0: self._print_line(np.array(line))
            if len(node) > 0: self._print_node(np.array(node))
            if len(text_nodes) > 0: self._print_text(np.array(text_strings), np.array(text_nodes))

        self.update()

    def _update_cam_data(self):
        """
        Updating all the data about the location and direction of the camera.
        """

        mat_camera_rot_y = TDE.get_rotation_y_matrix(self.CamYRot)
        mat_camera_rot_x = TDE.get_rotation_x_matrix(self.CamXRot)
        self.LookDir = (np.append(np.array([0, 0, 1]), 1) @ mat_camera_rot_y @ mat_camera_rot_x)[:-1]
        self.LookDir = self.LookDir / np.linalg.norm(self.LookDir)
        self.Target = self.Camera + self.LookDir
        self.Up = (np.append(np.array([0, 1, 0]), 1) @ mat_camera_rot_y @ mat_camera_rot_x)[:-1]
        self.Up = self.Up / np.linalg.norm(self.Up)

    def _print_tri(self, tri, tri_color):
        """
        Prints tris to the canvas.

        :param tri: Tris to be printed.
        :param tri_color: Single int value representing the shading of the tri. Only one value needed as gray scale
                         is used.
        """

        avg_distance = np.sum(tri[:, :, 2], axis=-1)
        idx = np.argsort(avg_distance)
        tri = tri[idx[::-1]]
        tri_color = tri_color[idx[::-1]]

        tri_print = []
        for i in range(len(tri)):
            tri_print.append([tri[i][0][0], tri[i][0][1],
                             tri[i][1][0], tri[i][1][1],
                             tri[i][2][0], tri[i][2][1]])
        for i in range(len(tri_print)):
            color = int(tri_color[i])
            self.graph.create_polygon(tri_print[i], outline='', fill="#%02x%02x%02x" % (color, color, color))

    def _print_line(self, line):
        """
        Prints lines to the canvas.

        :param line: Lines to be printed.
        """

        avg_distance = np.sum(line[:, :, 2], axis=-1)
        idx = np.argsort(avg_distance)
        line = line[idx[::-1]]

        line_print = []
        for i in range(len(line)):
            line_print.append([line[i][0][0], line[i][0][1],
                             line[i][1][0], line[i][1][1]])
        for i in range(len(line_print)): self.graph.create_line(line_print[i], width = 5, fill="red")

    def _print_node(self, node):
        """
        Prints nodes to the canvas.

        :param node: Nodes to be printed.
        """

        idx = np.argsort(node[:,2])
        node = node[idx[::-1]]

        for i in range(len(node)):
            node_p = [node[i][0], node[i][1]]
            self.graph.create_oval(node_p[0] - 4, node_p[1] - 4, node_p[0] + 4, node_p[1] + 4, fill="green")

    def _print_text(self, text_strings, text_nodes):
        """
        Prints text to the canvas.

        :param text_strings: Strings that are to be printed.
        :param text_nodes: Two nodes that the strings are to be printed on top of.
        """
        for i in range(len(text_strings)):
            dx = text_nodes[i*2+1][0]-text_nodes[i*2][0]
            dy = text_nodes[i*2+1][1]-text_nodes[i*2][1]
            angle = np.atan2(dy,dx)

            node = [text_nodes[i*2][0] + dx/2, text_nodes[i*2][1] + dy/2]
            self.graph.create_text(node[0],node[1], text=text_strings[i], angle=angle, font=("Arial", 16), fill="black")

def _select_file_gui(file_types):
    """
    Opens a file explorer dialogue using Tkinter and returns the selected file path.

    :param file_types: List of file types that can be selected.
    :return: File path or None.
    """

    file_path = filedialog.askopenfilename(title="Select a file", filetypes= file_types)
    if file_path: return file_path
    else: return None

def _add_data_to_frame(window, d):
    """
    Creates a 3D frame using the provided Data.

    :param window: Main window.
    :param d: Data to be used to create the 3D frame.
    :return: 3D frame object.
    """

    frame = Frame3D(window.Root)

    for i in range(len(d.Nodes[0])):
        frame.AddNode(d.Nodes[0][i], d.Nodes[1][i], d.Nodes[2][i])

    for i in range(len(d.Materials[0])):
        frame.AddMaterial(d.Materials[0][i], d.Materials[1][i], d.Materials[2][i], d.Materials[3][i], d.Materials[4][i])

    for i in range(len(d.Members[0])):
        frame.AddMember(d.Members[0][i], d.Members[1][i], d.Members[2][i], d.Members[3][i], d.Members[4][i],
                        d.Members[5][i], d.Members[6][i], d.Members[7][i])

    for i in range(len(d.Supports[0])):
        frame.AddSupport(d.Supports[0][i], bool(d.Supports[1][i]), bool(d.Supports[2][i]), bool(d.Supports[3][i]),
                         bool(d.Supports[4][i]),bool(d.Supports[5][i]), bool(d.Supports[6][i]))

    for i in range(len(d.Releases[0])):
        frame.AddReleases(int(d.Releases[0][i]), d.Releases[1][i], d.Releases[2][i], d.Releases[3][i], d.Releases[4][i],
                          d.Releases[5][i], d.Releases[6][i], d.Releases[7][i], d.Releases[8][i], d.Releases[9][i],
                          d.Releases[10][i], d.Releases[11][i], d.Releases[12][i])

    for i in range(len(d.NodeLoad[0])):
        frame.AddNodeLoad(d.NodeLoad[0][i], d.NodeLoad[1][i], d.NodeLoad[2][i], d.NodeLoad[3][i], d.NodeLoad[4][i],
                          d.NodeLoad[5][i], d.NodeLoad[6][i], d.NodeLoad[7][i])

    for i in range(len(d.MemberPointLoad[0])):
        frame.AddMemberPointLoad(int(d.MemberPointLoad[0][i]), d.MemberPointLoad[1][i], d.MemberPointLoad[2][i],
                                 d.MemberPointLoad[3][i], d.MemberPointLoad[4][i], d.MemberPointLoad[5][i],
                                 d.MemberPointLoad[6][i], d.MemberPointLoad[7][i], int(d.MemberPointLoad[8][i]))

    for i in range(len(d.MemberDistLoad[0])):
        frame.AddMemberDistLoad(int(d.MemberDistLoad[0][i]), d.MemberDistLoad[1][i], d.MemberDistLoad[2][i],
                                d.MemberDistLoad[3][i], d.MemberDistLoad[4][i], d.MemberDistLoad[5][i],
                                d.MemberDistLoad[6][i], d.MemberDistLoad[7][i], d.MemberDistLoad[8][i],
                                int(d.MemberDistLoad[9][i]))
    return frame
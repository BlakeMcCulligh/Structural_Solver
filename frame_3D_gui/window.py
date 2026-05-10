"""
Holds the class that handes everything to do with 3D frames gui
including thair windows, inputing, exporting, and saving.
"""

import tkinter as tk
from tkinter import ttk,  filedialog
from copy import copy
import pandas as pd
import numpy as np

from frame_3D_gui import export
from frame_3D_gui.data import Data
from frame_3D_gui.opening import open_frame, open_results
from frame_3D_gui.optimize_pop_up import OptimizationPopUp
from frame_3D_gui.results import Results
from frame_3D_gui.save import save_frame, save_results
from frame_3D_solver.main import Frame3D
# noinspection PyPep8Naming
import drawing_3D.engine_3D as TDE
import frame_3D_solver.helper_functions as hf

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
    Main window Object for 3D Frame Analysises and anything elce to do with 3D Frames.
    """

    def __init__(self, Root):
        """
        Constructor for the main window of the 3D frame.

        :param Root: Root of the window.
        """

        tk.Frame.__init__(self, Root)

        self.FilePath = None # File path for saving.

        self.Root = Root

        self._center_window()

        self._create_top_menu()

        self.Data = Data()

        # setting anlysis objects to none.
        self.Frame = None
        self.Results = None
        self.OptimizationResults = None

        # 3D rendering canvas
        self.graph = tk.Canvas(Root, bg="white")
        self.graph.pack(fill="both", expand=True)

        # 3D camera Data
        self.Camera = np.array([0.0, 0.0, 0.0])
        self.Up = np.array([0.0, 1.0, 0.0])
        self.LookDir = np.array([0.0, 0.0, 1.0])
        self.CamYRot = 0
        self.CamXRot = 0
        self.Target = np.array([0.0, 0.0, 1.0])
        self.FOV = 90
        self.Z_FAR, self.Z_NEAR = 1000, 0.1

        self.LIGHT_DIR = np.array([0, 0, -1]) # 3D rendering ligingting direction unit vector

        # arrays of geomitry being displayed in the 3D rendering
        self.PrintNodes: np.ndarray = np.empty((0, 3))
        self.PrintLines: np.ndarray = np.empty((0, 2, 3))
        self.PrintSurfaceTri: np.ndarray = np.empty((0, 3, 3))
        self.PrintSolidTri: np.ndarray = np.empty((0, 3, 3))

        # binding 3D rendering movement inputs
        self.Root.bind("<MouseWheel>", self._zoom)
        self._scrool_cords_last = []
        self._shift_state_last = False
        self.Root.bind("<Button-2>", self._scrool_down)
        self.Root.bind("<B2-Motion>", self._scrool_update)
        self.Root.bind("<ButtonRelease-2>", self._scrool_up)

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

        analysis_menu = tk.Menu(menubar, tearoff=False)
        analysis_menu.add_command(label='Linear Analysis', command=self._linear_analysis)
        analysis_menu.add_command(label='Global Optimization', command=self._optimization_window)
        menubar.add_cascade(label="Analysis", menu=analysis_menu)

    def _exit(self):
        """
        CLoses program
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

    def _open(self):
        """
        Opens files under specified name and location.
        """

        open_frame(self)
        open_results(self, self.FilePath)

    def _import_nodes(self):
        """
        Imports coordanits for nodes from the specified Excel file.
        Headers must be "X, Y, Z" and there must only be one sheet.
        """

        file_types = [("Excel Files", "*.xlsx")]
        file_path = _select_file_gui(file_types) # geting file path from user
        nodes_df = pd.read_excel(file_path)
        nodes = [nodes_df["X"].tolist(), nodes_df["Y"].tolist(), nodes_df["Z"].tolist()]
        self.Data.AddNodes(self, nodes, True, True)

    def _import_materials(self):
        """
        Imports material properties from the specified Excel file.
        Headers must be "E, G, nu, rho, fy" and there must only be one sheet.
        """

        file_types = [("Excel Files", "*.xlsx")]
        file_path = _select_file_gui(file_types) # geting file path from user
        m_df = pd.read_excel(file_path)
        materials = [m_df["E"].tolist(),m_df["G"].tolist(),m_df["nu"].tolist(),m_df["rho"].tolist(),m_df["fy"].tolist()]
        self.Data.AddMaterials(self, materials, True, True)

    def _import_members(self):
        """
        Imports member end node indecies, material index, if cross-sections properties are to be set,
        and the cross-section properties for the specified Excel file.
        Headers must be "i Node, j Node, Material, Set Cross-Section Properties, A, Iy, Iz, J"
        and there must only be one sheet.
        """

        file_types = [("Excel Files", "*.xlsx")]
        file_path = _select_file_gui(file_types) # geting file path from user
        m_df = pd.read_excel(file_path)
        members = [m_df["i Node"].tolist(), m_df["j Node"].tolist(), m_df["Material"].tolist(),
                   m_df["Set Cross-Section Properties"].tolist(), m_df["A"].tolist(), m_df["Iy"].tolist(),
                   m_df["Iz"].tolist(), m_df["J"].tolist()]
        self.Data.AddMembers(self, members, True, True)

    def _import_supports(self):
        """
        Imports supports node indeces, and what degrees of freedom are to be supported for the specififed Excel file.
        Headers must be "Node, DX, DY, DZ, RX, RY, RZ" and there must only be one sheet.
        """

        file_types = [("Excel Files", "*.xlsx")]
        file_path = _select_file_gui(file_types) # geting file path from user
        s_df = pd.read_excel(file_path)
        supports = [s_df["Node"].tolist(),s_df["DX"].tolist(),s_df["DY"].tolist(),s_df["DZ"].tolist(),
                    s_df["RX"].tolist(),s_df["RY"].tolist(),s_df["RZ"].tolist()]
        self.Data.AddSupports(self, supports, True, True)

    def _import_releases(self):
        """
        Inporst releced member indeces, and what ends and degrees of freedom are to be releced for the specified
        Excel file.
        Headers must be "Member, i DX, i DY, i DZ, i RX, i RY, i RZ, j DX, j DY, j DZ, j RX, j RY, j RZ" and there
        must only be one sheet.
        """

        file_types = [("Excel Files", "*.xlsx")]
        file_path = _select_file_gui(file_types) # geting file path from user
        r_df = pd.read_excel(file_path)
        releases = [r_df["Member"].tolist(), r_df["i DX"].tolist(), r_df["i DY"].tolist(), r_df["i DZ"].tolist(),
                    r_df["i RX"].tolist(), r_df["i RY"].tolist(), r_df["i RZ"].tolist(),r_df["j DX"].tolist(),
                    r_df["j DY"].tolist(), r_df["j DZ"].tolist(), r_df["j RX"].tolist(), r_df["j RY"].tolist(),
                    r_df["j RZ"].tolist()]
        self.Data.AddReleases(self, releases, True, True)

    def _import_node_loads(self):
        """
        Imports node load, node indeces, what directions and magnitueds the load is to be in,
         and the index of the load case.
        Headers must be "Node, PX, PY, PZ, MX, MY, MZ, Case" and there must only be one sheet.
        """

        file_types = [("Excel Files", "*.xlsx")]
        file_path = _select_file_gui(file_types) # geting file path from user
        n_df = pd.read_excel(file_path)
        node_loads = [n_df["Node"].tolist(),n_df["PX"].tolist(),n_df["PY"].tolist(),n_df["PZ"].tolist(),
                     n_df["MX"].tolist(),n_df["MY"].tolist(),n_df["MZ"].tolist(),n_df["Case"].tolist()]
        self.Data.AddNodeLoads(self, node_loads, True, True)

    def _import_member_point_loads(self):
        """
        Imports member point loads, member indeces, location, what directions and magnitueds the load is to be in
        and the index of the load case.
        Headers must be "Member, X, PX, PY, PZ, MX, MY, MZ, Case" and there must only be one sheet.
        """

        file_types = [("Excel Files", "*.xlsx")]
        file_path = _select_file_gui(file_types) # geting file path from user
        m_df = pd.read_excel(file_path)
        member_point_loads = [m_df["Member"].tolist(),m_df["X"].tolist(),m_df["PX"].tolist(),m_df["PY"].tolist(),
                            m_df["PZ"].tolist(),m_df["MX"].tolist(),m_df["MY"].tolist(),m_df["MZ"].tolist(),
                            m_df["Case"].tolist()]
        self.Data.AddMemberPointLoads(self, member_point_loads, True, True)

    def _import_member_dist_loads(self):
        """
        Imports member distributed loads: member indeces, start and end locations, start and end force magnitued and
        directions, and the index of the load case
        Headers must be "Member, X1, X2, WX1, WX2, WY1, WY2, WZ1, WZ2, Case" and there must only be one sheet.
        """

        file_types = [("Excel Files", "*.xlsx")]
        file_path = _select_file_gui(file_types) # geting file path from user
        m_df = pd.read_excel(file_path)
        member_dist_loads = [m_df["Member"].tolist(),m_df["X1"].tolist(),m_df["X2"].tolist(),m_df["WX1"].tolist(),
                           m_df["WX2"].tolist(),m_df["WY1"].tolist(),m_df["WY2"].tolist(),m_df["WZ1"].tolist(),
                           m_df["WZ2"].tolist(),m_df["Case"].tolist()]
        self.Data.AddMemberDistLoads(self, member_dist_loads, True, True)

    def _linear_analysis(self):
        """
         Runs linear analysis on the defined frame, and saves the results.
        """

        frame = _add_data_to_frame(self.Data)

        frame.PreAnalysisLinear()
        D, DX, DY, DZ, RX, RY, RZ, weight, reactions, internal_forces = frame.AnalysisLinear(get_weight= True,
                                                                                             get_reactions= True,
                                                                                             get_internal_forces= True)
        self.Frame = frame
        self.Results = Results()
        self.Results.AddNodalDeflections(DX, DY, DZ, RX, RY, RZ)
        self.Results.AddWeight(weight)
        self.Results.AddReactions(reactions)
        self.Results.AddInternalForces(internal_forces)

        self._save()
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

        :param GroupAssignments: Member groups assignments for members withou specified cross-sections.
        :param GroupTypes: What kind of cross-section each member group is.
        :param LowerBound: Lower bounds on optimization variables (dimentions on the cross-sections).
        :param UpperBound: Upper bounds on optimization variables (dimentions on the cross-sections).
        :param CostFunction: String that carrys the equation to get the cost of the frame.
                             The optimizer minimizes the value this equation returns.
        :param WeightRun: If the weight of the members are needed for the cost equation.
        :param ReactionRun: If the reactions of the frame are needed for the cost equation.
        :param InternalForcesRun: If the internal forces of the frame are needed for the cost equation.
        """

        save_frame(self.Data)

        frame = _add_data_to_frame(self.Data)

        frame.PreAnalysisLinear()
        results = frame.optimize(GroupAssignments, GroupTypes, LowerBound, UpperBound, CostFunction, WeightRun,
                                 ReactionRun, InternalForcesRun)
        self.Frame = frame
        self.OptimizationResults = results

        self.FilePath = None
        member_indices =  self._update_frame_to_optimization_results(GroupAssignments, GroupTypes)
        self._linear_analysis() # runs an analysis of the optimized frame to get the full analysis Results
        self._save_optimization_results(member_indices, GroupAssignments, GroupTypes)

    def _update_frame_to_optimization_results(self, group_assignments, group_types):
        """
        Updates the frame's cross-sections to the results found the optimizer.

        :param group_assignments: Member groups assignments for members withou specified cross-sections.
        :param group_types: What kind of cross-section each member group is.
        :return: member_indices: Indeces of the optimized member cross-sections.
        """

        x = self.OptimizationResults.x
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

        return member_indices

    def _save_optimization_results(self, member_indices, group_assignments, group_types):
        """
        Exports the optimization results to an Excel file.

        :param member_indices: Indeces of the optimized member cross-sections.
        :param group_assignments: Member groups assignments for members withou specified cross-sections.
        :param group_types: What kind of cross-section each member group is.
        """

        x = self.OptimizationResults.x
        cost = [self.OptimizationResults.fun]
        opt_cross_section_props = hf.get_cross_section_props(x, group_assignments, group_types)

        results = []
        curent_x_index = 0
        for i in range(len(opt_cross_section_props)):
            memeber_index = [member_indices[i]]
            group_type = [group_types[i]]

            dim = []
            if group_type == ["Angle"]:
                dim.append(x[curent_x_index])
                dim.append(x[curent_x_index + 1])
                dim.append(x[curent_x_index + 2])
                curent_x_index += 3
            elif group_type == ["RectHSS"]:
                dim.append(x[curent_x_index])
                dim.append(x[curent_x_index + 1])
                dim.append(x[curent_x_index + 2])
                curent_x_index += 3
            elif group_type == ["SquareHSS"]:
                dim.append(x[curent_x_index])
                dim.append("N/A")
                dim.append(x[curent_x_index + 1])
                curent_x_index += 2
            elif group_type == ["TubeHSS"]:
                dim.append(x[curent_x_index])
                dim.append("N/A")
                dim.append(x[curent_x_index + 1])
                curent_x_index += 2

            cs = opt_cross_section_props[i]

            results.append(memeber_index + group_type + dim + cs)

        export.export_optimization_results(results, cost)

    """ ----------------------------------------------------------------------------------------------"""
    """ ---------------------------------------- INPUT TABLES ----------------------------------------"""
    """ ----------------------------------------------------------------------------------------------"""

    def _create_table_window(self):
        """
        Creates window with all the input tables.
        """

        # TODO add lables to input boxes
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
        self._boxes[0].append(tk.Entry(self._node_tab))
        self._boxes[0].append(tk.Entry(self._node_tab))
        self._boxes[0].append(tk.Entry(self._node_tab))
        self._boxes[0][0].place(x=350, y=300) # X
        self._boxes[0][1].place(x=350, y=350) # Y
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
        self._boxes[1].append(tk.Entry(self._mat_tab))
        self._boxes[1].append(tk.Entry(self._mat_tab))
        self._boxes[1].append(tk.Entry(self._mat_tab))
        self._boxes[1].append(tk.Entry(self._mat_tab))
        self._boxes[1].append(tk.Entry(self._mat_tab))
        self._boxes[1][0].place(x=150, y=300)
        self._boxes[1][1].place(x=350, y=300)
        self._boxes[1][2].place(x=150, y=350)
        self._boxes[1][3].place(x=350, y=350)
        self._boxes[1][4].place(x=150, y=400)
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
        for i in range(8): self._boxes[2].append(tk.Entry(self._member_tab))
        self._boxes[2][0].place(x=150, y=300)
        self._boxes[2][1].place(x=350, y=300)
        self._boxes[2][2].place(x=150, y=350)
        self._boxes[2][3].place(x=350, y=350)
        self._boxes[2][4].place(x=150, y=400)
        self._boxes[2][5].place(x=350, y=400)
        self._boxes[2][6].place(x=150, y=450)
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
        for i in range(7): self._boxes[3].append(tk.Entry(self._support_tab))
        self._boxes[3][0].place(x=150, y=300)
        self._boxes[3][1].place(x=350, y=300)
        self._boxes[3][2].place(x=150, y=350)
        self._boxes[3][3].place(x=350, y=350)
        self._boxes[3][4].place(x=150, y=400)
        self._boxes[3][5].place(x=350, y=400)
        self._boxes[3][6].place(x=150, y=450)
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
        for i in range(13): self._boxes[4].append(tk.Entry(self._release_tab))
        self._boxes[4][0].place(x=150, y=300)
        self._boxes[4][1].place(x=350, y=300)
        self._boxes[4][2].place(x=150, y=350)
        self._boxes[4][3].place(x=350, y=350)
        self._boxes[4][4].place(x=150, y=400)
        self._boxes[4][5].place(x=350, y=400)
        self._boxes[4][6].place(x=150, y=450)
        self._boxes[4][7].place(x=350, y=450)
        self._boxes[4][8].place(x=150, y=500)
        self._boxes[4][9].place(x=350, y=500)
        self._boxes[4][10].place(x=150, y=550)
        self._boxes[4][11].place(x=350, y=550)
        self._boxes[4][12].place(x=150, y=600)
        self._buttons.append([])
        self._buttons[4].append(tk.Button(self._release_tab, text="Edit", command=self._edit_values))
        self._buttons[4].append(tk.Button(self._release_tab, text="Add", command=self._add_values))
        self._buttons[4][0].place(x=550, y=300)
        self._buttons[4][1].place(x=550, y=350)

        # Nodal Loads
        self._node_load_tab = ttk.Frame(self._table_tabs)
        self._table_tabs.add(self._node_load_tab, text="Node Loads")
        self.Tables.append(self._create_table(self._node_load_tab, ("Index", "i Node", "P X", "P Y", "P Z.", "M X",
                                                              "M Y", "M Z","Casee")))
        self._boxes.append([])
        for i in range(8): self._boxes[5].append(tk.Entry(self._node_load_tab))
        self._boxes[5][0].place(x=150, y=300)
        self._boxes[5][1].place(x=350, y=300)
        self._boxes[5][2].place(x=150, y=350)
        self._boxes[5][3].place(x=350, y=350)
        self._boxes[5][4].place(x=150, y=400)
        self._boxes[5][5].place(x=350, y=400)
        self._boxes[5][6].place(x=150, y=450)
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
                                                                      "P Z.", "M X", "M Y", "M Z", "Casee")))
        self._boxes.append([])
        for i in range(9): self._boxes[6].append(tk.Entry(self._member_point_load_tab))
        self._boxes[6][0].place(x=150, y=300)
        self._boxes[6][1].place(x=350, y=300)
        self._boxes[6][2].place(x=150, y=350)
        self._boxes[6][3].place(x=350, y=350)
        self._boxes[6][4].place(x=150, y=400)
        self._boxes[6][5].place(x=350, y=400)
        self._boxes[6][6].place(x=150, y=450)
        self._boxes[6][7].place(x=350, y=450)
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
                                                                    "wx 2", "wy 1", "wy 2", "wz 1", "wz 2", "Casee")))
        self._boxes.append([])
        for i in range(10): self._boxes[7].append(tk.Entry(self._member_dist_load_tab))
        self._boxes[7][0].place(x=150, y=300)
        self._boxes[7][1].place(x=350, y=300)
        self._boxes[7][2].place(x=150, y=350)
        self._boxes[7][3].place(x=350, y=350)
        self._boxes[7][4].place(x=150, y=400)
        self._boxes[7][5].place(x=350, y=400)
        self._boxes[7][6].place(x=150, y=450)
        self._boxes[7][7].place(x=350, y=450)
        self._boxes[7][8].place(x=150, y=500)
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

        table = ttk.Treeview(tab, columns=headings, show='headings', height=10)

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
        for i in range(num_col[t_id]): new_val.append(self._boxes[t_id][i].get()) # geting values from text boxes

        # Converting strings to floats and bools
        data_new = []
        for i in range(num_col[t_id]):
            if new_val[i + 1] == "True" or new_val[i + 1] == "False":
                data_new.append(bool(new_val[i + 1]))
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
        for i in range(num_col[t_id]): new_val.append(self._boxes[t_id][i].get()) # geting values from text boxes

        # Converting strings to floats and bools
        data_new = []
        for i in range(num_col[t_id]):
            if new_val[i + 1] == "True" or new_val[i + 1] == "False":
                data_new.append(bool(new_val[i + 1]))
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

    def AddPrintNode(self, Node):
        """
        Adds a node to be printed in the 3D rendering.

        :param Node: Node cordinates. [X,Y,Z].
        """

        self.PrintNodes = np.vstack((self.PrintNodes, Node))
        self.UpdateCanves()

    def AddPrintLine(self, Line):
        """
        Adds a line to be printed in the 3D rendering.

        :param Line: Line end Node indeces in the print Node array [i Node, j Node].
        """

        l = []
        for i in range(len(Line)):
            l.append(int(Line[i]))
        l = np.array(l)
        self.PrintLines = np.vstack((self.PrintLines, [self.PrintNodes[l]]))
        self.UpdateCanves()

    def AddPrintSurface(self, Surface):
        """
        Adds a surface to be printed in the 3D rendering.

        :param Surface: Serface's corner cordinates. [X,Y,Z].
        """

        tri = TDE.triangalize_surface(self, Surface, False)
        for i in range(len(tri)): self.PrintSurfaceTri = np.vstack((self.PrintSurfaceTri, tri[i]))
        self.UpdateCanves()

    def AddPrintSolid(self, Solid, FlipNormal: list | None = None):
        """
        Adds a solid to be printed in the 3D rendering.

        :param Solid: Solid's corner cordinates. [X,Y,Z].
        :param FlipNormal: List of weather each of the solid's face's normals are to be flipped or not.
                           Is not needed if the solid is already triangalized.
        """

        num_surfaces = len(Solid)
        for i in range(num_surfaces):
            if isinstance(FlipNormal, list): # needs triangalized
                tri = TDE.triangalize_surface(self, Solid[i], FlipNormal[i])
                for j in range(len(tri)): self.PrintSolidTri = np.vstack((self.PrintSolidTri, [tri[j]]))
            else: # already triangalized
                self.PrintSolidTri = np.vstack((self.PrintSolidTri, [Solid[i]]))

    def _scrool_down(self, event):
        """
        Ran when the scroll wheel is pressed down.
        Checks if the shift buttion is held down and if it is start a rotation movement,
        and if not start a pan movement.

        :param event: Event of the scrool wheel being down.
        """

        if bool(event.state & 0x0001):
            self._shift_state_last = True
            self._rot_start(event)
        else:
            self._shift_state_last = False
            self._pan_start(event)

    def _scrool_update(self, event):
        """
        Ran if the currser moves while the scrool wheel is held down.
        Updates the rotation or panning movement if the shift button is still in the same state, if not end the
        current movement and start the other one.

        :param event: Event of the scrool wheel being down.
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

    def _scrool_up(self, event):
        """
        Ran when the scrool wheel is released.
        Ends Roation or paning movement.

        :param event: Event of the scrool wheel being releced.
        """

        if bool(event.state & 0x0001):
            if not self._shift_state_last: self._pan_end(event)
            else: self._rot_end(event)
        else:
            if self._shift_state_last: self._rot_end(event)
            else: self._pan_end(event)

    def _zoom(self, event):
        """
        Ran when the scrool wheel is turned.
        Runs the zooming movement.

        :param event: Event of the scrool wheel being turned.
        """

        if event.delta > 0:
            self.Camera = self.Camera + self.LookDir * 0.1
        else:
            self.Camera = self.Camera - self.LookDir * 0.1
        self.UpdateCanves()

    def _pan_start(self, event):
        """
        Saves the codinats of the currser when a panning movement is started.

        :param event: Event of the scrool wheel being down.
        """

        self._scrool_cords_last = [event.x, event.y]

    def _pan_update(self, event):
        """
        Updates the coordinates of the currser when a pan movement is in progress and updates the rendering.

        :param event: Event of the scrool wheel being down.
        """

        cords = [event.x, event.y]
        self._pan_move(cords)

    def _pan_end(self, event):
        """
        Updates the coordinates of the currser when a pan movement is endded and updates the rendering.

        :param event: Event of the scrool wheel being releced.
        """

        cords = [event.x, event.y]
        self._pan_move(cords)

    def _pan_move(self, cords):
        """
        Updates the render when a panning movement is taking place.

        :param cords: Current cordinates of the currser.
        """

        change_cords = [cords[0] - self._scrool_cords_last[0], cords[1] - self._scrool_cords_last[1]]
        self.Camera = self.Camera - self.Up * change_cords[1] * 0.01
        self.Camera = self.Camera + np.cross(self.LookDir, self.Up) * change_cords[0] * 0.01
        self.UpdateCanves()
        self._scrool_cords_last = cords

    def _rot_start(self, event):
        """
        Saves the codinats of the currser when a rotation movement is started.

        :param event: Event of the scrool wheel being down.
        """

        self._scrool_cords_last = [event.x, event.y]

    def _rot_update(self, event):
        """
        Updates the coordinates of the currser when a rotation movement is in progress and updates the rendering.

        :param event: Event of the scrool wheel being down.
        """

        cords = [event.x, event.y]
        self._rot_move(cords)

    def _rot_end(self, event):
        """
        Updates the coordinates of the currser when a rotation movement is ended and updates the rendering.

        :param event: Event of the scrool wheel being releced.
        """

        cords = [event.x, event.y]
        self._rot_move(cords)

    def _rot_move(self, cords):
        """
        Updates the render when a rotation movement is taking place.

        :param cords: Current cordinates of the currser.
        """

        change_cords = [cords[0] - self._scrool_cords_last[0], cords[1] - self._scrool_cords_last[1]]
        self.CamYRot = self.CamYRot + change_cords[0] * 0.005
        self.CamXRot = self.CamXRot + change_cords[1] * 0.005
        self.UpdateCanves()
        self._scrool_cords_last = cords

    def UpdateCanves(self):
        """
        Updates the 3D rendering.
        """

        # creates a copy of all the geomitry to be shown
        node = copy(self.PrintNodes)
        line = copy(self.PrintLines)
        surf_tri = copy(self.PrintSurfaceTri)
        solid_tri = copy(self.PrintSolidTri)

        # updating all the data about the location and direction of the camera.
        self._update_cam_data()

        # gets all the normals of the tris used for solids being displeyed
        solid_tri_normals = TDE.get_normals(solid_tri)

        # removes solid tris that are facing away from the camera.
        solid_tri, solid_tri_normals = TDE.remove_tri_faceing_away(self, solid_tri, solid_tri_normals)

        # handels shadding and ilumination of solids
        tri_color = TDE.illumination(self, solid_tri_normals, len(surf_tri))

        # combing the surfice tri and the solid tri into one list
        for i in range(len(surf_tri)): solid_tri = np.append(solid_tri, surf_tri[i], axis=0)
        tri = solid_tri

        # transforming all the cordinates to the local cordinate system based on the camera.
        node, line, tri = TDE.transform_to_local(self, node, line, tri)

        # trims all nodes, lines, and tris that are eather to close to the camera or are behinde it.
        node, line, tri, tri_color = TDE.clip_close(node, line, tri, tri_color)

        self.graph.delete("all") # clearing canvise of old rendering

        if len(node) > 0 or len(line) > 0 or len(tri) > 0: # skip if nothing is in the display space.

            # projecting the nodes, lines, and tris and scales to window size.
            w = self.graph.winfo_width()
            h = self.graph.winfo_height()
            node, line, tri = TDE.project(self, node, line, tri)

            # trimming all nodes, lines, and tris that are outside the frame of the camera.
            if len(tri) > 0: tri = tri[:, :, :-1]
            if len(line) > 0: line = line[:, :, :-1]
            node, line, tri, tri_color = TDE.clip_eadges(node, line, tri, tri_color, h, w)

            # printing tris, lines, and nodes
            if len(tri) > 0: self._print_tri(np.array(tri), np.array(tri_color))
            if len(line) > 0: self._print_line(np.array(line))
            if len(node) > 0: self._print_node(np.array(node))

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
        Prints tris to the canves.

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
        Prints lines to the canves.

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
        Prints nodes to the canves.

        :param node: Nodes to be printed.
        """

        idx = np.argsort(node[:,2])
        node = node[idx[::-1]]

        for i in range(len(node)):
            node_p = [node[i][0], node[i][1]]
            self.graph.create_oval(node_p[0] - 4, node_p[1] - 4, node_p[0] + 4, node_p[1] + 4, fill="green")


def _select_file_gui(file_types):
    """
    Opens a file explorer dialog using Tkinter and returns the selected file path.

    :param file_types: List of file types that can be selected.
    :return: File path or None.
    """

    file_path = filedialog.askopenfilename(title="Select a file", filetypes= file_types)
    if file_path: return file_path
    else: return None

def _add_data_to_frame(d):
    """
    Creates a 3D frame using the provided Data.

    :param d: Data to be used to create the 3D frame
    :return: 3D frame object.
    """

    frame = Frame3D()

    for i in range(len(d.Nodes[0])):
        frame.AddNode(d.Nodes[0][i], d.Nodes[1][i], d.Nodes[2][i])

    for i in range(len(d.Materials[0])):
        frame.AddMaterial(d.Materials[0][i], d.Materials[1][i], d.Materials[2][i], d.Materials[3][i], d.Materials[4][i])

    for i in range(len(d.Members[0])):
        frame.AddMember(d.Members[0][i], d.Members[1][i], d.Members[2][i], d.Members[3][i], d.Members[4][i],
                        d.Members[5][i], d.Members[6][i], d.Members[7][i])

    for i in range(len(d.Supports[0])):
        frame.AddSupport(d.Supports[0][i], d.Supports[1][i], d.Supports[2][i], d.Supports[3][i], d.Supports[4][i],
                         d.Supports[5][i], d.Supports[6][i])

    for i in range(len(d.Releases[0])):
        frame.AddReleases(d.Releases[0][i], d.Releases[1][i], d.Releases[2][i], d.Releases[3][i], d.Releases[4][i],
                          d.Releases[5][i], d.Releases[6][i], d.Releases[7][i], d.Releases[8][i], d.Releases[9][i],
                          d.Releases[10][i], d.Releases[11][i], d.Releases[12][i])

    for i in range(len(d.NodeLoad[0])):
        frame.AddNodeLoad(d.NodeLoad[0][i], d.NodeLoad[1][i], d.NodeLoad[2][i], d.NodeLoad[3][i], d.NodeLoad[4][i],
                          d.NodeLoad[5][i], d.NodeLoad[6][i], d.NodeLoad[7][i])

    for i in range(len(d.MemberPointLoad[0])):
        frame.AddMemberPointLoad(d.MemberPointLoad[0][i], d.MemberPointLoad[1][i], d.MemberPointLoad[2][i],
                                 d.MemberPointLoad[3][i], d.MemberPointLoad[4][i], d.MemberPointLoad[5][i],
                                 d.MemberPointLoad[6][i], d.MemberPointLoad[7][i], d.MemberPointLoad[8][i])

    for i in range(len(d.MemberDistLoad[0])):
        frame.addMemberDistLoad(d.MemberDistLoad[0][i], d.MemberDistLoad[1][i], d.MemberDistLoad[2][i],
                                d.MemberDistLoad[3][i], d.MemberDistLoad[4][i], d.MemberDistLoad[5][i],
                                d.MemberDistLoad[6][i], d.MemberDistLoad[7][i], d.MemberDistLoad[8][i],
                                d.MemberDistLoad[9][i])
    return frame
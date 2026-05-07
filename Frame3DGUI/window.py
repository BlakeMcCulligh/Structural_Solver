import tkinter as tk
from tkinter import ttk,  filedialog
from copy import copy
import pandas as pd
import numpy as np

from Frame3DGUI.inputData import data
from Frame3DGUI.opening import openFrame
from Frame3DGUI.optimizePopUp import OptimizationPopUp
from Frame3DGUI.results import Results
from Frame3DGUI.save import saveFrame, saveResults
from StructuralAnalysis.frame3DSolver.__main__ import Frame3D
# noinspection PyPep8Naming
import ThreeDDrawing.ThreeDEngine as TDE

class MainWindow(tk.Frame):
    def __init__(self, root):
        tk.Frame.__init__(self, root)

        self.root = root

        self.create_top_menu()

        self.centerWindow()

        self.data = data()

        self.Frame = None
        self.results = None
        self.Optimization_Results = None

        self.graph = tk.Canvas(root, bg="white")
        self.graph.pack(fill="both", expand=True)
        self.camera = np.array([0.0, 0.0, 0.0])
        self.up = np.array([0.0, 1.0, 0.0])
        self.lookDir = np.array([0.0, 0.0, 1.0])
        self.cam_Y_rot = 0
        self.cam_X_rot = 0
        self.target = np.array([0.0, 0.0, 1.0])
        self.FOV = 90
        self.Z_far, self.Z_near = 1000, 0.1
        self.light_direction = np.array([0, 0, -1])

        self.printNodes: np.ndarray = np.empty((0, 3))
        self.printLines: np.ndarray = np.empty((0, 2, 3))
        self.surfaceTri: np.ndarray = np.empty((0, 3, 3))
        self.solidTri: np.ndarray = np.empty((0, 3, 3))

        self.root.bind("<MouseWheel>", self.zoom)

        self.scroolCordsLast = []
        self.shiftStateLast = False
        self.root.bind("<Button-2>", self.scroolDown)
        self.root.bind("<B2-Motion>", self.scroolUpdate)
        self.root.bind("<ButtonRelease-2>", self.scroolUp)

        self.tableWindow = None
        self.tableTabs = None
        self.tables = []
        self.Boxes = []
        self.Buttons = []

        self.nodeTab = None
        self.matTab = None
        self.memberTab = None
        self.supportTab = None
        self.releaseTab = None
        self.nodeLoadTab = None
        self.memberPointLoadTab = None
        self.memberDistLoadTab = None

        self.createTableWindow()

    def createTableWindow(self):
        self.tableWindow = tk.Toplevel(self.root)
        self.tableWindow.title("Sheets")
        self.tableWindow.geometry("800x700")
        self.tableWindow.resizable(False, True)

        self.tableTabs = ttk.Notebook(self.tableWindow)
        self.tableTabs.pack(expand=True, fill="both")

        self.nodeTab = ttk.Frame(self.tableTabs)

        # Nodes
        self.tableTabs.add(self.nodeTab, text="Nodes")
        self.tables.append(self.createTable(self.nodeTab, ("Index", "X", "Y", "Z")))
        self.Boxes.append([])
        self.Boxes[0].append(tk.Entry(self.nodeTab))
        self.Boxes[0].append(tk.Entry(self.nodeTab))
        self.Boxes[0].append(tk.Entry(self.nodeTab))
        self.Boxes[0][0].place(x=350, y=300) # X
        self.Boxes[0][1].place(x=350, y=350) # Y
        self.Boxes[0][2].place(x=350, y=400) # Z
        self.Buttons.append([])
        self.Buttons[0].append(tk.Button(self.nodeTab, text="Edit", command=self.editValues))
        self.Buttons[0].append(tk.Button(self.nodeTab, text="Add", command=self.addValues))
        self.Buttons[0][0].place(x=550, y=300)
        self.Buttons[0][1].place(x=550, y=350)
        # todo add delete button all tabs

        # Materials
        self.matTab = ttk.Frame(self.tableTabs)
        self.tableTabs.add(self.matTab, text="Materials")
        self.tables.append(self.createTable(self.matTab, ("Index", "E", "G", "nu", "rho", "fy")))
        self.Boxes.append([])
        self.Boxes[1].append(tk.Entry(self.matTab))
        self.Boxes[1].append(tk.Entry(self.matTab))
        self.Boxes[1].append(tk.Entry(self.matTab))
        self.Boxes[1].append(tk.Entry(self.matTab))
        self.Boxes[1].append(tk.Entry(self.matTab))
        self.Boxes[1][0].place(x=150, y=300)
        self.Boxes[1][1].place(x=350, y=300)
        self.Boxes[1][2].place(x=150, y=350)
        self.Boxes[1][3].place(x=350, y=350)
        self.Boxes[1][4].place(x=150, y=400)
        self.Buttons.append([])
        self.Buttons[1].append(tk.Button(self.matTab, text="Edit", command=self.editValues))
        self.Buttons[1].append(tk.Button(self.matTab, text="Add", command=self.addValues))
        self.Buttons[1][0].place(x=550, y=300)
        self.Buttons[1][1].place(x=550, y=350)

        # Members
        self.memberTab = ttk.Frame(self.tableTabs)
        self.tableTabs.add(self.memberTab, text="Members")
        self.tables.append(self.createTable(self.memberTab, ("Index", "i Node", "j Node", "Material Id", "Set C.S.", "A", "Iy", "Iz", "J")))
        self.Boxes.append([])
        for i in range(8): self.Boxes[2].append(tk.Entry(self.memberTab))
        self.Boxes[2][0].place(x=150, y=300)
        self.Boxes[2][1].place(x=350, y=300)
        self.Boxes[2][2].place(x=150, y=350)
        self.Boxes[2][3].place(x=350, y=350)
        self.Boxes[2][4].place(x=150, y=400)
        self.Boxes[2][5].place(x=350, y=400)
        self.Boxes[2][6].place(x=150, y=450)
        self.Boxes[2][7].place(x=350, y=450)
        self.Buttons.append([])
        self.Buttons[2].append(tk.Button(self.memberTab, text="Edit", command=self.editValues))
        self.Buttons[2].append(tk.Button(self.memberTab, text="Add", command=self.addValues))
        self.Buttons[2][0].place(x=550, y=300)
        self.Buttons[2][1].place(x=550, y=350)

        # Supports
        self.supportTab = ttk.Frame(self.tableTabs)
        self.tableTabs.add(self.supportTab, text="Supports")
        self.tables.append(self.createTable(self.supportTab,("Index", "i Node", "D X", "D Y", "D Z.", "R X", "R Y", "R Z")))
        self.Boxes.append([])
        for i in range(7): self.Boxes[3].append(tk.Entry(self.supportTab))
        self.Boxes[3][0].place(x=150, y=300)
        self.Boxes[3][1].place(x=350, y=300)
        self.Boxes[3][2].place(x=150, y=350)
        self.Boxes[3][3].place(x=350, y=350)
        self.Boxes[3][4].place(x=150, y=400)
        self.Boxes[3][5].place(x=350, y=400)
        self.Boxes[3][6].place(x=150, y=450)
        self.Buttons.append([])
        self.Buttons[3].append(tk.Button(self.supportTab, text="Edit", command=self.editValues))
        self.Buttons[3].append(tk.Button(self.supportTab, text="Add", command=self.addValues))
        self.Buttons[3][0].place(x=550, y=300)
        self.Buttons[3][1].place(x=550, y=350)

        # Releases
        self.releaseTab = ttk.Frame(self.tableTabs)
        self.tableTabs.add(self.releaseTab, text="Releases")
        self.tables.append(self.createTable(self.releaseTab,("Index", "i Member", "i DX", "i DY", "i DZ", "i RX", "i RY", "i RZ", "j DX", "j DY", "j DZ", "j RX", "j RY", "j RZ")))
        self.Boxes.append([])
        for i in range(13): self.Boxes[4].append(tk.Entry(self.releaseTab))
        self.Boxes[4][0].place(x=150, y=300)
        self.Boxes[4][1].place(x=350, y=300)
        self.Boxes[4][2].place(x=150, y=350)
        self.Boxes[4][3].place(x=350, y=350)
        self.Boxes[4][4].place(x=150, y=400)
        self.Boxes[4][5].place(x=350, y=400)
        self.Boxes[4][6].place(x=150, y=450)
        self.Boxes[4][7].place(x=350, y=450)
        self.Boxes[4][8].place(x=150, y=500)
        self.Boxes[4][9].place(x=350, y=500)
        self.Boxes[4][10].place(x=150, y=550)
        self.Boxes[4][11].place(x=350, y=550)
        self.Boxes[4][12].place(x=150, y=600)
        self.Buttons.append([])
        self.Buttons[4].append(tk.Button(self.releaseTab, text="Edit", command=self.editValues))
        self.Buttons[4].append(tk.Button(self.releaseTab, text="Add", command=self.addValues))
        self.Buttons[4][0].place(x=550, y=300)
        self.Buttons[4][1].place(x=550, y=350)

        # Nodal Loads
        self.nodeLoadTab = ttk.Frame(self.tableTabs)
        self.tableTabs.add(self.nodeLoadTab, text="Node Loads")
        self.tables.append(self.createTable(self.nodeLoadTab,("Index", "i Node", "P X", "P Y", "P Z.", "M X", "M Y", "M Z","Casee")))
        self.Boxes.append([])
        for i in range(8): self.Boxes[5].append(tk.Entry(self.nodeLoadTab))
        self.Boxes[5][0].place(x=150, y=300)
        self.Boxes[5][1].place(x=350, y=300)
        self.Boxes[5][2].place(x=150, y=350)
        self.Boxes[5][3].place(x=350, y=350)
        self.Boxes[5][4].place(x=150, y=400)
        self.Boxes[5][5].place(x=350, y=400)
        self.Boxes[5][6].place(x=150, y=450)
        self.Boxes[5][7].place(x=350, y=450)
        self.Buttons.append([])
        self.Buttons[5].append(tk.Button(self.nodeLoadTab, text="Edit", command=self.editValues))
        self.Buttons[5].append(tk.Button(self.nodeLoadTab, text="Add", command=self.addValues))
        self.Buttons[5][0].place(x=550, y=300)
        self.Buttons[5][1].place(x=550, y=350)

        # Member Point Loads
        self.memberPointLoadTab = ttk.Frame(self.tableTabs)
        self.tableTabs.add(self.memberPointLoadTab, text="Member Point Loads")
        self.tables.append(self.createTable(self.memberPointLoadTab, ("Index", "i Member", "x", "P X", "P Y", "P Z.", "M X", "M Y", "M Z", "Casee")))
        self.Boxes.append([])
        for i in range(9): self.Boxes[6].append(tk.Entry(self.memberPointLoadTab))
        self.Boxes[6][0].place(x=150, y=300)
        self.Boxes[6][1].place(x=350, y=300)
        self.Boxes[6][2].place(x=150, y=350)
        self.Boxes[6][3].place(x=350, y=350)
        self.Boxes[6][4].place(x=150, y=400)
        self.Boxes[6][5].place(x=350, y=400)
        self.Boxes[6][6].place(x=150, y=450)
        self.Boxes[6][7].place(x=350, y=450)
        self.Boxes[6][8].place(x=150, y=500)
        self.Buttons.append([])
        self.Buttons[6].append(tk.Button(self.memberPointLoadTab, text="Edit", command=self.editValues))
        self.Buttons[6].append(tk.Button(self.memberPointLoadTab, text="Add", command=self.addValues))
        self.Buttons[6][0].place(x=550, y=300)
        self.Buttons[6][1].place(x=550, y=350)

        # Member Distributed Loads
        self.memberDistLoadTab = ttk.Frame(self.tableTabs)
        self.tableTabs.add(self.memberDistLoadTab, text="Member Distributed Loads")
        self.tables.append(self.createTable(self.memberDistLoadTab,("Index", "i Member", "x 1", "x 2", "wx 1", "wx 2", "wy 1", "wy 2", "wz 1", "wz 2", "Casee")))
        self.Boxes.append([])
        for i in range(10): self.Boxes[7].append(tk.Entry(self.memberDistLoadTab))
        self.Boxes[7][0].place(x=150, y=300)
        self.Boxes[7][1].place(x=350, y=300)
        self.Boxes[7][2].place(x=150, y=350)
        self.Boxes[7][3].place(x=350, y=350)
        self.Boxes[7][4].place(x=150, y=400)
        self.Boxes[7][5].place(x=350, y=400)
        self.Boxes[7][6].place(x=150, y=450)
        self.Boxes[7][7].place(x=350, y=450)
        self.Boxes[7][8].place(x=150, y=500)
        self.Boxes[7][9].place(x=350, y=500)
        self.Buttons.append([])
        self.Buttons[7].append(tk.Button(self.memberDistLoadTab, text="Edit", command=self.editValues))
        self.Buttons[7].append(tk.Button(self.memberDistLoadTab, text="Add", command=self.addValues))
        self.Buttons[7][0].place(x=550, y=300)
        self.Buttons[7][1].place(x=550, y=350)

    def createTable(self, tab, Headings):
        table = ttk.Treeview(tab, columns=Headings, show='headings', height=10)

        scrollbar = ttk.Scrollbar(tab, orient="vertical", command=table.yview)
        scrollbar.pack(side="right", fill="y")

        table.configure(yscrollcommand=scrollbar.set)

        for i in range(len(Headings)):
            table.heading(Headings[i], text=Headings[i])
            table.column(Headings[i], width=50)

        table.bind("<<TreeviewSelect>>", self.getSelectedRow)

        table.pack()

        return table

    # noinspection PyUnusedLocal
    def getSelectedRow(self, event):
        numCol = [3,5,8,7,13,8,9,10]
        t_id = self.tableTabs.index("current")

        row_id = self.tables[t_id].focus()
        row_info = self.tables[t_id].item(row_id).get('values')
        for i in range(numCol[t_id]):
            self.Boxes[t_id][i].delete(0, 'end')
            self.Boxes[t_id][i].insert(0, row_info[i+1])

    def editValues(self):
        numCol = [3, 5, 8, 7, 13, 8, 9, 10]
        t_id = self.tableTabs.index("current")

        row_id = self.tables[t_id].focus()
        row_info = self.tables[t_id].item(row_id).get('values')
        newVal = [row_info[0]]
        for i in range(numCol[t_id]): newVal.append(self.Boxes[t_id][i].get()) # geting values from text boxes
        newVal = np.array(newVal)

        # Updating data
        if t_id == 0:   self.data.editnode(self, newVal[1:], newVal[[0]], True, True, row_id)
        elif t_id == 1: self.data.editmaterials(self, newVal[1:], newVal[[0]], True, True, row_id)
        elif t_id == 2: self.data.editmembers(self, newVal[1:], newVal[[0]], True, True, row_id)
        elif t_id == 3: self.data.editsupports(self, newVal[1:], newVal[[0]], True, True, row_id)
        elif t_id == 4: self.data.editreleases(self, newVal[1:], newVal[[0]], True, True, row_id)
        elif t_id == 5: self.data.editnodeLoads(self, newVal[1:], newVal[[0]], True, True, row_id)
        elif t_id == 6: self.data.editmemberPointLoad(self, newVal[1:], newVal[[0]], True, True, row_id)
        elif t_id == 7: self.data.editmemberDistLoad(self, newVal[1:], newVal[[0]], True, True, row_id)

    def addValues(self):
        numCol = [3, 5, 8, 7, 13, 8, 9, 10]
        t_id = self.tableTabs.index("current")
        numRow = len(self.tables[t_id].get_children())
        newVal = [numRow]
        for i in range(numCol[t_id]): newVal.append(self.Boxes[t_id][i].get()) # geting values from text boxes

        data_new = []
        for i in range(numCol[t_id]):
            if newVal[i + 1] == "True" or newVal[i + 1] == "False":
                data_new.append(bool(newVal[i + 1]))
            else:
                data_new.append(float(newVal[i + 1]))

        # Updating data
        if t_id == 0: self.data.addnodes(self, data_new, True, True)
        elif t_id == 1: self.data.addmaterials(self, data_new, True, True)
        elif t_id == 2: self.data.addmembers(self, data_new, True, True)
        elif t_id == 3: self.data.addsupports(self, data_new, True, True)
        elif t_id == 4: self.data.addreleases(self, data_new, True, True)
        elif t_id == 5: self.data.addnodeLoads(self, data_new, True, True)
        elif t_id == 6: self.data.addmemberPointLoads(self, data_new, True, True)
        elif t_id == 7: self.data.addmemberDistLoads(self, data_new, True, True)

    def centerWindow(self):
        width = 1000
        height = 800

        screen_width = self.root.winfo_screenwidth()
        screen_height = self.root.winfo_screenheight()

        if width > screen_width or height > screen_height:
            width = screen_width * 0.8
            height = screen_height * 0.8
        x = (screen_width // 2) - (width // 2)
        y = (screen_height // 2) - (height // 2) - 50

        self.root.geometry(f"{width}x{height}+{x}+{y}")

    def create_top_menu(self):
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)

        File_menu = tk.Menu(menubar, tearoff=False)
        File_menu.add_command(label='Save As', command=self.saveAs)
        File_menu.add_command(label='Open Frame', command=self.openFrame)
        menubar.add_cascade(label="File", menu=File_menu)

        import_menu = tk.Menu(menubar, tearoff=False)
        import_menu.add_command(label='Nodes', command= self.InportNodes)
        import_menu.add_command(label='Members', command=self.InportMembers)
        import_menu.add_command(label='Materials', command=self.InportMaterials)
        import_menu.add_command(label='Supports', command=self.ImportSupports)
        import_menu.add_command(label='Releases', command=self.ImportReleases)
        import_menu.add_command(label='Node Loads', command=self.ImportNodeLoads)
        import_menu.add_command(label='Members Point Loads', command=self.ImportMemberPointLoads)
        import_menu.add_command(label='Members Distributed Loads', command=self.ImportMemberDistLoads)
        menubar.add_cascade(label="Import", menu=import_menu)

        Analysis_menu = tk.Menu(menubar, tearoff=False)
        Analysis_menu.add_command(label='Linear Analysis', command=self.LinearAnalysis)
        Analysis_menu.add_command(label='Global Optimization', command=self.otimizationWindow)
        menubar.add_cascade(label="Analysis", menu=Analysis_menu)

    def saveAs(self):
        saveFrame(self.data)

    def openFrame(self):
        openFrame(self)

    def otimizationWindow(self):
        OptimizationPopUp(self, self.root, self.data)

    def closeProgram(self):
        self.root.destroy()

    def InportNodes(self):
        file_Types = [("Excel Files", "*.xlsx")]
        filePath = select_file_gui(file_Types)
        Nodes_df = pd.read_excel(filePath)
        Nodes = [Nodes_df["X"].tolist(), Nodes_df["Y"].tolist(), Nodes_df["Z"].tolist()]
        self.data.addnodes(self, Nodes, True, True)

    def InportMaterials(self):
        file_Types = [("Excel Files", "*.xlsx")]
        filePath = select_file_gui(file_Types)
        M_df = pd.read_excel(filePath)
        Materials = [M_df["E"].tolist(),M_df["G"].tolist(),M_df["nu"].tolist(),M_df["rho"].tolist(),M_df["fy"].tolist()]
        self.data.addmaterials(self, Materials, True, True)

    def InportMembers(self):
        file_Types = [("Excel Files", "*.xlsx")]
        filePath = select_file_gui(file_Types)
        M_df = pd.read_excel(filePath)
        Members = [M_df["i Node"].tolist(), M_df["j Node"].tolist(), M_df["Material"].tolist(),
                   M_df["Set Cross-Section Properties"].tolist(), M_df["A"].tolist(), M_df["Iy"].tolist(),
                   M_df["Iz"].tolist(), M_df["J"].tolist()]
        self.data.addmembers(self, Members, True, True)

    def ImportSupports(self):
        file_Types = [("Excel Files", "*.xlsx")]
        filePath = select_file_gui(file_Types)
        S_df = pd.read_excel(filePath)
        Supports = [S_df["Node"].tolist(),S_df["DX"].tolist(),S_df["DY"].tolist(),S_df["DZ"].tolist(),
                    S_df["RX"].tolist(),S_df["RY"].tolist(),S_df["RZ"].tolist()]
        self.data.addsupports(self, Supports, True, True)

    def ImportReleases(self):
        file_Types = [("Excel Files", "*.xlsx")]
        filePath = select_file_gui(file_Types)
        R_df = pd.read_excel(filePath)
        Releases = [R_df["Member"].tolist(), R_df["i DX"].tolist(), R_df["i DY"].tolist(), R_df["i DZ"].tolist(),
                    R_df["i RX"].tolist(), R_df["i RY"].tolist(), R_df["i RZ"].tolist(),R_df["j DX"].tolist(),
                    R_df["j DY"].tolist(), R_df["j DZ"].tolist(), R_df["j RX"].tolist(), R_df["j RY"].tolist(),
                    R_df["j RZ"].tolist()]
        self.data.addreleases(self, Releases, True, True)

    def ImportNodeLoads(self):
        file_Types = [("Excel Files", "*.xlsx")]
        filePath = select_file_gui(file_Types)
        N_df = pd.read_excel(filePath)
        NodeLoads = [N_df["Node"].tolist(),N_df["PX"].tolist(),N_df["PY"].tolist(),N_df["PZ"].tolist(),
                     N_df["MX"].tolist(),N_df["MY"].tolist(),N_df["MZ"].tolist()]
        self.data.addnodeLoads(self, NodeLoads, True, True)

    def ImportMemberPointLoads(self):
        file_Types = [("Excel Files", "*.xlsx")]
        filePath = select_file_gui(file_Types)
        M_df = pd.read_excel(filePath)
        MemberPointLoads = [M_df["Member"].tolist(),M_df["X"].tolist(),M_df["PX"].tolist(),M_df["PY"].tolist(),
                            M_df["PZ"].tolist(),M_df["MX"].tolist(),M_df["MY"].tolist(),M_df["MZ"].tolist()]
        self.data.addmemberPointLoads(self, MemberPointLoads, True, True)

    def ImportMemberDistLoads(self):
        file_Types = [("Excel Files", "*.xlsx")]
        filePath = select_file_gui(file_Types)
        M_df = pd.read_excel(filePath)
        MemberDistLoads = [M_df["Member"].tolist(),M_df["X1"].tolist(),M_df["X2"].tolist(),M_df["WX1"].tolist(),
                           M_df["WX2"].tolist(),M_df["WY1"].tolist(),M_df["WY2"].tolist(),M_df["WZ1"].tolist(),
                           M_df["WZ2"].tolist()]
        self.data.addmemberDistLoads(self, MemberDistLoads, True, True)

    def LinearAnalysis(self):
        Frame = addDataToFrame(self.data)

        Frame.preAnalysis_linear()
        D, DX, DY, DZ, RX, RY, RZ, weight, reactions, internalForces = Frame.analysis_linear(getWeight = True, getReactions = True, getInternalForces = True)
        self.Frame = Frame
        self.results = Results()
        self.results.addNodalDeflections(D, DX, DY, DZ, RX, RY, RZ)
        self.results.addWeight(weight)
        self.results.addReactions(reactions)
        self.results.addInternalForces(internalForces)
        saveFrame(self.data)
        saveResults(self.results)

    def GlobalOptimization(self, GroupAssignments, GroupTypes, lowerBound, upperBound, costFunction, weightRun, reactionRun, internalForcesRun):
        saveFrame(self.data)

        Frame = addDataToFrame(self.data)

        Frame.preAnalysis_linear()
        results = Frame.optimize(GroupAssignments, GroupTypes, lowerBound, upperBound, costFunction, weightRun, reactionRun, internalForcesRun)
        self.Frame = Frame
        self.Optimization_Results = results

        # todo: update frame object to optimization results
        self.LinearAnalysis()
        # todo: save optimization results

    def addPrintNode(self, node):
        self.printNodes = np.vstack((self.printNodes, node))
        self.updateCanves()

    def addPrintLine(self, line):
        l = []
        for i in range(len(line)):
            l.append(int(line[i]))
        l = np.array(l)
        self.printLines = np.vstack((self.printLines, [self.printNodes[l]]))
        self.updateCanves()

    def addPrintSurface(self, surface):
        tri = TDE.triangalizeSurface(self, surface, False)
        for i in range(len(tri)): self.surfaceTri = np.vstack((self.surfaceTri, tri[i]))
        self.updateCanves()

    def addPrintSolid(self, solid, flipNormal: list | None = None):
        numSurfaces = len(solid)
        for i in range(numSurfaces):
            if isinstance(flipNormal, list): # needs triangalized
                tri = TDE.triangalizeSurface(self, solid[i], flipNormal[i])
                for j in range(len(tri)): self.solidTri = np.vstack((self.solidTri, [tri[j]]))
            else: # already triangalized
                self.solidTri = np.vstack((self.solidTri, [solid[i]]))

    def scroolDown(self, event):
        if bool(event.state & 0x0001):
            self.shiftStateLast = True
            self.rotStart(event)
        else:
            self.shiftStateLast = False
            self.panStart(event)

    def scroolUpdate(self, event):
        if bool(event.state & 0x0001):
            if not self.shiftStateLast:
                self.shiftStateLast = True
                self.panEnd(event)
                self.rotStart(event)
            else: self.rotUpdate(event)
        else:
            if self.shiftStateLast:
                self.shiftStateLast = False
                self.rotEnd(event)
                self.panStart(event)
            else: self.panUpdate(event)

    def scroolUp(self, event):
        if bool(event.state & 0x0001):
            if not self.shiftStateLast: self.panEnd(event)
            else: self.rotEnd(event)
        else:
            if self.shiftStateLast: self.rotEnd(event)
            else: self.panEnd(event)

    def zoom(self, event):
        if event.delta > 0:
            self.camera = self.camera + self.lookDir * 0.1
        else:
            self.camera = self.camera - self.lookDir * 0.1
        self.updateCanves()

    def panStart(self, event):
        self.scroolCordsLast = [event.x, event.y]

    def panUpdate(self, event):
        cords = [event.x, event.y]
        self.panMove(cords)

    def panEnd(self, event):
        cords = [event.x, event.y]
        self.panMove(cords)

    def panMove(self, cords):
        changeCords = [cords[0] - self.scroolCordsLast[0], cords[1] - self.scroolCordsLast[1]]
        self.camera = self.camera - self.up * changeCords[1] * 0.01
        self.camera = self.camera + np.cross(self.lookDir, self.up) * changeCords[0] * 0.01
        self.updateCanves()
        self.scroolCordsLast = cords

    def rotStart(self, event):
        self.scroolCordsLast = [event.x, event.y]

    def rotUpdate(self, event):
        cords = [event.x, event.y]
        self.rotMove(cords)

    def rotEnd(self, event):
        cords = [event.x, event.y]
        self.rotMove(cords)

    def rotMove(self, cords):
        changeCords = [cords[0] - self.scroolCordsLast[0], cords[1] - self.scroolCordsLast[1]]
        self.cam_Y_rot = self.cam_Y_rot + changeCords[0] * 0.005
        self.cam_X_rot = self.cam_X_rot + changeCords[1] * 0.005
        self.updateCanves()
        self.scroolCordsLast = cords

    def updateCanves(self):
        node = copy(self.printNodes)
        line = copy(self.printLines)
        surfTri = copy(self.surfaceTri)
        solidTri = copy(self.solidTri)

        self.updateCamData()
        solidTriNormals = TDE.getNormals(solidTri)
        solidTri, solidTriNormals = TDE.removeTriFaceingAway(self, solidTri, solidTriNormals)
        triColor = TDE.illumination(self, solidTriNormals, len(surfTri))

        for i in range(len(surfTri)): solidTri = np.append(solidTri, surfTri[i], axis=0)
        tri = solidTri

        node, line, tri = TDE.transformToLocal(self, node, line, tri)

        node, line, tri, triColor = TDE.clipClose(node, line, tri, triColor)

        self.graph.delete("all")
        if len(node) > 0 or len(line) > 0 or len(tri) > 0:
            w = self.graph.winfo_width()
            h = self.graph.winfo_height()
            node, line, tri = TDE.project(self, node, line, tri)

            if len(tri) > 0: tri = tri[:, :, :-1]
            if len(line) > 0: line = line[:, :, :-1]
            node, line, tri, triColor = TDE.clipEadges(node, line, tri, triColor, h, w)

            if len(tri) > 0: self.printTri(np.array(tri), np.array(triColor))
            if len(line) > 0: self.printLine(np.array(line))
            if len(node) > 0: self.printNode(np.array(node))

        self.update()

    def updateCamData(self):
        matCameraRotY = TDE.getRotationYMatrix(self.cam_Y_rot)
        matCameraRotX = TDE.getRotationXMatrix(self.cam_X_rot)
        self.lookDir = (np.append(np.array([0, 0, 1]), 1) @ matCameraRotY @ matCameraRotX)[:-1]
        self.lookDir = self.lookDir / np.linalg.norm(self.lookDir)
        self.target = self.camera + self.lookDir
        self.up = (np.append(np.array([0, 1, 0]), 1) @ matCameraRotY @ matCameraRotX)[:-1]
        self.up = self.up / np.linalg.norm(self.up)

    def printTri(self, tri, triColor):
        avgDistance = np.sum(tri[:, :, 2], axis=-1)
        idx = np.argsort(avgDistance)
        tri = tri[idx[::-1]]
        triColor = triColor[idx[::-1]]

        triPrint = []
        for i in range(len(tri)):
            triPrint.append([tri[i][0][0], tri[i][0][1],
                             tri[i][1][0], tri[i][1][1],
                             tri[i][2][0], tri[i][2][1]])
        for i in range(len(triPrint)):
            color = int(triColor[i])
            self.graph.create_polygon(triPrint[i], outline='', fill="#%02x%02x%02x" % (color, color, color))

    def printLine(self, line):
        avgDistance = np.sum(line[:, :, 2], axis=-1)
        idx = np.argsort(avgDistance)
        line = line[idx[::-1]]

        linePrint = []
        for i in range(len(line)):
            linePrint.append([line[i][0][0], line[i][0][1],
                             line[i][1][0], line[i][1][1]])
        for i in range(len(linePrint)): self.graph.create_line(linePrint[i], width = 5, fill="red")

    def printNode(self, node):
        idx = np.argsort(node[:,2])
        node = node[idx[::-1]]

        for i in range(len(node)):
            nodeP = [node[i][0], node[i][1]]
            self.graph.create_oval(nodeP[0] - 4, nodeP[1] - 4, nodeP[0] + 4, nodeP[1] + 4, fill="green")


def select_file_gui(file_Types):
    """
    Opens a file explorer dialog using Tkinter and returns the selected file path.
    """
    # Create a Tk root window (but hide it)
    root = tk.Tk()
    root.withdraw()  # Hide the main window

    # Open the file dialog
    file_path = filedialog.askopenfilename(
        title="Select a file",
        filetypes= file_Types
    )

    if file_path:
        print(f"Selected file: {file_path}")
        return file_path
    else:
        print("No file selected.")
        return None

def addDataToFrame(d):
    Frame = Frame3D()

    for i in range(len(d.nodes[0])):
        Frame.addNode(d.nodes[0][i], d.nodes[1][i], d.nodes[2][i])

    for i in range(len(d.materials[0])):
        Frame.addMaterial(d.materials[0][i], d.materials[1][i], d.materials[2][i], d.materials[3][i], d.materials[4][i])

    for i in range(len(d.members[0])):
        Frame.addMember(d.members[0][i], d.members[1][i], d.members[2][i], d.members[3][i], d.members[4][i],
                        d.members[5][i], d.members[6][i], d.members[7][i])

    for i in range(len(d.supports[0])):
        Frame.defSupport(d.supports[0][i], d.supports[1][i], d.supports[2][i], d.supports[3][i], d.supports[4][i],
                         d.supports[5][i], d.supports[6][i])

    for i in range(len(d.releases[0])):
        Frame.defReleases(d.releases[0][i], d.releases[1][i], d.releases[2][i], d.releases[3][i], d.releases[4][i],
                          d.releases[5][i], d.releases[6][i], d.releases[7][i], d.releases[8][i],d.releases[9][i],
                          d.releases[10][i], d.releases[11][i], d.releases[12][i])

    for i in range(len(d.nodeLoad[0])):
        Frame.addNodeLoad(d.nodeLoad[0][i], d.nodeLoad[1][i], d.nodeLoad[2][i], d.nodeLoad[3][i], d.nodeLoad[4][i],
                          d.nodeLoad[5][i], d.nodeLoad[6][i], d.nodeLoad[7][i])

    for i in range(len(d.memberPointLoad[0])):
        Frame.addMemberPointLoad(d.memberPointLoad[0][i], d.memberPointLoad[1][i], d.memberPointLoad[2][i],
                                 d.memberPointLoad[3][i], d.memberPointLoad[4][i], d.memberPointLoad[5][i],
                                 d.memberPointLoad[6][i], d.memberPointLoad[7][i], d.memberPointLoad[8][i])

    for i in range(len(d.memberDistLoad[0])):
        Frame.addMemberDistLoad(d.memberDistLoad[0][i], d.memberDistLoad[1][i], d.memberDistLoad[2][i],
                                d.memberDistLoad[3][i], d.memberDistLoad[4][i], d.memberDistLoad[5][i],
                                d.memberDistLoad[6][i], d.memberDistLoad[7][i], d.memberDistLoad[8][i],
                                d.memberDistLoad[9][i])
    return Frame
import tkinter as tk
from tkinter import filedialog

import pandas as pd

from Frame3DGUI.inputData import data
from StructuralAnalysis.frame3DSolver.__main__ import Frame3D

class MainWindow(tk.Frame):
    def __init__(self, root):
        tk.Frame.__init__(self, root)

        self.root = root

        self.create_top_menu()

        self.centerWindow()

        exitButton = tk.Button(root, text="Complete Sketch", command=self.closeProgram)

        self.data = data()

        self.Frame = None
        self.Results = None

        #self.canvas.create_window(100, 100, window=exitButton)

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

        import_menu = tk.Menu(menubar, tearoff="off")
        import_menu.add_command(label='Nodes', command= self.InportNodes)
        import_menu.add_command(label='Members', command=self.InportMembers)
        import_menu.add_command(label='Materials', command=self.InportMaterials)
        import_menu.add_command(label='Supports', command=self.ImportSupports)
        import_menu.add_command(label='Releases', command=self.ImportReleases)
        import_menu.add_command(label='Node Loads', command=self.ImportNodeLoads)
        import_menu.add_command(label='Members Point Loads', command=self.ImportMemberPointLoads)
        import_menu.add_command(label='Members Distributed Loads', command=self.ImportMemberDistLoads)
        menubar.add_cascade(label="Import", menu=import_menu)

        Analysis_menu = tk.Menu(menubar, tearoff="off")
        Analysis_menu.add_command(label='Linear Analysis', command=self.LinearAnalysis)
        menubar.add_cascade(label="Analysis", menu=Analysis_menu)

    def closeProgram(self):
        self.root.destroy()

    def InportNodes(self):
        file_Types = [("Excel Files", "*.xlsx")]
        filePath = select_file_gui(file_Types)
        Nodes_df = pd.read_excel(filePath)
        Nodes = [Nodes_df["X"].tolist(), Nodes_df["Y"].tolist(), Nodes_df["Z"].tolist()]
        self.data.addnodes(Nodes)

        self.Frame = None
        self.Results = None

    def InportMaterials(self):
        file_Types = [("Excel Files", "*.xlsx")]
        filePath = select_file_gui(file_Types)
        M_df = pd.read_excel(filePath)
        Materials = [M_df["E"].tolist(),M_df["G"].tolist(),M_df["nu"].tolist(),M_df["rho"].tolist(),M_df["fy"].tolist()]
        self.data.addmaterials(Materials)

        self.Frame = None
        self.Results = None

    def InportMembers(self):
        file_Types = [("Excel Files", "*.xlsx")]
        filePath = select_file_gui(file_Types)
        M_df = pd.read_excel(filePath)
        Members = [M_df["i Node"].tolist(), M_df["j Node"].tolist(), M_df["Material"].tolist(),
                   M_df["Set Cross-Section Properties"].tolist(), M_df["A"].tolist(), M_df["Iy"].tolist(),
                   M_df["Iz"].tolist(), M_df["J"].tolist()]
        self.data.addmembers(Members)

        self.Frame = None
        self.Results = None

    def ImportSupports(self):
        file_Types = [("Excel Files", "*.xlsx")]
        filePath = select_file_gui(file_Types)
        S_df = pd.read_excel(filePath)
        Supports = [S_df["Node"].tolist(),S_df["DX"].tolist(),S_df["DY"].tolist(),S_df["DZ"].tolist(),
                    S_df["RX"].tolist(),S_df["RY"].tolist(),S_df["RZ"].tolist()]
        self.data.addsupports(Supports)

        self.Frame = None
        self.Results = None

    def ImportReleases(self):
        file_Types = [("Excel Files", "*.xlsx")]
        filePath = select_file_gui(file_Types)
        R_df = pd.read_excel(filePath)
        Releases = [R_df["Member"].tolist(), R_df["i DX"].tolist(), R_df["i DY"].tolist(), R_df["i DZ"].tolist(),
                    R_df["i RX"].tolist(), R_df["i RY"].tolist(), R_df["i RZ"].tolist(),R_df["j DX"].tolist(),
                    R_df["j DY"].tolist(), R_df["j DZ"].tolist(), R_df["j RX"].tolist(), R_df["j RY"].tolist(),
                    R_df["j RZ"].tolist()]
        self.data.addreleases(Releases)

        self.Frame = None
        self.Results = None

    def ImportNodeLoads(self):
        file_Types = [("Excel Files", "*.xlsx")]
        filePath = select_file_gui(file_Types)
        N_df = pd.read_excel(filePath)
        NodeLoads = [N_df["Node"].tolist(),N_df["PX"].tolist(),N_df["PY"].tolist(),N_df["PZ"].tolist(),
                     N_df["MX"].tolist(),N_df["MY"].tolist(),N_df["MZ"].tolist()]
        self.data.addnodeLoads(NodeLoads)

        self.Frame = None
        self.Results = None

    def ImportMemberPointLoads(self):
        file_Types = [("Excel Files", "*.xlsx")]
        filePath = select_file_gui(file_Types)
        M_df = pd.read_excel(filePath)
        MemberPointLoads = [M_df["Member"].tolist(),M_df["X"].tolist(),M_df["PX"].tolist(),M_df["PY"].tolist(),
                            M_df["PZ"].tolist(),M_df["MX"].tolist(),M_df["MY"].tolist(),M_df["MZ"].tolist()]
        self.data.addmemberPointLoads(MemberPointLoads)

        self.Frame = None
        self.Results = None

    def ImportMemberDistLoads(self):
        file_Types = [("Excel Files", "*.xlsx")]
        filePath = select_file_gui(file_Types)
        M_df = pd.read_excel(filePath)
        MemberDistLoads = [M_df["Member"].tolist(),M_df["X1"].tolist(),M_df["X2"].tolist(),M_df["WX1"].tolist(),
                           M_df["WX2"].tolist(),M_df["WY1"].tolist(),M_df["WY2"].tolist(),M_df["WZ1"].tolist(),
                           M_df["WZ2"].tolist()]
        self.data.addmemberDistLoads(MemberDistLoads)

        self.Frame = None
        self.Results = None


    def LinearAnalysis(self):
        Frame = addDataToFrame(self.data)

        Frame.preAnalysis_linear()
        results = Frame.analysis_linear(getWeight = True, getReactions = True, getInternalForces = True)
        self.Frame = Frame
        self.results = results

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
        Frame.addMaterial(d.materials[0][i], d.materials[1][i], d.materials[2][i],
                          d.materials[3][i], d.materials[4][i])

    for i in range(len(d.members[0])):
        Frame.addMember(d.members[0][i], d.members[1][i], d.members[2][i],
                        d.members[3][i], d.members[4][i], d.members[5][i],
                        d.members[6][i], d.members[7][i])

    for i in range(len(d.supports[0])):
        Frame.defSupport(d.supports[0][i], d.supports[1][i], d.supports[2][i],
                         d.supports[3][i], d.supports[4][i], d.supports[5][i],
                         d.supports[6][i])

    for i in range(len(d.releases[0])):
        Frame.defReleases(d.releases[0][i], d.releases[1][i], d.releases[2][i],
                          d.releases[3][i], d.releases[4][i], d.releases[5][i],
                          d.releases[6][i], d.releases[7][i], d.releases[8][i],
                          d.releases[9][i], d.releases[10][i], d.releases[11][i],
                          d.releases[12][i])

    for i in range(len(d.nodeLoad[0])):
        Frame.addNodeLoad(d.nodeLoad[0][i], d.nodeLoad[1][i], d.nodeLoad[2][i],
                          d.nodeLoad[3][i], d.nodeLoad[4][i], d.nodeLoad[5][i],
                          d.nodeLoad[6][i], d.nodeLoad[7][i])

    for i in range(len(d.memberPointLoad[0])):
        Frame.addMemberPointLoad(d.memberPointLoad[0][i], d.memberPointLoad[1][i],
                                 d.memberPointLoad[2][i], d.memberPointLoad[3][i],
                                 d.memberPointLoad[4][i], d.memberPointLoad[5][i],
                                 d.memberPointLoad[6][i], d.memberPointLoad[7][i],
                                 d.memberPointLoad[8][i])

    for i in range(len(d.memberDistLoad[0])):
        Frame.addMemberDistLoad(d.memberDistLoad[0][i], d.memberDistLoad[1][i],
                                d.memberDistLoad[2][i], d.memberDistLoad[3][i],
                                d.memberDistLoad[4][i], d.memberDistLoad[5][i],
                                d.memberDistLoad[6][i], d.memberDistLoad[7][i],
                                d.memberDistLoad[8][i], d.memberDistLoad[9][i])
    return Frame
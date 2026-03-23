import numpy as np
import pandas as pd

from TopologyOptimization.Truss.optimize import OptimizeTruss
from CrossSectionOptimization.Truss.optimize import TrussMain

def openTrussTopologyOptimizationExcel(filePath):
    Boundary_df = pd.read_excel(filePath, sheet_name='Boundary')
    Noeds_df = pd.read_excel(filePath, sheet_name='Nodes')

    LoadCases_df = pd.read_excel(filePath, sheet_name='LoadCases')
    Supports_df = pd.read_excel(filePath, sheet_name='Supports')

    Boundary = [Boundary_df["X"].tolist(), Boundary_df["Y"].tolist()]
    Nodes = [Noeds_df["X"].tolist(), Noeds_df["Y"].tolist()]
    Supports = [Supports_df["X"].tolist(), Supports_df["Y"].tolist(), Supports_df["Support x"].tolist(), Supports_df["Support y"].tolist()]

    numLoadCases = LoadCases_df["Number Load Casses"].tolist()[0]

    LoadCasses = []
    for i in range(int(numLoadCases)):
        x = LoadCases_df[f"LoadCase{i+1} x"].tolist()
        y = LoadCases_df[f"LoadCase{i+1} y"].tolist()
        fx = LoadCases_df[f"LoadCase{i+1} fx"].tolist()
        fy = LoadCases_df[f"LoadCase{i+1} fy"].tolist()

        dis1x = LoadCases_df[f"LoadCase{i+1} dis1x"].tolist()
        dis1y = LoadCases_df[f"LoadCase{i+1} dis1y"].tolist()
        dis2x = LoadCases_df[f"LoadCase{i + 1} dis2x"].tolist()
        dis2y = LoadCases_df[f"LoadCase{i + 1} dis2y"].tolist()
        disfx = LoadCases_df[f"LoadCase{i + 1} disfx"].tolist()
        disfy = LoadCases_df[f"LoadCase{i + 1} disfy"].tolist()

        Case = np.array([x,y,fx,fy, dis1x, dis1y, dis2x, dis2y, disfx, disfy]).T.tolist()
        print(Case)
        LoadCasses.append(Case)

    Boundary = np.array(Boundary).T.tolist()
    Nodes = np.array(Nodes).T.tolist()
    Supports = np.array(Supports).T.tolist()

    OptimizeTruss(filePath, Boundary, LoadCasses, Supports, nodes = Nodes)

def openTrussCrossSectionOptimizationExcel(filePath):
    Members_df = pd.read_excel(filePath, sheet_name='Members')
    Nodes_df = pd.read_excel(filePath, sheet_name='Nodes')

    LoadCases_df = pd.read_excel(filePath, sheet_name='LoadCases')
    Supports_df = pd.read_excel(filePath, sheet_name='Supports')

    Nodes = [Nodes_df["X"].tolist(), Nodes_df["Y"].tolist()]
    Members = [Members_df["Node1"].tolist(), Members_df["Node2"].tolist()]
    Supports = [Supports_df["X"].tolist(), Supports_df["Y"].tolist(), Supports_df["Support x"].tolist(),
                Supports_df["Support y"].tolist()]

    numLoadCases = LoadCases_df["Number Load Casses"].tolist()[0]

    loadCasses = []
    for i in range(int(numLoadCases)):
        x = LoadCases_df[f"LoadCase{i + 1} x"].tolist()
        y = LoadCases_df[f"LoadCase{i + 1} y"].tolist()
        fx = LoadCases_df[f"LoadCase{i + 1} fx"].tolist()
        fy = LoadCases_df[f"LoadCase{i + 1} fy"].tolist()

        Case = np.array([x, y, fx, fy]).T.tolist()
        loadCasses.append(Case)

    nodes = np.array(Nodes).T.tolist()
    members = np.array(Members).T.tolist()
    supports = np.array(Supports).T.tolist()

    TrussMain(filePath, nodes, members, loadCasses, supports)

def openFrameCrossSectionOptimization(filePath):
    print("Opening")
    Members_df = pd.read_excel(filePath, sheet_name='Members')
    Nodes_df = pd.read_excel(filePath, sheet_name='Nodes')

    Loads_df = pd.read_excel(filePath, sheet_name='Loads')
    Supports_df = pd.read_excel(filePath, sheet_name='Supports')

    CrossSections_df = pd.read_excel(filePath, sheet_name='CrossSections')

    Other_df = pd.read_excel(filePath, sheet_name='Other')

    Nodes = [Nodes_df["X"].tolist(), Nodes_df["Y"].tolist()]
    Members = [Members_df["Node1"].tolist(), Members_df["Node2"].tolist()]
    Supports = [Supports_df["Node"].tolist(), Supports_df["Support x"].tolist(),
                Supports_df["Support y"].tolist(), Supports_df["Support m"].tolist()]
    Loads = [Loads_df["Node"].tolist(), Loads_df["x"].tolist(), Loads_df["y"].tolist(), Loads_df["m"].tolist()]

    CrossSectionNames = CrossSections_df["Name"].tolist()
    CrossSections = [CrossSections_df["A"].tolist(), CrossSections_df["I"].tolist(), CrossSections_df["Density"].tolist()]

    E = Other_df["E"].tolist()[0]
    A_0 = Other_df["Small A"].tolist()[0]
    I_0 = Other_df["Small I"].tolist()[0]

    Nodes = np.array(Nodes).T
    Members = np.array(Members).T
    Supports = np.array(Supports).T
    Loads = np.array(Loads).T
    CrossSections = np.array(CrossSections).T
    CrossSectionNames = np.array(CrossSectionNames)

    Length = ((Nodes[Members[:,1]][:,0] - Nodes[Members[:,0]][:,0])**2+(Nodes[Members[:,1]][:,1] - Nodes[Members[:,0]][:,1])**2)**0.5
    Members = np.column_stack((Members, Length))

    FrameMain(Nodes, Members, Loads, Supports, CrossSections, CrossSectionNames, E, A_0, I_0, 0)
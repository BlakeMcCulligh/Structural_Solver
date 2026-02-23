import numpy as np
import pandas as pd

from TopologyOptimization.Truss.optimize import OptimizeTruss
from CrossSectionOptimization.Truss.optimize import main

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

    main(filePath, nodes, members, loadCasses, supports)






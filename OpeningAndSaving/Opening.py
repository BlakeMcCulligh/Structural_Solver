import numpy as np
import pandas as pd

from TopologyOptimization.__main__ import OptimizeTruss

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

    OptimizeTruss(Boundary, LoadCasses, Supports, nodes = Nodes)
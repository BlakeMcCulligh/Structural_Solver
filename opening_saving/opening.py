"""
Opens Excel files to run analyses and optimizations on the structure they store.
"""

import numpy as np
import pandas as pd

from TopologyOptimization.Truss.optimize import OptimizeTruss
from CrossSectionOptimization.Truss.optimize import TrussMain

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

def open_truss_topology_optimization_excel(file_path):
    """
    Opens Excel file constaining a truss and runs a truss topology optimization on the truss.

    :param file_path: File path to the Excel file.
    """

    boundary_df = pd.read_excel(file_path, sheet_name='Boundary')
    nodes_df = pd.read_excel(file_path, sheet_name='Nodes')

    load_cases_df = pd.read_excel(file_path, sheet_name='LoadCases')
    supports_df = pd.read_excel(file_path, sheet_name='Supports')

    boundary = [boundary_df["X"].tolist(), boundary_df["Y"].tolist()]
    nodes = [nodes_df["X"].tolist(), nodes_df["Y"].tolist()]
    supports = [supports_df["X"].tolist(), supports_df["Y"].tolist(), supports_df["Support x"].tolist(), supports_df["Support y"].tolist()]

    num_load_cases = load_cases_df["Number Load Casses"].tolist()[0]

    LoadCasses = []
    for i in range(int(num_load_cases)):
        x = load_cases_df[f"LoadCase{i+1} x"].tolist()
        y = load_cases_df[f"LoadCase{i+1} y"].tolist()
        fx = load_cases_df[f"LoadCase{i+1} fx"].tolist()
        fy = load_cases_df[f"LoadCase{i+1} fy"].tolist()

        dis1x = load_cases_df[f"LoadCase{i+1} dis1x"].tolist()
        dis1y = load_cases_df[f"LoadCase{i+1} dis1y"].tolist()
        dis2x = load_cases_df[f"LoadCase{i + 1} dis2x"].tolist()
        dis2y = load_cases_df[f"LoadCase{i + 1} dis2y"].tolist()
        disfx = load_cases_df[f"LoadCase{i + 1} disfx"].tolist()
        disfy = load_cases_df[f"LoadCase{i + 1} disfy"].tolist()

        Case = np.array([x,y,fx,fy, dis1x, dis1y, dis2x, dis2y, disfx, disfy]).T.tolist()
        LoadCasses.append(Case)

    boundary = np.array(boundary).T.tolist()
    nodes = np.array(nodes).T.tolist()
    supports = np.array(supports).T.tolist()

    OptimizeTruss(file_path, boundary, LoadCasses, supports, nodes = nodes)

def open_truss_cross_section_optimization_excel(file_path):
    """
    Opens Excel file constaining a truss and runs a truss cross-section optimization on the truss.

    :param file_path: File path to the Excel file.
    """

    members_df = pd.read_excel(file_path, sheet_name='Members')
    nodes_df = pd.read_excel(file_path, sheet_name='Nodes')

    load_cases_df = pd.read_excel(file_path, sheet_name='LoadCases')
    supports_df = pd.read_excel(file_path, sheet_name='Supports')

    nodes = [nodes_df["X"].tolist(), nodes_df["Y"].tolist()]
    members = [members_df["Node1"].tolist(), members_df["Node2"].tolist()]
    supports = [supports_df["X"].tolist(), supports_df["Y"].tolist(), supports_df["Support location"].tolist(),
                supports_df["Support y"].tolist()]

    num_load_cases = load_cases_df["Number Load Casses"].tolist()[0]

    load_casses = []
    for i in range(int(num_load_cases)):
        x = load_cases_df[f"LoadCase{i + 1} location"].tolist()
        y = load_cases_df[f"LoadCase{i + 1} y"].tolist()
        fx = load_cases_df[f"LoadCase{i + 1} fx"].tolist()
        fy = load_cases_df[f"LoadCase{i + 1} fy"].tolist()

        Case = np.array([x, y, fx, fy]).T.tolist()
        load_casses.append(Case)

    nodes = np.array(nodes).T.tolist()
    members = np.array(members).T.tolist()
    supports = np.array(supports).T.tolist()

    TrussMain(file_path, nodes, members, load_casses, supports)
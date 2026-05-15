"""
Saves the results from the Excel import analyses.
"""

import numpy as np
import pandas as pd

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

def save_truss_cross_section_optimization_excel(file_path, nodes, members, load_cases, supports, areas, deflections,
                                                forces, volume):
    """
    Saves the results form a truss cross-section optimization to the provided Excel file.

    :param file_path: File path to the Excel file to save to.
    :param nodes: Locations of the nodes of the truss. Shape: (# Nodes, 2)
    :param members: Indices of the Nodes that the members run between. Shape: (# members, 3: [node i, node j, length])
    :param load_cases: Load cases applied to the truss. Shape: (# Load Cases, # Loads per Case, 4: [x,y,fx,fy] )
    :param supports: Nodes that are supported and in what directions. Shape: (# Supported Nodes, 4: [x,y,sx,sy])
    :param areas: Areas of each member. Shape: (# Members)
    :param deflections: Deflections at each node. Shape: (2 * # Nodes * # Load Cases)
    :param forces: Axial forces within each member. Shape: (# Members * # Load Cases)
    :param volume: Overall volume of the truss.
    """

    nodes_df = _convert_nodes_to_df(nodes, deflections)
    members_df = _convert_members_to_df(members, areas, forces)
    supports_df = _convert_supports_to_df(supports)
    volume_df = pd.DataFrame({"Volume": [volume]})
    load_cases_df = _convert_load_cases_to_df(load_cases)
    _write_to_excel(file_path, nodes_df, members_df, supports_df, load_cases_df, volume_df)

def save_truss_topology_optimization_excel(file_path, nodes, members, load_cases, supports, areas, deflections,
                                           forces, volume):
    """
    Saves the results form a truss topology optimization to the provided Excel file.

    :param file_path: File path to the Excel file to save to.
    :param nodes: Locations of the nodes of the truss. shape: (# Nodes, 2)
    :param members: Indices of the Nodes that the members run between. Shape: (# members, 3: [node i, node j, length])
    :param load_cases: Load cases applied to the truss. Shape: (# Load Cases, # Loads per Case, 4: [x1,y1,fx1,fy1] )
    :param supports: Nodes that are supported and in what directions. Shape: (# Supported Nodes, 4: [x,y,sx,sy])
    :param areas: Areas of each member. Shape: (# Members)
    :param deflections: Deflections at each node. Shape: (2 * # Nodes * # Load Cases)
    :param forces: Axial forces within each member. Shape: (# Members * # Load Cases)
    :param volume: Overall volume of the truss.
    """

    nodes_df = _convert_nodes_to_df(nodes, deflections)
    members_df = _convert_members_to_df(members, areas, forces, remove_small_members=True)
    supports_df = _convert_supports_to_df(supports)
    volume_df = pd.DataFrame({"Volume": [volume]})
    load_cases_df = _convert_load_cases_to_df(load_cases)
    _write_to_excel(file_path, nodes_df, members_df, supports_df, load_cases_df, volume_df)


def _convert_nodes_to_df(nodes, deflections):
    """
    Gets the data frame for the nodes sheet to be printed to the Excel.

    :param nodes: Locations of the nodes of the truss. shape: (# Nodes, 2)
    :param deflections: Deflections at each node. Shape: (2 * # Nodes * # Load Cases)
    :return: Nodes Data Frame
    """

    nodes = np.array(nodes).T.tolist()
    deflections = deflections.reshape(-1, len(nodes[0])).tolist()

    nodes_df = {"X": nodes[0], "Y": nodes[1]}
    for i, u in enumerate(deflections): nodes_df = nodes_df | {f"Delection{i + 1}": u}

    nodes_df = pd.DataFrame(nodes_df)
    return nodes_df

def _convert_members_to_df(members, areas, forces, remove_small_members=False):
    """
    Gets the data frame for the members sheet to be printed to the Excel.

    :param members: Indices of the Nodes that the members run between. Shape: (# members, 3: [node i, node j, length])
    :param areas: Areas of each member. Shape: (# Members)
    :param forces: Axial forces within each member. Shape: (# Members * # Load Cases)
    :param remove_small_members: If members under 0.001 times the cross-section area of the largest member should be
                                 removed.
    :return: Members Data Frame
    """

    max_area = max(areas)
    if remove_small_members:
        areas = areas.tolist()
        members = members.tolist()
        forces = forces.T.tolist()
        for i in reversed(range(len(members))):
            if areas[i] < max_area * 0.001:
                areas.pop(i)
                members.pop(i)
                forces.pop(i)
        forces = np.array(forces).T.tolist()

    members = np.array(members).T.tolist()

    members_df = {"Node1": members[0],"Node2": members[1],"Areas": areas}
    for i, q in enumerate(forces): members_df = members_df | {f"Force{i + 1}": q}

    members_df = pd.DataFrame(members_df)
    return members_df

def _convert_supports_to_df(supports):
    """
    Gets the data frame for the supports sheet to be printed to the Excel.

    :param supports: Nodes that are supported and in what directions. Shape: (# Supported Nodes, 4: [x,y,sx,sy])
    :return: Supports Data Frame
    """

    supports = np.array(supports).T.tolist()

    supports_df = {"X": supports[0],"Y": supports[1],"Support location": supports[2],"Support y": supports[3]}

    supports_df = pd.DataFrame(supports_df)
    return supports_df

def _convert_load_cases_to_df(load_cases):
    """
    Gets the data frame for the load cases sheet to be printed to the Excel.

    :param load_cases: Load cases applied to the truss. Shape: (# Load Cases, # Loads per Case, 4: [x,y,fx,fy] )
    :return: Load cases Data Frame
    """

    # TODO only works if all load cases have the same number of loads

    num_load_cases = len(load_cases)
    load_cases_df = {"Number Load cases": [num_load_cases] * len(load_cases[0])}
    for i, load_case in enumerate(load_cases):
        load_case = np.array(load_case).T.tolist()
        load = {f"LoadCase{i + 1} location": load_case[0],
                f"LoadCase{i + 1} y": load_case[1],
                f"LoadCase{i + 1} fx": load_case[2],
                f"LoadCase{i + 1} fy": load_case[3]}

        load_cases_df = load_cases_df | load
    load_cases_df = pd.DataFrame(load_cases_df)

    return load_cases_df

def _write_to_excel(file_path, nodes_df, members_df, supports_df, load_cases_df, volume_df):
    """
    Writes all the data frames to the Excel.

    :param file_path: File path to the Excel file to save to.
    :param nodes_df: Nodes Data Frame
    :param members_df: Members Data Frame
    :param supports_df: Supports Data Frame
    :param load_cases_df:  Load cases Data Frame
    :param volume_df: Volume Data Frame
    """

    with pd.ExcelWriter(file_path) as writer:
        nodes_df.to_excel(writer, sheet_name="Nodes", index=False)
        members_df.to_excel(writer, sheet_name="Members", index=False)
        supports_df.to_excel(writer, sheet_name="Supports", index=False)
        load_cases_df.to_excel(writer, sheet_name="LoadCases", index=False)
        volume_df.to_excel(writer, sheet_name="Volume", index=False)
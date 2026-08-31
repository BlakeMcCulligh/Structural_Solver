"""
Handels Analysis done by importing an Excel file and printing the results to an Excel file.
"""

from __future__ import annotations
import numpy as np
import pandas as pd
from typing import TYPE_CHECKING

from rebuild_program_layout.data_objects.material import Material
from rebuild_program_layout.data_objects.member import Member
from rebuild_program_layout.data_objects.node import Node
#from rebuild_program_layout.frame_3D import Frame3D

if TYPE_CHECKING:
    from rebuild_program_layout import __main__

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

def frame_3d_global_opt(file_path: str, controller: __main__.Structural_Solver, frame) -> None:
    """
    Optimizes the cross-sections of a 3D frame in the global space.

    :param file_path: File path to the Excel file to import
    :type file_path: String
    :param controller: Object containing the main program.
    :type controller: __main__.Structural_Solver
    """

    #frame = Frame3D(controller, None)

    # Nodes
    nodes_df = pd.read_excel(file_path, sheet_name='Nodes')
    nodes = [nodes_df["X"].tolist(), nodes_df["Y"].tolist(), nodes_df["Z"].tolist()]
    for i in range(len(nodes[0])): frame.nodes.append(Node(x=nodes[0][i], y=nodes[1][i], z=nodes[2][i]))

    # Materials
    mat_df = pd.read_excel(file_path, sheet_name='Materials')
    m = [mat_df["E"].tolist(), mat_df["G"].tolist(), mat_df["nu"].tolist(), mat_df["rho"].tolist(),
         mat_df["fy"].tolist()]
    for i in range(len(m[0])): frame.materials.append(Material(m[0][i], m[1][i], m[2][i], m[3][i], m[4][i]))

    # Members
    mem_df = pd.read_excel(file_path, sheet_name='Members')
    m = [mem_df["i Node"].tolist(), mem_df["j Node"].tolist(), mem_df["Material"].tolist(),
         mem_df["Set Cross-Section Properties"].tolist(), mem_df["A"].tolist(), mem_df["Iy"].tolist(),
         mem_df["Iz"].tolist(), mem_df["J"].tolist()]
    for i in range(len(m[0])): frame.members.append(Member(m[0][i], m[1][i], m[2][i], m[3][i],
                                              [m[4][i], m[5][i], m[6][i], m[7][i]], frame))

    # Supports
    s_df = pd.read_excel(file_path, sheet_name='Supports')
    supports = [s_df["Node"].tolist(), s_df["DX"].tolist(), s_df["DY"].tolist(), s_df["DZ"].tolist(),
                s_df["RX"].tolist(), s_df["RY"].tolist(), s_df["RZ"].tolist()]
    supports = np.array(supports).T
    for i in range(len(supports)):
        frame.nodes[supports[i, 0]].set_support(supports[i, 1:].tolist())

    # Releases
    r_df = pd.read_excel(file_path, sheet_name='Releases')
    releases = [r_df["Member"].tolist(), r_df["i DX"].tolist(), r_df["i DY"].tolist(), r_df["i DZ"].tolist(),
                r_df["i RX"].tolist(), r_df["i RY"].tolist(), r_df["i RZ"].tolist(), r_df["j DX"].tolist(),
                r_df["j DY"].tolist(), r_df["j DZ"].tolist(), r_df["j RX"].tolist(), r_df["j RY"].tolist(),
                r_df["j RZ"].tolist()]
    releases = np.array(releases).T
    for i in range(len(releases)):
        frame.members[releases[i, 0]].set_releces(releases[i, 1:].tolist())

    # Node Loads
    n_df = pd.read_excel(file_path, sheet_name='Node_Loads')
    node_loads = [n_df["Node"].tolist(), n_df["Case"].tolist(), n_df["PX"].tolist(), n_df["PY"].tolist(),
                  n_df["PZ"].tolist(), n_df["MX"].tolist(), n_df["MY"].tolist(), n_df["MZ"].tolist()]
    node_loads = np.array(node_loads).T
    for i in range(len(node_loads)):
        frame.nodes[node_loads[i, 0]].add_load(node_loads[i, 1].tolist(), node_loads[i, 2:].tolist())

    # Member Point Loads
    mp_df = pd.read_excel(file_path, sheet_name = 'Member_Point_Loads')
    member_point_loads = [mp_df["Member"].tolist(), mp_df["Case"].tolist(), mp_df["X"].tolist(), mp_df["PX"].tolist(),
                          mp_df["PY"].tolist(), mp_df["PZ"].tolist(), mp_df["MX"].tolist(), mp_df["MY"].tolist(),
                          mp_df["MZ"].tolist()]
    member_point_loads = np.array(member_point_loads).T
    for i in range(len(member_point_loads)):
        frame.members[int(member_point_loads[i, 0])].add_point_load(member_point_loads[i, 1].tolist(),
                                                              member_point_loads[i, 2].tolist(),
                                                              member_point_loads[i, 3:].tolist())

    # Member Distributed Loads
    md_df = pd.read_excel(file_path, sheet_name='Member_Dist_Loads')
    member_dist_loads = [md_df["Member"].tolist(), md_df["Case"].tolist(), md_df["X1"].tolist(), md_df["X2"].tolist(),
                         md_df["WX1"].tolist(), md_df["WX2"].tolist(), md_df["WY1"].tolist(), md_df["WY2"].tolist(),
                         md_df["WZ1"].tolist(), md_df["WZ2"].tolist()]
    member_dist_loads = np.array(member_dist_loads).T
    for i in range(len(member_dist_loads)):
        frame.members[int(member_dist_loads[i, 0])].add_dist_load(member_dist_loads[i, 1].tolist(),
                                                            member_dist_loads[i, 2:4].tolist(),
                                                            member_dist_loads[i, 4:].tolist())

    # Member Groups
    mga_df = pd.read_excel(file_path, sheet_name='Member_Group_Assignments')
    group_assignments = mga_df["Group_Number"].tolist()

    mgt_df = pd.read_excel(file_path, sheet_name='Member_Group_Types')
    group_types = mgt_df["Group Cross-Section Type"].tolist()
    d_low_bounds = mgt_df["Min d"].tolist()
    d_high_bounds = mgt_df["Max d"].tolist()
    b_low_bounds = mgt_df["Min b"].tolist()
    b_high_bounds = mgt_df["Max b"].tolist()
    t_low_bounds = mgt_df["Min t"].tolist()
    t_high_bounds = mgt_df["Max t"].tolist()

    lower_bounds = []
    upper_bounds = []
    for i in range(len(group_types)):
        lower_bounds.append(float(d_low_bounds[i]))
        upper_bounds.append(float(d_high_bounds[i]))
        if group_types[i] == "Angle" or group_types[i] == "RectHSS":
            lower_bounds.append(float(b_low_bounds[i]))
            upper_bounds.append(float(b_high_bounds[i]))
        lower_bounds.append(float(t_low_bounds[i]))
        upper_bounds.append(float(t_high_bounds[i]))

    # Cost Function
    cost_df = pd.read_excel(file_path, sheet_name='Cost_Function')
    method_of_defination = cost_df["Cost Function Name"].tolist()[0]

    if method_of_defination == "Other":
        cost_function = cost_df["Function"].tolist()[0]
        weight_run = cost_df["Weight Run"].tolist()[0]
        reaction_run = cost_df["Reaction Run"].tolist()[0]
        internal_forces_run = cost_df["Internal Forces Run"].tolist()[0]
    elif method_of_defination == "2027_Steel_Bridge":
        cost_function = "max(max(DZ)) + sum(Weight)"
        weight_run = True
        reaction_run = False
        internal_forces_run = False
    else:
        cost_function = ""
        weight_run = False
        reaction_run = False
        internal_forces_run = False

    import_data = [file_path, nodes_df, mem_df, mat_df, s_df, r_df, n_df, mp_df, md_df, mga_df, mgt_df, cost_df]

    frame.global_optimization(group_assignments, group_types, lower_bounds, upper_bounds, cost_function, weight_run,
                              reaction_run, internal_forces_run, import_data = import_data)

def results_frame_3d_global_opt(cost: float, results: list, import_data: list) -> None:
    """
    Writes the results of the 3D frame optimization of the global space to the Excel file.

    :param cost: Best cost found.
    :type cost: float
    :param results: The Cross-sections that resulted in the best cost.
    :type results: list.
    :param import_data: List of information that was imported.
    :type import_data: list
    """

    [file_path, nodes_df, mem_df, mat_df, s_df, r_df, n_df, mp_df, md_df, mga_df, mgt_df, cost_df] = import_data
    df1 = pd.DataFrame(results,
                       columns=["Member Index", "Cross-Section _type", "d", "b", "t", "A", "Iy", "Iz", "J"])
    df2 = pd.DataFrame({'Cost': [cost]})

    with pd.ExcelWriter(file_path) as writer:
        nodes_df.to_excel(writer, sheet_name="Nodes", index=False)
        mem_df.to_excel(writer, sheet_name="Members", index=False)
        mat_df.to_excel(writer, sheet_name="Materials", index=False)
        s_df.to_excel(writer, sheet_name="Supports", index=False)
        r_df.to_excel(writer, sheet_name="Releases", index=False)
        n_df.to_excel(writer, sheet_name="Node_Loads", index=False)
        mp_df.to_excel(writer, sheet_name="Member_Point_Loads", index=False)
        md_df.to_excel(writer, sheet_name="Member_Dist_Loads", index=False)
        mga_df.to_excel(writer, sheet_name="Member_Group_Assignments", index=False)
        mgt_df.to_excel(writer, sheet_name="Member_Group_Types", index=False)
        cost_df.to_excel(writer, sheet_name="Cost_Function", index=False)
        df1.to_excel(writer, sheet_name="Results", index=False)
        df2.to_excel(writer, sheet_name="Cost", index=False)

def frame_3d_global_opt_temp(file_path: str) -> None:
    """
    Exports a template Excel file forff the 3D frame global analysis.
    """

    nodes_df = pd.DataFrame([["","",""]], columns=["X", "Y", "Z"])
    mem_df = pd.DataFrame([["","","","","","","",""]],
                          columns=["i Node", "j Node", "Material", "Set Cross-Section Properties", "A", "Iy", "Iz",
                                   "J"])
    mat_df = pd.DataFrame([["","","","",""]], columns=["E", "G", "nu", "rho", "fy"])
    s_df = pd.DataFrame([["","","","","","",""]], columns=["Node", "DX", "DY", "DZ", "RX", "RY", "RZ"])
    r_df = pd.DataFrame([["","","","","","","","","","","","",""]],
                        columns=["Member","i DX", "i DY", "i DZ", "i RX", "i RY", "i RZ",
                                 "j DX", "j DY", "j DZ", "j RX", "j RY", "j RZ"])
    n_df = pd.DataFrame([["","","","","","","",""]], columns=["Node", "Case", "PX", "PY", "PZ", "MX", "MY", "MZ"])
    mp_df = pd.DataFrame([["","","","","","","","",""]],
                         columns=["Member", "Case", "X", "PX", "PY", "PZ", "MX", "MY", "MZ"])
    md_df = pd.DataFrame([["","","","","","","","","",""]],
                         columns=["Member", "Case", "X1", "X2", "WX1", "WX2", "WY1", "WY2", "WZ1", "WZ2"])
    mga_df = pd.DataFrame([[""]], columns=["Group_Number"])
    mgt_df = pd.DataFrame([["","","","","","",""]],
                          columns=["Group Cross-Section Type", "Min d", "Max d", "Min b", "Max b", "Min t", "Max t"])
    cost_df = pd.DataFrame([["","","","",""]],
                           columns=["Cost Function Name", "Function", "Weight Run", "Reaction Run",
                                    "Internal Forces Run"])

    with pd.ExcelWriter(file_path) as writer:
        nodes_df.to_excel(writer, sheet_name="Nodes", index=False)
        mem_df.to_excel(writer, sheet_name="Members", index=False)
        mat_df.to_excel(writer, sheet_name="Materials", index=False)
        s_df.to_excel(writer, sheet_name="Supports", index=False)
        r_df.to_excel(writer, sheet_name="Releases", index=False)
        n_df.to_excel(writer, sheet_name="Node_Loads", index=False)
        mp_df.to_excel(writer, sheet_name="Member_Point_Loads", index=False)
        md_df.to_excel(writer, sheet_name="Member_Dist_Loads", index=False)
        mga_df.to_excel(writer, sheet_name="Member_Group_Assignments", index=False)
        mgt_df.to_excel(writer, sheet_name="Member_Group_Types", index=False)
        cost_df.to_excel(writer, sheet_name="Cost_Function", index=False)
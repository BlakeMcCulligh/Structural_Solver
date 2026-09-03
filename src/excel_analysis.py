"""
Handels Analysis done by importing an Excel file and printing the results to an Excel file.
"""

from __future__ import annotations
import numpy as np
import pandas as pd
import tkinter as tk
from typing import TYPE_CHECKING

from data_objects.material import Material
from data_objects.member import Member
from data_objects.node import Node

if TYPE_CHECKING:
    from src import __main__
    from src.frame_3D import Frame3D
    from src.frame_2D import Frame2D

__author__ = "Blake McCulligh"
__copyright__ = "Copyright 2026 Blake McCulligh"
__credits__ = ["Blake McCulligh"]

__license__ = "MIT"
__version__ = "0.1.0b1"
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = "Beta"

def frame_3d_global_opt(file_path: str, controller: __main__.Structural_Solver, frame: Frame3D) -> None:
    """
    Optimizes the cross-sections of a 3D frame in the global space.

    :param file_path: File path to the Excel file to import
    :type file_path: String
    :param controller: Object containing the main program.
    :type controller: __main__.Structural_Solver
    :param frame: Frame to be optimized.
    :type frame: Frame3D
    """

    (group_assignments, group_types, lower_bounds, upper_bounds, inital, cost_function, weight_run,
    reaction_run, internal_forces_run, import_data) = _import_3D_opt_data_from_excel(file_path, frame)

    if group_assignments is not None:

        # noinspection PyBroadException
        try:
            # noinspection PyUnboundLocalVariable
            frame.global_optimization(group_assignments, group_types, lower_bounds, upper_bounds, cost_function,
                                      weight_run, reaction_run, internal_forces_run, import_data = import_data)
        except Exception as e:
            tk.messagebox.showinfo("Error", "An error ha occurred while running the analysis: {e}")

def frame_3d_local_opt(file_path: str, controller: __main__.Structural_Solver, frame: Frame3D) -> None:
    """
    Optimizes the cross-sections of a 3D frame in the local space.

    :param file_path: File path to the Excel file to import
    :type file_path: String
    :param controller: Object containing the main program.
    :type controller: __main__.Structural_Solver
    :param frame: Frame to be optimized.
    :type frame: Frame3D
    """

    (group_assignments, group_types, lower_bounds, upper_bounds, inital, cost_function, weight_run,
    reaction_run, internal_forces_run, import_data) = _import_3D_opt_data_from_excel(file_path, frame)

    if group_assignments is not None:

        # noinspection PyBroadException
        try:
            # noinspection PyUnboundLocalVariable
            frame.local_optimization(group_assignments, group_types, lower_bounds, upper_bounds, inital, cost_function,
                                      weight_run, reaction_run, internal_forces_run, import_data = import_data)
        except Exception as e:
            tk.messagebox.showinfo("Error", "An error ha occurred while running the analysis: {e}")

def results_frame_3d_opt(cost: float, results: list, import_data: list, optimization_type: str) -> None:
    """
    Writes the results of the 3D frame optimization of the global space to the Excel file.

    :param cost: Best cost found.
    :type cost: float
    :param results: The Cross-sections that resulted in the best cost.
    :type results: list.
    :param import_data: List of information that was imported.
    :type import_data: list
    :param optimization_type: Type of optimization that was performed.
    :type optimization_type: str
    """

    [file_path, nodes_df, mem_df, mat_df, s_df, r_df, n_df, mp_df, md_df, mga_df, mgt_df, cost_df] = import_data
    df1 = pd.DataFrame(results,
                       columns=["Member Index", "Cross-Section _type", "d", "b", "t", "A", "Iy", "Iz", "J"])

    if optimization_type == "Local":
        df2 = pd.DataFrame({'Cost': [cost], 'Optimization Type': [optimization_type]})
    else:
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

    tk.messagebox.showinfo("Analysis Complete", "Analysis Complete")

def _import_3D_opt_data_from_excel(file_path: str, frame: Frame3D):
    """
    Imports information form the specified Excel file

    :param file_path: File path to the Excel file to import
    :type file_path: String
    :param frame: Frame to be optimized.
    :type frame: Frame3D
    :return: group_assignments, group_types, lower_bounds, upper_bounds, inital, cost_function, weight_run,
                reaction_run, internal_forces_run, import_data
    """

    fail = False

    # Nodes
    nodes_df = pd.read_excel(file_path, sheet_name='Nodes')
    nodes = [nodes_df["Name"].tolist(), nodes_df["X"].tolist(), nodes_df["Y"].tolist(), nodes_df["Z"].tolist()]
    for i in range(len(nodes[0])): frame.nodes.append(Node(x=nodes[1][i], y=nodes[2][i], z=nodes[3][i]))

    # Materials
    mat_df = pd.read_excel(file_path, sheet_name='Materials')
    mat = [mat_df["Name"].tolist(), mat_df["E"].tolist(), mat_df["G"].tolist(), mat_df["nu"].tolist(),
           mat_df["rho"].tolist(), mat_df["fy"].tolist()]
    for i in range(len(mat[0])): frame.materials.append(Material(mat[1][i], mat[2][i], mat[3][i], mat[4][i], mat[5][i]))

    # Members
    mem_df = pd.read_excel(file_path, sheet_name='Members')
    mem = [mem_df["Name"].tolist(), mem_df["i Node"].tolist(), mem_df["j Node"].tolist(), mem_df["Material"].tolist(),
           mem_df["Set Cross-Section Properties"].tolist(), mem_df["A"].tolist(), mem_df["Iy"].tolist(),
           mem_df["Iz"].tolist(), mem_df["J"].tolist()]
    for i in range(len(mem[0])):
        nodei = None
        nodej = None
        material = None
        for j in range(len(nodes[0])):
            if nodes[0][j] == mem[1][i]:
                nodei = j
            elif nodes[0][j] == mem[2][i]:
                nodej = j
        for j in range(len(mat[0])):
            if mat[0][j] == mem[3][i]:
                material = j
        if nodei is not None and nodej is not None and material is not None:
            frame.members.append(Member(nodei, nodej, material, mem[4][i],
                                        [mem[5][i], mem[6][i], mem[7][i], mem[8][i]], frame))
        else:
            print("Node or Material not found. For member: ", mem[0][i])

    # Supports
    s_df = pd.read_excel(file_path, sheet_name='Supports')
    supports = [s_df["Node"].tolist(), s_df["DX"].tolist(), s_df["DY"].tolist(), s_df["DZ"].tolist(),
                s_df["RX"].tolist(), s_df["RY"].tolist(), s_df["RZ"].tolist()]
    supports = np.array(supports).T
    for i in range(len(supports)):
        node = None
        for j in range(len(nodes[0])):
            if nodes[0][j] == supports[i, 0]:
                node = j
        if node is not None:
            frame.nodes[node].set_support(supports[i, 1:].tolist())
        else:
            print("Support not found. For node: ", supports[i, 0])

    # Releases
    r_df = pd.read_excel(file_path, sheet_name='Releases')
    releases = [r_df["Member"].tolist(), r_df["i DX"].tolist(), r_df["i DY"].tolist(), r_df["i DZ"].tolist(),
                r_df["i RX"].tolist(), r_df["i RY"].tolist(), r_df["i RZ"].tolist(), r_df["j DX"].tolist(),
                r_df["j DY"].tolist(), r_df["j DZ"].tolist(), r_df["j RX"].tolist(), r_df["j RY"].tolist(),
                r_df["j RZ"].tolist()]
    releases = np.array(releases).T
    for i in range(len(releases)):
        member = None
        for j in range(len(mem[0])):
            if mem[0][j] == releases[i, 0]:
                member = j
        if member is not None:
            frame.members[member].set_releces(releases[i, 1:].tolist())
        else:
            print("Mmeber not found. For Release: ", releases[i, 0])

    # Node Loads
    n_df = pd.read_excel(file_path, sheet_name='Node_Loads')
    node_loads = [n_df["Node"].tolist(), n_df["Case"].tolist(), n_df["PX"].tolist(), n_df["PY"].tolist(),
                  n_df["PZ"].tolist(), n_df["MX"].tolist(), n_df["MY"].tolist(), n_df["MZ"].tolist()]
    node_loads = np.array(node_loads).T
    for i in range(len(node_loads)):
        node = None
        for j in range(len(nodes[0])):
            if nodes[0][j] == node_loads[i, 0]:
                node = j
        if node is not None:
            frame.nodes[node].add_load(node_loads[i, 1].tolist(), node_loads[i, 2:].tolist())
        else:
            print("Support not found. For Node Load: ", node_loads[i, 0])

    # Member Point Loads
    mp_df = pd.read_excel(file_path, sheet_name='Member_Point_Loads')
    member_point_loads = [mp_df["Member"].tolist(), mp_df["Case"].tolist(), mp_df["X"].tolist(), mp_df["PX"].tolist(),
                          mp_df["PY"].tolist(), mp_df["PZ"].tolist(), mp_df["MX"].tolist(), mp_df["MY"].tolist(),
                          mp_df["MZ"].tolist()]
    member_point_loads = np.array(member_point_loads).T
    for i in range(len(member_point_loads)):
        member = None
        for j in range(len(mem[0])):
            if mem[0][j] == int(member_point_loads[i, 0]):
                member = j
        if member is not None:
            frame.members[member].add_point_load(member_point_loads[i, 1].tolist(),
                                                 member_point_loads[i, 2].tolist(),
                                                 member_point_loads[i, 3:].tolist())
        else:
            print("Mmeber not found. For Member Point Load: ", int(member_point_loads[i, 0]))

    # Member Distributed Loads
    md_df = pd.read_excel(file_path, sheet_name='Member_Dist_Loads')
    member_dist_loads = [md_df["Member"].tolist(), md_df["Case"].tolist(), md_df["X1"].tolist(), md_df["X2"].tolist(),
                         md_df["WX1"].tolist(), md_df["WX2"].tolist(), md_df["WY1"].tolist(), md_df["WY2"].tolist(),
                         md_df["WZ1"].tolist(), md_df["WZ2"].tolist()]
    member_dist_loads = np.array(member_dist_loads).T
    for i in range(len(member_dist_loads)):
        member = None
        for j in range(len(mem[0])):
            if mem[0][j] == int(member_dist_loads[i, 0]):
                member = j
        if member is not None:
            frame.members[member].add_dist_load(member_dist_loads[i, 1].tolist(),
                                                member_dist_loads[i, 2:4].tolist(),
                                                member_dist_loads[i, 4:].tolist())
        else:
            print("Mmeber not found. For Member Distributed Load: ", int(member_dist_loads[i, 0]))

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
    inital = None
    try:
        inital = []
        d_inital = mgt_df["Inital d"].tolist()
        b_inital = mgt_df["Inital b"].tolist()
        t_inital = mgt_df["Inital t"].tolist()
        for i in range(len(group_types)):
            inital.append(float(d_inital[i]))
            if group_types[i] == "Angle" or group_types[i] == "RectHSS":
                inital.append(float(b_inital[i]))
            inital.append(float(t_inital[i]))
    except KeyError:
        pass

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

    try:
        method_of_defination = cost_df["Cost Function Name"].tolist()[0]
    except IndexError:
        fail = True
        tk.messagebox.showinfo("Undefined Information", "Please define a valid 'Cost Function Name'.")

    if not fail:
        # noinspection PyUnboundLocalVariable
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
            fail = True
            tk.messagebox.showinfo("Undefined Information", "Please define a valid 'Cost Function Name'.")

    if not fail:
        import_data = [file_path, nodes_df, mem_df, mat_df, s_df, r_df, n_df, mp_df, md_df, mga_df, mgt_df, cost_df]

        # noinspection PyUnboundLocalVariable
        return (group_assignments, group_types, lower_bounds, upper_bounds, inital, cost_function, weight_run,
                reaction_run, internal_forces_run, import_data)
    return None, None, None, None, None, None, None, None, None, None

def frame_3d_opt_temp(file_path: str, opt_type: str) -> None:
    """
    Exports a template Excel file forff the 3D frame global analysis.
    """

    nodes_df = pd.DataFrame([["","","",""]], columns=["Name","X", "Y", "Z"])
    mem_df = pd.DataFrame([["","","","","","","","",""]],
                          columns=["Name","i Node", "j Node", "Material", "Set Cross-Section Properties", "A", "Iy", "Iz",
                                   "J"])
    mat_df = pd.DataFrame([["","","","","",""]], columns=["Name","E", "G", "nu", "rho", "fy"])
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
    if opt_type == "Global":
        mgt_df = pd.DataFrame([["","","","","","",""]],
                          columns=["Group Cross-Section Type", "Min d", "Max d", "Min b", "Max b", "Min t", "Max t"])
    else:
        mgt_df = pd.DataFrame([["", "", "", "", "", "", "","","",""]],
                              columns=["Group Cross-Section Type", "Min d", "Max d", "Min b", "Max b", "Min t",
                                       "Max t","Inital d", "Inital b", "Inital t"])
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

    tk.messagebox.showinfo("File Made", "Your excel tmeplate has been made.")

def frame_2d_global_opt(file_path: str, controller: __main__.Structural_Solver, frame: Frame2D) -> None:
    """
    Optimizes the cross-sections of a 2D frame in the global space.

    :param file_path: File path to the Excel file to import
    :type file_path: String
    :param controller: Object containing the main program.
    :type controller: __main__.Structural_Solver
    :param frame: Frame to be optimized.
    :type frame: Frame2D
    """

    (group_assignments, group_types, lower_bounds, upper_bounds, inital, cost_function, weight_run,
    reaction_run, import_data) = _import_2D_opt_data_from_excel(file_path, frame)

    if group_assignments is not None:

        # noinspection PyBroadException
        # try:
            # noinspection PyUnboundLocalVariable
            frame.global_optimization(group_assignments, group_types, lower_bounds, upper_bounds, cost_function,
                                      weight_run, reaction_run, import_data = import_data)
        # except Exception as e:
        #     tk.messagebox.showinfo("Error", "An error ha occurred while running the analysis: {e}")

def frame_2d_local_opt(file_path: str, controller: __main__.Structural_Solver, frame: Frame2D) -> None:
    """
    Optimizes the cross-sections of a 2D frame in the local space.

    :param file_path: File path to the Excel file to import
    :type file_path: String
    :param controller: Object containing the main program.
    :type controller: __main__.Structural_Solver
    :param frame: Frame to be optimized.
    :type frame: Frame2D
    """

    (group_assignments, group_types, lower_bounds, upper_bounds, inital, cost_function, weight_run,
    reaction_run, import_data) = _import_2D_opt_data_from_excel(file_path, frame)

    if group_assignments is not None:

        # noinspection PyBroadException
        try:
            # noinspection PyUnboundLocalVariable
            frame.local_optimization(group_assignments, group_types, lower_bounds, upper_bounds, inital, cost_function,
                                      weight_run, reaction_run, import_data = import_data)
        except Exception as e:
            tk.messagebox.showinfo("Error", "An error ha occurred while running the analysis: {e}")

def results_frame_2d_opt(cost: float, results: list, import_data: list, optimization_type: str) -> None:
    """
    Writes the results of the 3D frame optimization of the global space to the Excel file.

    :param cost: Best cost found.
    :type cost: float
    :param results: The Cross-sections that resulted in the best cost.
    :type results: list.
    :param import_data: List of information that was imported.
    :type import_data: list
    :param optimization_type: Type of optimization that was performed.
    :type optimization_type: str
    """

    [file_path, nodes_df, mem_df, mat_df, s_df, r_df, n_df, mp_df, md_df, mga_df, mgt_df, cost_df] = import_data
    df1 = pd.DataFrame(results,
                       columns=["Member Index", "Cross-Section _type", "d", "b", "t", "A", "Iy", "Iz", "J"])

    if optimization_type == "Local":
        df2 = pd.DataFrame({'Cost': [cost], 'Optimization Type': [optimization_type]})
    else:
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

    tk.messagebox.showinfo("Analysis Complete", "Analysis Complete")

def _import_2D_opt_data_from_excel(file_path: str, frame: Frame2D):
    """
    Imports information form the specified Excel file

    :param file_path: File path to the Excel file to import
    :type file_path: String
    :param frame: Frame to be optimized.
    :type frame: Frame2D
    :return: group_assignments, group_types, lower_bounds, upper_bounds, inital, cost_function, weight_run,
                reaction_run, internal_forces_run, import_data
    """

    fail = False

    # Nodes
    nodes_df = pd.read_excel(file_path, sheet_name='Nodes')
    nodes = [nodes_df["Name"].tolist(), nodes_df["X"].tolist(), nodes_df["Y"].tolist()]
    for i in range(len(nodes[0])): frame.nodes.append(Node(x=nodes[1][i], y=nodes[2][i], z=0))

    # Materials
    mat_df = pd.read_excel(file_path, sheet_name='Materials')
    mat = [mat_df["Name"].tolist(), mat_df["E"].tolist(), mat_df["G"].tolist(), mat_df["nu"].tolist(),
           mat_df["rho"].tolist(), mat_df["fy"].tolist()]
    for i in range(len(mat[0])): frame.materials.append(Material(mat[1][i], mat[2][i], mat[3][i], mat[4][i], mat[5][i]))

    # Members
    mem_df = pd.read_excel(file_path, sheet_name='Members')
    mem = [mem_df["Name"].tolist(), mem_df["i Node"].tolist(), mem_df["j Node"].tolist(), mem_df["Material"].tolist(),
           mem_df["Set Cross-Section Properties"].tolist(), mem_df["A"].tolist(), mem_df["I"].tolist()]
    for i in range(len(mem[0])):
        nodei = None
        nodej = None
        material = None
        for j in range(len(nodes[0])):
            if nodes[0][j] == mem[1][i]:
                nodei = j
            elif nodes[0][j] == mem[2][i]:
                nodej = j
        for j in range(len(mat[0])):
            if mat[0][j] == mem[3][i]:
                material = j
        if nodei is not None and nodej is not None and material is not None:
            frame.members.append(Member(nodei, nodej, material, mem[4][i],
                                        [mem[5][i], mem[6][i], 0, 0], frame))
        else:
            print("Node or Material not found. For member: ", mem[0][i])

    # Supports
    s_df = pd.read_excel(file_path, sheet_name='Supports')
    supports = [s_df["Node"].tolist(), s_df["DX"].tolist(), s_df["DY"].tolist(), s_df["R"].tolist()]
    supports = np.array(supports).T
    for i in range(len(supports)):
        node = None
        for j in range(len(nodes[0])):
            if nodes[0][j] == supports[i, 0]:
                node = j
        if node is not None:
            frame.nodes[node].set_support(supports[i, 1:].tolist())
        else:
            print("Support not found. For node: ", supports[i, 0])

    # Releases
    r_df = pd.read_excel(file_path, sheet_name='Releases')
    releases = [r_df["Member"].tolist(), r_df["i DX"].tolist(), r_df["i DY"].tolist(),
                r_df["i R"].tolist(), r_df["j DX"].tolist(), r_df["j DY"].tolist(), r_df["j R"].tolist()]
    releases = np.array(releases).T
    for i in range(len(releases)):
        member = None
        for j in range(len(mem[0])):
            if mem[0][j] == releases[i, 0]:
                member = j
        if member is not None:
            frame.members[member].set_releces(releases[i, 1:].tolist())
        else:
            print("Mmeber not found. For Release: ", releases[i, 0])

    # Node Loads
    n_df = pd.read_excel(file_path, sheet_name='Node_Loads')
    node_loads = [n_df["Node"].tolist(), n_df["Case"].tolist(), n_df["PX"].tolist(), n_df["PY"].tolist(),
                  n_df["M"].tolist()]
    node_loads = np.array(node_loads).T
    for i in range(len(node_loads)):
        node = None
        for j in range(len(nodes[0])):
            if nodes[0][j] == node_loads[i, 0]:
                node = j
        if node is not None:
            frame.nodes[node].add_load(node_loads[i, 1].tolist(), node_loads[i, 2:].tolist())
        else:
            print("Support not found. For Node Load: ", node_loads[i, 0])

    # Member Point Loads
    mp_df = pd.read_excel(file_path, sheet_name='Member_Point_Loads')
    member_point_loads = [mp_df["Member"].tolist(), mp_df["Case"].tolist(), mp_df["X"].tolist(), mp_df["PX"].tolist(),
                          mp_df["PY"].tolist(), mp_df["MX"].tolist()]
    member_point_loads = np.array(member_point_loads).T
    for i in range(len(member_point_loads)):
        member = None
        for j in range(len(mem[0])):
            if mem[0][j] == int(member_point_loads[i, 0]):
                member = j
        if member is not None:
            frame.members[member].add_point_load(member_point_loads[i, 1].tolist(),
                                                 member_point_loads[i, 2].tolist(),
                                                 member_point_loads[i, 3:].tolist())
        else:
            print("Mmeber not found. For Member Point Load: ", int(member_point_loads[i, 0]))

    # Member Distributed Loads
    md_df = pd.read_excel(file_path, sheet_name='Member_Dist_Loads')
    member_dist_loads = [md_df["Member"].tolist(), md_df["Case"].tolist(), md_df["X1"].tolist(), md_df["X2"].tolist(),
                         md_df["WX1"].tolist(), md_df["WX2"].tolist(), md_df["WY1"].tolist(), md_df["WY2"].tolist()]
    member_dist_loads = np.array(member_dist_loads).T
    for i in range(len(member_dist_loads)):
        member = None
        for j in range(len(mem[0])):
            if mem[0][j] == int(member_dist_loads[i, 0]):
                member = j
        if member is not None:
            frame.members[member].add_dist_load(member_dist_loads[i, 1].tolist(),
                                                member_dist_loads[i, 2:4].tolist(),
                                                member_dist_loads[i, 4:].tolist())
        else:
            print("Mmeber not found. For Member Distributed Load: ", int(member_dist_loads[i, 0]))

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
    inital = None
    try:
        inital = []
        d_inital = mgt_df["Inital d"].tolist()
        b_inital = mgt_df["Inital b"].tolist()
        t_inital = mgt_df["Inital t"].tolist()
        for i in range(len(group_types)):
            inital.append(float(d_inital[i]))
            if group_types[i] == "Angle" or group_types[i] == "RectHSS":
                inital.append(float(b_inital[i]))
            inital.append(float(t_inital[i]))
    except KeyError:
        pass

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

    try:
        method_of_defination = cost_df["Cost Function Name"].tolist()[0]
    except IndexError:
        fail = True
        tk.messagebox.showinfo("Undefined Information", "Please define a valid 'Cost Function Name'.")

    if not fail:
        # noinspection PyUnboundLocalVariable
        if method_of_defination == "Other":
            cost_function = cost_df["Function"].tolist()[0]
            weight_run = cost_df["Weight Run"].tolist()[0]
            reaction_run = cost_df["Reaction Run"].tolist()[0]
        elif method_of_defination == "2027_Steel_Bridge":
            cost_function = "max(max(DY)) + sum(Weight)"
            weight_run = True
            reaction_run = False
        else:
            fail = True
            tk.messagebox.showinfo("Undefined Information", "Please define a valid 'Cost Function Name'.")

    if not fail:
        import_data = [file_path, nodes_df, mem_df, mat_df, s_df, r_df, n_df, mp_df, md_df, mga_df, mgt_df, cost_df]

        # noinspection PyUnboundLocalVariable
        return (group_assignments, group_types, lower_bounds, upper_bounds, inital, cost_function, weight_run,
                reaction_run, import_data)
    return None, None, None, None, None, None, None, None, None, None

def frame_2d_opt_temp(file_path: str, opt_type: str) -> None:
    """
    Exports a template Excel file forff the 3D frame global analysis.
    """

    nodes_df = pd.DataFrame([["","",""]], columns=["Name","X", "Y"])
    mem_df = pd.DataFrame([["","","","","","",""]],
                          columns=["Name","i Node", "j Node", "Material", "Set Cross-Section Properties", "A", "I"])
    mat_df = pd.DataFrame([["","","","","",""]], columns=["Name","E", "G", "nu", "rho", "fy"])
    s_df = pd.DataFrame([["","","",""]], columns=["Node", "DX", "DY", "R"])
    r_df = pd.DataFrame([["","","","","","",""]],
                        columns=["Member","i DX", "i DY", "i R", "j DX", "j DY", "j R"])
    n_df = pd.DataFrame([["","","","",""]], columns=["Node", "Case", "PX", "PY", "M"])
    mp_df = pd.DataFrame([["","","","","",""]],
                         columns=["Member", "Case", "X", "PX", "PY", "MX"])
    md_df = pd.DataFrame([["","","","","","","",""]],
                         columns=["Member", "Case", "X1", "X2", "WX1", "WX2", "WY1", "WY2"])
    mga_df = pd.DataFrame([[""]], columns=["Group_Number"])
    if opt_type == "Global":
        mgt_df = pd.DataFrame([["","","","","","",""]],
                          columns=["Group Cross-Section Type", "Min d", "Max d", "Min b", "Max b", "Min t", "Max t"])
    else:
        mgt_df = pd.DataFrame([["", "", "", "", "", "", "","","",""]],
                              columns=["Group Cross-Section Type", "Min d", "Max d", "Min b", "Max b", "Min t",
                                       "Max t","Inital d", "Inital b", "Inital t"])
    cost_df = pd.DataFrame([["","","",""]],
                           columns=["Cost Function Name", "Function", "Weight Run", "Reaction Run"])

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

    tk.messagebox.showinfo("File Made", "Your excel tmeplate has been made.")
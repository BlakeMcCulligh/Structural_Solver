import numpy as np
import pandas as pd


def saveTrussCrossSectionOptimizationExcel(filePath, nodes, members, loadCasses, supports, areas, deflections, forces, volume):

    Nodes_df = convertNodesToDF(nodes, deflections)
    Members_df = convertMembersToDF(members, areas, forces)
    Supports_df = convertSupportsToDF(supports)
    Volume_df = pd.DataFrame({"Volume": [volume]})
    LoadCases_df = convertLoadCassesToDF(loadCasses)
    writeToExcel(filePath, Nodes_df, Members_df, Supports_df, LoadCases_df, Volume_df)

def saveTrussTopologyOptimizationExcel(filePath, nodes, members, loadCasses, supports, areas, deflections, forces, volume):

    Nodes_df = convertNodesToDF(nodes, deflections)
    Members_df = convertMembersToDF(members, areas, forces, removeSmallMembers=True)
    Supports_df = convertSupportsToDF(supports)
    Volume_df = pd.DataFrame({"Volume": [volume]})
    LoadCases_df = convertLoadCassesToDF(loadCasses)
    writeToExcel(filePath, Nodes_df, Members_df, Supports_df, LoadCases_df, Volume_df)


def convertNodesToDF(nodes, deflections):
    nodes = np.array(nodes).T.tolist()
    deflections = deflections.reshape(-1, len(nodes[0])).tolist()

    Nodes_df = {"X": nodes[0],
                "Y": nodes[1]}
    for i, u in enumerate(deflections):
        Nodes_df = Nodes_df | {f"Delection{i + 1}": u}

    Nodes_df = pd.DataFrame(Nodes_df)
    return Nodes_df

def convertMembersToDF(members, areas, forces, removeSmallMembers=False):
    maxA = max(areas)
    if removeSmallMembers:
        areas = areas.tolist()
        members = members.tolist()
        forces = forces.T.tolist()
        for i in reversed(range(len(members))):
            print("I:", i)
            print("areas:", areas)
            if areas[i] < maxA * 0.001:
                areas.pop(i)
                members.pop(i)
                forces.pop(i)
        forces = np.array(forces).T.tolist()

    members = np.array(members).T.tolist()

    Members_df = {"Node1": members[0],
                 "Node2": members[1],
                 "Areas": areas}
    for i, q in enumerate(forces):
        Members_df = Members_df | {f"Force{i + 1}": q}
    Members_df = pd.DataFrame(Members_df)
    return Members_df

def convertSupportsToDF(supports):
    supports = np.array(supports).T.tolist()

    Supports_df = {"X": supports[0],
                   "Y": supports[1],
                   "Support x": supports[2],
                   "Support y": supports[3]}
    Supports_df = pd.DataFrame(Supports_df)

    return Supports_df

def convertLoadCassesToDF(loadCasses):
    # TODO only works if all load cases have the same number of nodes
    numLoadCases = len(loadCasses)
    LoadCases_df = {"Number Load Casses": [numLoadCases] * len(loadCasses[0])}
    print("loadCasses:", loadCasses)
    for i, loadCase in enumerate(loadCasses):
        LoadCase = np.array(loadCase).T.tolist()
        load = {f"LoadCase{i + 1} x": LoadCase[0],
                f"LoadCase{i + 1} y": LoadCase[1],
                f"LoadCase{i + 1} fx": LoadCase[2],
                f"LoadCase{i + 1} fy": LoadCase[3]}

        LoadCases_df = LoadCases_df | load
    LoadCases_df = pd.DataFrame(LoadCases_df)

    return LoadCases_df

def writeToExcel(filePath, Nodes_df, Members_df, Supports_df, LoadCases_df, Volume_df):
    with pd.ExcelWriter(filePath) as writer:
        Nodes_df.to_excel(writer, sheet_name="Nodes", index=False)
        Members_df.to_excel(writer, sheet_name="Members", index=False)
        Supports_df.to_excel(writer, sheet_name="Supports", index=False)
        LoadCases_df.to_excel(writer, sheet_name="LoadCases", index=False)
        Volume_df.to_excel(writer, sheet_name="Volume", index=False)


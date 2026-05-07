import tkinter as tk
from tkinter import filedialog

import numpy as np

from Frame3DGUI.inputData import data
from Frame3DGUI.results import Results


def openFrame(window, file_path = None):
    if file_path is None:
        fileTypes = [("Struct Frame files", "*.structframe")]
        file_path = select_file_gui(window, fileTypes)
        window.filePath = file_path.replace(".structframe","")

    if file_path is not None:
        with open(file_path, "r") as f:
            lines = f.readlines()
            f.close()

        # clearing current file
        window.data = data()
        window.Results = Results()
        window.printNodes = np.empty((0, 3))
        window.printLines = np.empty((0, 2, 3))
        window.surfaceTri = np.empty((0, 3, 3))
        window.solidTri = np.empty((0, 3, 3))
        for i in range(len(window.tables)):
            try:
                window.tables[i].delete(*window.tables[i].get_children())
            except AttributeError:
                print("failed to delete table data")

        # nodes
        nodes = []
        numNodes = int(lines[0])
        for i in range(numNodes):
            nodeLine = lines[i+1].replace(" ","")
            n = []
            for x in nodeLine.split(','):
                try:
                    n.append(float(x))
                except ValueError:
                    pass
            nodes.append(n)
        for i in range(len(nodes)):
            window.data.addnodes(window, nodes[i], True, True)
        nodes = np.array(nodes)
        index = numNodes+1

        # materials
        materials = []
        numMaterials = int(lines[index])
        for i in range(numMaterials):
            materialLine = lines[i+index+1].replace(" ","")
            n = []
            for x in materialLine.split(','):
                try:
                    n.append(float(x))
                except ValueError:
                    pass
            materials.append(n)
        for i in range(len(materials)):
            window.data.addmaterials(window, materials[i], True, True)
        index += numMaterials + 1

        # members
        members = []
        numMembers = int(lines[index])
        for i in range(numMembers):
            memberLine = lines[i+index+1].replace(" ","")
            n = []
            for x in memberLine.split(','):
                try:
                    if x == "True" or x == "False":
                        n.append(bool(x))
                    else:
                        n.append(float(x))
                except ValueError:
                    pass
            members.append(n)
        for i in range(len(members)):
            window.data.addmembers(window, members[i], True, True)
        index += numMembers + 1
        members = np.array(members)

        #Supports
        supports = []
        numSupports = int(lines[index])
        for i in range(numSupports):
            supportLine = lines[i+index+1].replace(" ","")
            n = []
            for x in supportLine.split(','):
                try:
                    if x == "True" or x == "False":
                        n.append(bool(x))
                    else:
                        n.append(float(x))
                except ValueError:
                    pass
            supports.append(n)
        for i in range(len(supports)):
            window.data.addsupports(window, supports[i], True, True)
        index += numSupports + 1

        # releases
        releases = []
        numReleases = int(lines[index])
        for i in range(numReleases):
            releaseLine = lines[i+index+1].replace(" ","")
            n = []
            for x in releaseLine.split(','):
                try:
                    if x == "True" or x == "False":
                        n.append(bool(x))
                    else:
                        n.append(float(x))
                except ValueError:
                    pass
            releases.append(n)
        for i in range(len(releases)):
            window.data.addreleases(window, releases[i], True, True)
        index += numReleases + 1

        # nodeLoad
        nodeLoads = []
        numNodeLoads = int(lines[index])
        for i in range(numNodeLoads):
            nodeLoadLine = lines[i+index+1].replace(" ","")
            n = []
            for x in nodeLoadLine.split(','):
                try:
                    if x == "True" or x == "False":
                        n.append(bool(x))
                    else:
                        n.append(float(x))
                except ValueError:
                    pass
            nodeLoads.append(n)
        for i in range(len(nodeLoads)):
            window.data.addnodeLoads(window, nodeLoads[i], True, True)
        index += numNodeLoads + 1

        # memberPointLoad
        memberPointLoads = []
        numMemberPointLoads = int(lines[index])
        for i in range(numMemberPointLoads):
            memberPointLine = lines[i+index+1].replace(" ","")
            n = []
            for x in memberPointLine.split(','):
                try:
                    if x == "True" or x == "False":
                        n.append(bool(x))
                    else:
                        n.append(float(x))
                except ValueError:
                    pass
            memberPointLoads.append(n)
        for i in range(len(memberPointLoads)):
            window.data.addmemberPointLoads(window, memberPointLoads[i], True, True)
        index += numMemberPointLoads + 1

        # memberDistLoad
        memberDistLoads = []
        numMemberDistLoads = int(lines[index])
        for i in range(numMemberDistLoads):
            memberDistLine = lines[i+index+1].replace(" ","")
            n = []
            for x in memberDistLine.split(','):
                try:
                    if x == "True" or x == "False":
                        n.append(bool(x))
                    else:
                        n.append(float(x))
                except ValueError:
                    pass
            memberDistLoads.append(n)
        for i in range(len(memberDistLoads)):
            window.data.addmemberDistLoads(window, memberDistLoads[i], True, True)

def select_file_gui(window, file_Types):
    """
    Opens a file explorer dialog using Tkinter and returns the selected file path.
    """

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

def openResults(window, filePath):
    if filePath is not None:
        filePath = filePath + ".structresult"

        try:
            with open(filePath, "r") as f:
                lines = f.readlines()
                f.close()

            numCasses = int(lines[0])
            numNodes = int(lines[1])

            d = [[],[],[],[],[],[]]
            for i in range(numCasses):
                deflection_sub = [[],[],[],[],[],[]]
                for j in range(numNodes):
                    deflectionLine = lines[i*numNodes+j+2].replace(" ","")
                    k = 0
                    for x in deflectionLine.split(','):
                        try:
                            deflection_sub[k].append(float(x))
                            k += 1
                        except ValueError:
                            if x != ",": print("Error Reading Deflections: ", x)
                for k in range(6): d[k].append(deflection_sub[k])
            window.Results.addNodalDeflections(d[0],d[1],d[2],d[3],d[4],d[5])
            index = numCasses * numNodes + 2

            numMembers = int(lines[index])

            weight = []
            for i in range(numMembers):
                weight.append(float(lines[index+i+1].replace(" ","")))
            window.Results.addWeight(weight)
            index += numMembers + 2

            numCasses = int(lines[index])
            numNodes = int(lines[index + 1])

            reactions = []
            for i in range(numCasses):
                subReactions = []
                for j in range(numNodes):
                    reactionLine = lines[index + 2 + i * numNodes + j].replace(" ","")
                    n = []
                    for x in reactionLine.split(','):
                        try:
                            n.append(float(x))
                        except ValueError:
                            pass
                    subReactions.append(n)
                reactions.append(subReactions)
            window.Results.addReactions(reactions)
            index += numCasses * numNodes + 2

            F = []
            for i in range(6):
                fLine = lines[index + i].replace(" ","")
                f = []
                for x in fLine.split(','):
                    try: f.append(float(x))
                    except ValueError: pass
                F.append(f)

            F_formated = [[F[0],F[1],F[2]],[F[3],F[4],F[5]]]
            window.Results.addInternalForces(F_formated)

        except FileNotFoundError:
            pass


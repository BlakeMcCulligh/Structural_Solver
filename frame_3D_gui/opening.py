"""
Handles the opening of 3D frames, and results.
"""

from tkinter import filedialog
import numpy as np

from frame_3D_gui.data import Data
from frame_3D_gui.results import Results

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

def open_frame(window, file_path = None):
    """
    Opens a .structframe file eather from the provided file path or the selected one and reads in all the frame data.

    :param window: Object storing the main window.
    :param file_path: Potential file Path to get the frame from.
    """

    if file_path is None:
        fileTypes = [("Struct Frame files", "*.structframe")]
        file_path = _select_file_gui(window, fileTypes)
        window.FilePath = file_path.replace(".structframe", "")

    if file_path is not None:
        with open(file_path, "r") as f:
            lines = f.readlines()
            f.close()

        # clearing current file
        window.Data = Data()
        window.Results = Results()
        window.PrintNodes = np.empty((0, 3))
        window.PrintLines = np.empty((0, 2, 3))
        window.PrintSurfaceTri = np.empty((0, 3, 3))
        window.PrintSolidTri = np.empty((0, 3, 3))
        for i in range(len(window.Tables)):
            try:
                window.Tables[i].delete(*window.Tables[i].get_children())
            except AttributeError:
                print("failed to delete table Data")

        # Nodes
        nodes = []
        num_nodes = int(lines[0])
        for i in range(num_nodes):
            node_line = lines[i+1].replace(" ","")
            n = []
            for x in node_line.split(','):
                try:
                    n.append(float(x))
                except ValueError:
                    pass
            nodes.append(n)
        for i in range(len(nodes)):
            window.Data.AddNodes(window, nodes[i], True, True)
        nodes = np.array(nodes)
        index = num_nodes+1

        # Materials
        materials = []
        num_materials = int(lines[index])
        for i in range(num_materials):
            material_line = lines[i+index+1].replace(" ","")
            n = []
            for x in material_line.split(','):
                try:
                    n.append(float(x))
                except ValueError:
                    pass
            materials.append(n)
        for i in range(len(materials)):
            window.Data.AddMaterials(window, materials[i], True, True)
        index += num_materials + 1

        # Members
        members = []
        num_members = int(lines[index])
        for i in range(num_members):
            member_line = lines[i+index+1].replace(" ","")
            n = []
            for x in member_line.split(','):
                try:
                    if x == "True" or x == "False":
                        n.append(bool(x))
                    else:
                        n.append(float(x))
                except ValueError:
                    pass
            members.append(n)
        for i in range(len(members)):
            window.Data.AddMembers(window, members[i], True, True)
        index += num_members + 1
        members = np.array(members)

        #Supports
        supports = []
        num_supports = int(lines[index])
        for i in range(num_supports):
            support_line = lines[i+index+1].replace(" ","")
            n = []
            for x in support_line.split(','):
                try:
                    if x == "True" or x == "False":
                        n.append(bool(x))
                    else:
                        n.append(float(x))
                except ValueError:
                    pass
            supports.append(n)
        for i in range(len(supports)):
            window.Data.AddSupports(window, supports[i], True, True)
        index += num_supports + 1

        # Releases
        releases = []
        num_releases = int(lines[index])
        for i in range(num_releases):
            release_line = lines[i+index+1].replace(" ","")
            n = []
            for x in release_line.split(','):
                try:
                    if x == "True" or x == "False":
                        n.append(bool(x))
                    else:
                        n.append(float(x))
                except ValueError:
                    pass
            releases.append(n)
        for i in range(len(releases)):
            window.Data.AddReleases(window, releases[i], True, True)
        index += num_releases + 1

        # NodeLoad
        node_loads = []
        num_node_loads = int(lines[index])
        for i in range(num_node_loads):
            node_load_line = lines[i+index+1].replace(" ","")
            n = []
            for x in node_load_line.split(','):
                try:
                    if x == "True" or x == "False":
                        n.append(bool(x))
                    else:
                        n.append(float(x))
                except ValueError:
                    pass
            node_loads.append(n)
        for i in range(len(node_loads)):
            window.Data.AddNodeLoads(window, node_loads[i], True, True)
        index += num_node_loads + 1

        # MemberPointLoad
        member_point_loads = []
        num_member_point_loads = int(lines[index])
        for i in range(num_member_point_loads):
            member_point_line = lines[i+index+1].replace(" ","")
            n = []
            for x in member_point_line.split(','):
                try:
                    if x == "True" or x == "False":
                        n.append(bool(x))
                    else:
                        n.append(float(x))
                except ValueError:
                    pass
            member_point_loads.append(n)
        for i in range(len(member_point_loads)):
            window.Data.AddMemberPointLoads(window, member_point_loads[i], True, True)
        index += num_member_point_loads + 1

        # MemberDistLoad
        member_dist_loads = []
        num_member_dist_loads = int(lines[index])
        for i in range(num_member_dist_loads):
            member_dist_line = lines[i+index+1].replace(" ","")
            n = []
            for x in member_dist_line.split(','):
                try:
                    if x == "True" or x == "False":
                        n.append(bool(x))
                    else:
                        n.append(float(x))
                except ValueError:
                    pass
            member_dist_loads.append(n)
        for i in range(len(member_dist_loads)):
            window.Data.AddMemberDistLoads(window, member_dist_loads[i], True, True)

def _select_file_gui(window, file_types):
    """
    Opens a file explorer dialog using Tkinter and returns the selected file path.

    :param window:  Object storing the main window.
    :param file_types: What type of file is being opened.
    :return: Selected file path or None.
    """

    # Open the file dialog
    file_path = filedialog.askopenfilename(title="Select a file", filetypes= file_types)

    if file_path:return file_path
    else: return None

def open_results(window, file_path):
    """
    Opens a .structresult file from the provided file path and reads in all the structural results.

    :param window:  Object storing the main window.
    :param file_path: File path to get the structural results from.
    """

    if file_path is not None:
        file_path = file_path + ".structresult"

        try:
            with open(file_path, "r") as f:
                lines = f.readlines()
                f.close()

            num_cases = int(lines[0])
            num_nodes = int(lines[1])

            d = [[],[],[],[],[],[]]
            for i in range(num_cases):
                deflection_sub = [[],[],[],[],[],[]]
                for j in range(num_nodes):
                    deflection_line = lines[i*num_nodes+j+2].replace(" ","")
                    k = 0
                    for x in deflection_line.split(','):
                        try:
                            deflection_sub[k].append(float(x))
                            k += 1
                        except ValueError:
                            if x != ",": print("Error Reading Deflections: ", x)
                for k in range(6): d[k].append(deflection_sub[k])
            window.Results.AddNodalDeflections(d[0], d[1], d[2], d[3], d[4], d[5])
            index = num_cases * num_nodes + 2

            num_members = int(lines[index])

            weight = []
            for i in range(num_members):
                weight.append(float(lines[index+i+1].replace(" ","")))
            window.Results.add_weight(weight)
            index += num_members + 2

            num_cases = int(lines[index])
            num_nodes = int(lines[index + 1])

            reactions = []
            for i in range(num_cases):
                sub_reactions = []
                for j in range(num_nodes):
                    reaction_line = lines[index + 2 + i * num_nodes + j].replace(" ","")
                    n = []
                    for x in reaction_line.split(','):
                        try:
                            n.append(float(x))
                        except ValueError:
                            pass
                    sub_reactions.append(n)
                reactions.append(sub_reactions)
            window.Results.AddReactions(reactions)
            index += num_cases * num_nodes + 2

            F = []
            for i in range(6):
                fLine = lines[index + i].replace(" ","")
                f = []
                for x in fLine.split(','):
                    try: f.append(float(x))
                    except ValueError: pass
                F.append(f)

            F_formated = [[F[0],F[1],F[2]],[F[3],F[4],F[5]]]
            window.Results.AddInternalForces(F_formated)

        except FileNotFoundError:
            pass
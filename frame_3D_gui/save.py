"""
Handels saving of 3D frames and results
"""

from tkinter import filedialog

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

def save_frame(data, file_path = None):
    """
    Saves frames to .structframe file using the provided file path or props the user to chose one.

    :param data: Info that defines the frame.
    :param file_path: File path to save the file to.
    :return: file path.
    """

    if file_path is None:
        file_path = _get_file_save_path_frame()
    else:
        file_path = file_path + ".structframe"

    if file_path is not None:
        with open(file_path, "w") as f:

            f.write(f"{len(data.Nodes[0])}\n")
            for i in range(len(data.Nodes[0])):
                f.write(f"{data.Nodes[0][i]},{data.Nodes[1][i]},{data.Nodes[2][i]},\n")

            f.write(f"{len(data.Materials[0])}\n")
            for i in range(len(data.Materials[0])):
                string = ""
                for j in range(len(data.Materials)): string = string + f"{data.Materials[j][i]},"
                f.write(string)
                f.write(f"\n")


            f.write(f"{len(data.Members[0])}\n")
            for i in range(len(data.Members[0])):
                string = ""
                for j in range(len(data.Members)): string = string + f"{data.Members[j][i]},"
                f.write(string)
                f.write(f"\n")

            f.write(f"{len(data.Supports[0])}\n")
            for i in range(len(data.Supports[0])):
                string = ""
                for j in range(len(data.Supports)): string = string + f"{data.Supports[j][i]},"
                f.write(string)
                f.write(f"\n")

            f.write(f"{len(data.Releases[0])}\n")
            for i in range(len(data.Releases[0])):
                string = ""
                for j in range(len(data.Releases)): string = string + f"{data.Releases[j][i]},"
                f.write(string)
                f.write(f"\n")

            f.write(f"{len(data.NodeLoad[0])}\n")
            for i in range(len(data.NodeLoad[0])):
                string = ""
                for j in range(len(data.NodeLoad)): string = string + f"{data.NodeLoad[j][i]},"
                f.write(string)
                f.write(f"\n")

            f.write(f"{len(data.MemberPointLoad[0])}\n")
            for i in range(len(data.MemberPointLoad[0])):
                string = ""
                for j in range(len(data.MemberPointLoad)): string = string + f"{data.MemberPointLoad[j][i]},"
                f.write(string)
                f.write(f"\n")

            f.write(f"{len(data.MemberDistLoad[0])}\n")
            for i in range(len(data.MemberDistLoad[0])):
                string = ""
                for j in range(len(data.MemberDistLoad)): string = string + f"{data.MemberDistLoad[j][i]},"
                f.write(string)
                f.write(f"\n")

    return file_path.replace(".structframe", "")

def _get_file_save_path_frame():
    """
    Gets the file path to save the frame to.

    :return: File path or None.
    """

    file_path = filedialog.asksaveasfilename(defaultextension=".structframe",
                                             filetypes=[("Struct Frame files", "*.structframe")])
    if file_path: return file_path
    return None

def save_results(results, file_path):
    """
    Saves results to .structresult file using the provided filepath.

    :param results: Results to be saved.
    :param file_path: File path to save the file to.
    """

    file_path = file_path + ".structresult"

    if file_path is not None:
        with open(file_path, "w") as f:

            f.write(f"{len(results.NodalDeflections)}\n")
            f.write(f"{len(results.NodalDeflections.DX)}\n")
            for i in range(len(results.NodalDeflections)):
                for j in range(len(results.NodalDeflections.DX)):
                    f.write(f"{results.NodalDeflections[i].DX[j]},{results.NodalDeflections[i].DY[j]},"
                            f"{results.NodalDeflections[i].DZ[j]},{results.NodalDeflections[i].RX[j]},"
                            f"{results.NodalDeflections[i].RY[j]},{results.NodalDeflections[i].RZ[j]},\n")

            f.write(f"{len(results.Weight)}\n")
            for i in range(len(results.Weight)):
                f.write(f"{results.Weight[i]}\n")

            f.write(f"{results.OverallWeight}\n")

            f.write(f"{len(results.Reactions)}\n")
            f.write(f"{len(results.Reactions.RX)}\n")
            for i in range(len(results.Reactions)):
                for j in range(len(results.Reactions.RX)):
                    f.write(f"{results.Reactions[i].RX[j]},{results.Reactions[i].RY[j]},{results.Reactions[i].RZ[j]},"
                            f"{results.Reactions[i].MX[j]},{results.Reactions[i].MY[j]},{results.Reactions[i].MZ[j]}\n")

            f.write(f"{results.MaxInternalForces[0][0][0]},{results.MaxInternalForces[0][0][1]},\n")
            f.write(f"{results.MaxInternalForces[0][1][0]},{results.MaxInternalForces[0][1][1]},\n")
            f.write(f"{results.MaxInternalForces[0][2][0]},{results.MaxInternalForces[0][2][1]},\n")
            f.write(f"{results.MaxInternalForces[1][0][0]},{results.MaxInternalForces[1][0][1]},\n")
            f.write(f"{results.MaxInternalForces[1][1][0]},{results.MaxInternalForces[1][1][1]},\n")
            f.write(f"{results.MaxInternalForces[1][2][0]},{results.MaxInternalForces[1][2][1]},\n")
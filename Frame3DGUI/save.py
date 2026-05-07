from tkinter import filedialog
import tkinter as tk

def saveFrame(data, filePath = None):

    if filePath is None:
        filePath = getFileSavePathFrame()
    else:
        filePath = filePath + ".structframe"

    if filePath is not None:
        with open(filePath, "w") as f:

            f.write(f"{len(data.nodes[0])}\n")
            for i in range(len(data.nodes[0])):
                f.write(f"{data.nodes[0][i]},{data.nodes[1][i]},{data.nodes[2][i]},\n")

            f.write(f"{len(data.materials[0])}\n")
            for i in range(len(data.materials[0])):
                string = ""
                for j in range(len(data.materials)): string = string + f"{data.materials[j][i]},"
                f.write(string)
                f.write(f"\n")


            f.write(f"{len(data.members[0])}\n")
            for i in range(len(data.members[0])):
                string = ""
                for j in range(len(data.members)): string = string + f"{data.members[j][i]},"
                f.write(string)
                f.write(f"\n")

            f.write(f"{len(data.supports[0])}\n")
            for i in range(len(data.supports[0])):
                string = ""
                for j in range(len(data.supports)): string = string + f"{data.supports[j][i]},"
                f.write(string)
                f.write(f"\n")

            f.write(f"{len(data.releases[0])}\n")
            for i in range(len(data.releases[0])):
                string = ""
                for j in range(len(data.releases)): string = string + f"{data.releases[j][i]},"
                f.write(string)
                f.write(f"\n")

            f.write(f"{len(data.nodeLoad[0])}\n")
            for i in range(len(data.nodeLoad[0])):
                string = ""
                for j in range(len(data.nodeLoad)): string = string + f"{data.nodeLoad[j][i]},"
                f.write(string)
                f.write(f"\n")

            f.write(f"{len(data.memberPointLoad[0])}\n")
            for i in range(len(data.memberPointLoad[0])):
                string = ""
                for j in range(len(data.memberPointLoad)): string = string + f"{data.memberPointLoad[j][i]},"
                f.write(string)
                f.write(f"\n")

            f.write(f"{len(data.memberDistLoad[0])}\n")
            for i in range(len(data.memberDistLoad[0])):
                string = ""
                for j in range(len(data.memberDistLoad)): string = string + f"{data.memberDistLoad[j][i]},"
                f.write(string)
                f.write(f"\n")

    return filePath.replace(".structframe","")

def getFileSavePathFrame():
    file_path = filedialog.asksaveasfilename(defaultextension=".structframe",
                                             filetypes=[("Struct Frame files", "*.structframe")])
    if file_path:return file_path
    return None

def saveResults(results, filepath):

    filePath = filepath + ".structresult"

    if filePath is not None:
        with open(filePath, "w") as f:

            f.write(f"{len(results.nodalDeflections)}\n")
            f.write(f"{len(results.nodalDeflections.DX)}\n")
            for i in range(len(results.nodalDeflections)):
                for j in range(len(results.nodalDeflections.DX)):
                    f.write(f"{results.nodalDeflections[i].DX[j]},{results.nodalDeflections[i].DY[j]},"
                            f"{results.nodalDeflections[i].DZ[j]},{results.nodalDeflections[i].RX[j]},"
                            f"{results.nodalDeflections[i].RY[j]},{results.nodalDeflections[i].RZ[j]},\n")

            f.write(f"{len(results.weight)}\n")
            for i in range(len(results.weight)):
                f.write(f"{results.weight[i]}\n")

            f.write(f"{results.overallWeight}\n")

            f.write(f"{len(results.reactions)}\n")
            f.write(f"{len(results.reactions.RX)}\n")
            for i in range(len(results.reactions)):
                for j in range(len(results.reactions.RX)):
                    f.write(f"{results.reactions[i].RX[j]},{results.reactions[i].RY[j]},{results.reactions[i].RZ[j]},"
                            f"{results.reactions[i].MX[j]},{results.reactions[i].MY[j]},{results.reactions[i].MZ[j]}\n")

            f.write(f"{results.maxInternalForces[0][0][0]},{results.maxInternalForces[0][0][1]},\n")
            f.write(f"{results.maxInternalForces[0][1][0]},{results.maxInternalForces[0][1][1]},\n")
            f.write(f"{results.maxInternalForces[0][2][0]},{results.maxInternalForces[0][2][1]},\n")
            f.write(f"{results.maxInternalForces[1][0][0]},{results.maxInternalForces[1][0][1]},\n")
            f.write(f"{results.maxInternalForces[1][1][0]},{results.maxInternalForces[1][1][1]},\n")
            f.write(f"{results.maxInternalForces[1][2][0]},{results.maxInternalForces[1][2][1]},\n")

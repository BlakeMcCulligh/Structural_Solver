from tkinter import filedialog
import pandas as pd

def exportOptimizationResults(ressults, cost):
    """
    Exports optimization results to the specified Excel file.

    :param ressults: Member Cross-Section Results.
    :param cost: Cost of the best setup.
    """

    filePath = getExcelSavePath()

    if filePath is not None:

        df1 = pd.DataFrame(ressults, columns=["Member Index", "Cross-Section Type", "d", "b", "t", "A", "Iy", "Iz",
                                              "J"])
        df2 = pd.DataFrame({'Cost': cost})

        df1.to_excel(filePath, index=False, sheet_name="Results")
        df2.to_excel(filePath, index=False, sheet_name="Cost")

    else:
        print("Failed to export optimization results")

def exportResults(results):
    """
    Exports structural analysis results to the specified Excel file.

    :param results: Strucural Analysis Results Object to be exported.
    """

    filePath = getExcelSavePath()

    if filePath is not None:

        # Node Deflections
        nodeDeflectionArray = []
        ND = results.nodeDeflections
        for i in range(len(results.nodalDeflections)):
            for j in range(len(results.nodalDeflections.DX)):
                n = [j,i,ND[i].DX[j],ND[i].DY[j],ND[i].DZ[j],ND[i].RX[j],ND[i].RY[j],ND[i].RZ[j],]
                nodeDeflectionArray.append(n)
        nodeDeflection_df = pd.DataFrame(nodeDeflectionArray, columns=["MemberIndex","Load Case Index",
                                                                       "DX", "DY", "DZ", "RX", "RY", "RZ",])

        # Weight
        weight = {'Member Weight': results.weight}
        weight_df = pd.DataFrame(weight)
        overallWeight = {'Overall Weight': results.weight}
        overallWeight_df = pd.DataFrame(overallWeight)

        # Reactions
        reactionsArray = []
        R = results.reactions
        for i in range(len(results.reactions)):
            for j in range(len(results.reactions.RX)):
                n = [j,i,R[i].RX[j],R[i].RY[j],R[i].RZ[j],R[i].MX[j],R[i].MY[j],R[i].MZ[j]]
                reactionsArray.append(n)
        reactions_df = pd.DataFrame(reactionsArray, columns=["MemberIndex","Load Case Index",
                                                             "RX", "RY", "RZ", "MX", "MY", "MZ",])

        # Internal Maximum forces in each member
        MIF = results.maxInternalForces
        internalLoadsArray = [[MIF.FX,MIF.FX_case],[MIF.FY,MIF.FY_case],[MIF.FZ,MIF.FZ_case],[MIF.MX,MIF.MX_case],
                              [MIF.MY,MIF.MY_case],[MIF.MZ,MIF.MZ_case]]
        internalLoads_df = pd.DataFrame(internalLoadsArray, columns=["Load", "Governing Load Case"])

        # Writing to excel file
        nodeDeflection_df.to_excel(filePath, index=False, sheet_name="Node Deflections")
        weight_df.to_excel(filePath, index=True, sheet_name="Weight")
        overallWeight_df.to_excel(filePath, index=False, sheet_name="Overall Weight")
        reactions_df.to_excel(filePath, index=False, sheet_name="Reactions")
        internalLoads_df.to_excel(filePath, index=False, sheet_name="Internal Loads")

    else:
        print("Failed to export results")


def getExcelSavePath():
    """
    Returns specified path to Excel file.

    :return: Excel file path or None.
    """

    file_path = filedialog.asksaveasfilename(defaultextension=".xlsx",
                                             filetypes=[("Excel Files", "*.xlsx")])
    if file_path:return file_path
    return None
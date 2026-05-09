"""
Handels the exporting of 3D frame analysis results and optimization results to Excel files.
"""

from tkinter import filedialog
import pandas as pd

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

def export_optimization_results(results, cost):
    """
    Exports optimization results to the specified Excel file.

    :param results: Member cross-section results.
    :param cost: Cost of the best setup.
    """

    file_path = _get_excel_save_path()

    if file_path is not None:

        df1 = pd.DataFrame(results, columns=["Member Index", "Cross-Section _type", "d", "b", "t", "A", "Iy", "Iz",
                                              "J"])
        df2 = pd.DataFrame({'Cost': cost})

        with pd.ExcelWriter(file_path) as writer:
            df1.to_excel(writer, sheet_name="Results", index=False)
            df2.to_excel(writer, sheet_name="Cost", index=False)

    else:
        print("Failed to export optimization Results")

def export_results(results):
    """
    Exports structural analysis results to the specified Excel file.

    :param results: Strucural analysis results object to be exported.
    """

    file_path = _get_excel_save_path()

    if file_path is not None:

        # Node Deflections
        node_deflection_array = []
        ND = results.nodeDeflections
        for i in range(len(results.NodalDeflections)):
            for j in range(len(results.NodalDeflections.DX)):
                n = [j,i,ND[i].DX[j],ND[i].DY[j],ND[i].DZ[j],ND[i].RX[j],ND[i].RY[j],ND[i].RZ[j],]
                node_deflection_array.append(n)
        node_deflection_df = pd.DataFrame(node_deflection_array, columns=["MemberIndex","Load Case Index",
                                                                       "DX", "DY", "DZ", "RX", "RY", "RZ",])

        # Weight
        weight = {'Member Weight': results.Weight}
        weight_df = pd.DataFrame(weight)
        overall_weight = {'Overall Weight': results.Weight}
        overall_weight_df = pd.DataFrame(overall_weight)

        # Reactions
        reactions_array = []
        R = results.Reactions
        for i in range(len(results.Reactions)):
            for j in range(len(results.Reactions.RX)):
                n = [j,i,R[i].RX[j],R[i].RY[j],R[i].RZ[j],R[i].MX[j],R[i].MY[j],R[i].MZ[j]]
                reactions_array.append(n)
        reactions_df = pd.DataFrame(reactions_array, columns=["MemberIndex","Load Case Index",
                                                             "RX", "RY", "RZ", "MX", "MY", "MZ",])

        # Internal Maximum forces in each member
        MIF = results.MaxInternalForces
        internal_loads_array = [[MIF.FX,MIF.FX_case],[MIF.FY,MIF.FY_case],[MIF.FZ,MIF.FZ_case],[MIF.MX,MIF.MX_case],
                              [MIF.MY,MIF.MY_case],[MIF.MZ,MIF.MZ_case]]
        internal_loads_df = pd.DataFrame(internal_loads_array, columns=["Load", "Governing Load Case"])

        # Writing to Excel file
        with pd.ExcelWriter(file_path) as writer:
            node_deflection_df.to_excel(writer, index=False, sheet_name="Node Deflections")
            weight_df.to_excel(writer, index=True, sheet_name="Weight")
            overall_weight_df.to_excel(writer, index=False, sheet_name="Overall Weight")
            reactions_df.to_excel(writer, index=False, sheet_name="Reactions")
            internal_loads_df.to_excel(writer, index=False, sheet_name="Internal Loads")

    else:
        print("Failed to export Results")


def _get_excel_save_path():
    """
    Returns specified path to Excel file.

    :return: Excel file path or None.
    """

    file_path = filedialog.asksaveasfilename(defaultextension=".xlsx",
                                             filetypes=[("Excel Files", "*.xlsx")])
    if file_path:return file_path
    return None
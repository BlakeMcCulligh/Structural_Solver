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

        df1 = pd.DataFrame(ressults, columns=["Member Index", "Cross-Section Type", "d", "b", "t", "A", "Iy", "Iz", "J"])
        df2 = pd.DataFrame({'Cost': cost})

        df1.to_excel(filePath, index=False, sheet_name="Results")
        df2.to_excel(filePath, index=False, sheet_name="Cost")

    else:
        print("Failed to export optimization results")

def getExcelSavePath():
    """
    Returns specified path to Excel file.

    :return: Excel file path or None.
    """

    file_path = filedialog.asksaveasfilename(defaultextension=".xlsx",
                                             filetypes=[("Excel Files", "*.xlsx")])
    if file_path:return file_path
    return None
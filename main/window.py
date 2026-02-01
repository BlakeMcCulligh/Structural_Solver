
import tkinter as tk
from tkinter import filedialog

from CrossSectionAnalysis.main_CrossSectionAnalysis import getSectionProperties
from OpeningAndSaving.Opening import openTrussTopologyOptimizationExcel
from Sketch.main import startSketch


class MainWindow(tk.Frame):
    def __init__(self, root, openSketch):
        tk.Frame.__init__(self, root)

        self.root = root
        self.sketch = openSketch

        self.create_top_menu()

        exitButton = tk.Button(root, text="Complete Sketch", command=self.closeProgram)
        #self.canvas.create_window(100, 100, window=exitButton)

    def create_top_menu(self):
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)

        new_menu = tk.Menu(menubar, tearoff="off")
        new_menu.add_command(label='Truss Topology Optimization', command= self.newTrussTopologyOptimization)
        new_menu.add_command(label='Cross Section Analysis', command= self.startCrossSectionAnalysis)
        menubar.add_cascade(label="New", menu=new_menu)


        open_menu = tk.Menu(menubar, tearoff="off")
        open_menu.add_command(label='Truss Topology Optimization', command=self.openTrussTopologyOptimization)
        menubar.add_cascade(label="Open", menu=open_menu)

    def closeProgram(self):
        self.root.destroy()

    def newTrussTopologyOptimization(self):
        startSketch(self.root, self, "Closed Shape", self.setNodesTrussTopologyOptimization)

    def setNodesTrussTopologyOptimization(self, poly):
        print(poly)

        #TODO
        # lock polygon in grayed out background
        # draw Functions:
            # nodes
            # have it so lines can be drawn with nodes auto generated in a specified spasing along it
            # fill polgon with specified spacing
        # when sketch complete run truss iptimization
        # display results

    def startCrossSectionAnalysis(self):
        startSketch(self.root, self, "Closed Shape", self.returnCrossSectionAnalysis)

    def returnCrossSectionAnalysis(self, poly):
        section, properties = getSectionProperties(poly)
        print(properties)

    def openTrussTopologyOptimization(self):
        fileTypes = [("Excel Files", "*.xlsx")]
        selected_file = select_file_gui(fileTypes)
        openTrussTopologyOptimizationExcel(selected_file)


def select_file_gui(file_Types):
    """
    Opens a file explorer dialog using Tkinter and returns the selected file path.
    """
    # Create a Tk root window (but hide it)
    root = tk.Tk()
    root.withdraw()  # Hide the main window

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





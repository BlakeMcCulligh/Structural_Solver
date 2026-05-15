"""
Holds the object for the main window of the program.
"""

import tkinter as tk
from tkinter import filedialog
from pathlib import Path

from frame_3D_gui.opening import open_frame
from opening_saving.opening import open_truss_topology_optimization_excel, open_truss_cross_section_optimization_excel
from main.new_structure_pop_up import NewStructurePopUp
import frame_3D_gui.window

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

class MainWindow(tk.Frame):
    """
    Object containing the main window of the program.
    """

    def __init__(self, Root):
        """
        Constructor for the main window.

        :param Root: Root of the main window.
        """

        tk.Frame.__init__(self, Root)

        self.Root = Root

        self._create_top_menu()

        self._center_window()

    def _center_window(self):
        """
        Centers the window on the screen and sets its dimensions.
        """

        WIDTH = 1000
        HEIGHT = 800

        screen_width = self.Root.winfo_screenwidth()
        screen_height = self.Root.winfo_screenheight()

        if WIDTH > screen_width or HEIGHT > screen_height:
            WIDTH = screen_width * 0.8
            HEIGHT = screen_height * 0.8
        x = (screen_width // 2) - (WIDTH // 2)
        y = (screen_height // 2) - (HEIGHT // 2) - 50

        self.Root.geometry(f"{WIDTH}x{HEIGHT}+{x}+{y}")

    def _create_top_menu(self):
        """
        Creates the top menu with all the dropdowns for the window and handles all of their buttons
        """

        menubar = tk.Menu(self.Root)
        self.Root.config(menu=menubar)

        file_menu = tk.Menu(menubar, tearoff=False)
        file_menu.add_command(label='Open', command= self._open)
        menubar.add_cascade(label='File', menu=file_menu)

        new_menu = tk.Menu(menubar, tearoff=False)
        new_menu.add_command(label='Structure', command= self._new_structure)
        menubar.add_cascade(label="New", menu=new_menu)


        open_menu = tk.Menu(menubar, tearoff=False)
        open_menu.add_command(label='Truss Topology Optimization', command=self._open_truss_topology_optimization)
        open_menu.add_command(label='Truss Cross-Section Optimization', command=self.open_truss_cross_section_optimization)
        menubar.add_cascade(label="Analysis Excel File", menu=open_menu)

    def _open(self):
        """
        Lets the user select a save file and opens it.
        """

        file_types = [("Struct Frame files", "*.structframe")]
        file_path = select_file_gui(file_types)

        if Path(file_path).suffix == ".structframe":
            self.Root.destroy()

            root_widget = tk.Tk()
            root_widget.title('3D Frame')
            window = frame_3D_gui.window.MainWindow(root_widget)
            open_frame(window, file_path)
            window.Root.mainloop()

    @staticmethod
    def _open_truss_topology_optimization():
        """
        Opens an Excel file and runs a truss topology optimization from the excels Data.
        """

        fileTypes = [("Excel Files", "*.xlsx")]
        selected_file = select_file_gui(fileTypes)
        open_truss_topology_optimization_excel(selected_file)

    @staticmethod
    def open_truss_cross_section_optimization():
        """
        Opens an Excel file and runs a truss cross-section optimization from the excels Data.
        """

        file_types = [("Excel Files", "*.xlsx")]
        selected_file = select_file_gui(file_types)
        open_truss_cross_section_optimization_excel(selected_file)

    def _new_structure(self):
        """
        When a new structure is selected, create a pop-Up to get the info on what kind of pop-Up is wanted.
        """
        NewStructurePopUp(self.Root)

def select_file_gui(file_types):
    """
    Lets the user select a file of the type: file_types.

    :param file_types: List of file types that can be opened.
    """

    file_path = filedialog.askopenfilename(title="Select a file", filetypes= file_types)
    if file_path:return file_path
    else:return None
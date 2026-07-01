"""
Handels display frame for the 3D frame structure.
"""

from __future__ import annotations
import tkinter as tk
from tkinter import filedialog
from typing import TYPE_CHECKING

from rebuild_program_layout.frame_3D import Frame3D
from rebuild_program_layout.frame_3D_popups import AddingTablesWindow, OptimizationPopUp

if TYPE_CHECKING:
    from rebuild_program_layout import __main__

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

class Frame3DFrame(tk.Frame):
    """
    Object containing the frame for the 3D frame structure.
    """

    def __init__(self, parent: tk.Tk, controller: __main__.Structural_Solver):
        """
        Constructor: initializes the frame for the 3D frame structure

        :param parent: Object containing the Tkinter for the program.
        :type parent: tkinter.Tk
        :param controller: Object containing the main program.
        :type controller: __main__.Structural_Solver
        """

        super().__init__(parent)
        self._controller = controller
        self.root_window = self.winfo_toplevel()

        self.menubar = self._create_top_menu()

        self.file_path = None
        self.frame_3d = Frame3D(self._controller)

    def _create_top_menu(self) -> tk.Menu:
        """
        Creates the top menu of the window.

        return The top menu.
        :rtype: tkinter.Menu
        """

        menubar = tk.Menu(self.root_window)

        file_menu = tk.Menu(menubar, tearoff=False)
        file_menu.add_command(label='Save', command=self._save)
        file_menu.add_command(label='Save As', command=self._save_as)
        file_menu.add_command(label='Export Results', command=self._export_results)
        menubar.add_cascade(label="File", menu=file_menu)

        import_menu = tk.Menu(menubar, tearoff=False)
        import_menu.add_command(label='Nodes', command=lambda: self._handel_excel_imports("import_nodes"))
        import_menu.add_command(label='Members', command=lambda: self._handel_excel_imports("import_members"))
        import_menu.add_command(label='Materials', command=lambda: self._handel_excel_imports("import_materials"))
        import_menu.add_command(label='Supports', command=lambda: self._handel_excel_imports("import_supports"))
        import_menu.add_command(label='Releases', command=lambda: self._handel_excel_imports("import_releases"))
        import_menu.add_command(label='Node Loads', command=lambda: self._handel_excel_imports("import_node_loads"))
        import_menu.add_command(label='Members Point Loads',
                                command=lambda: self._handel_excel_imports("import_member_point_loads"))
        import_menu.add_command(label='Members Distributed Loads',
                                command=lambda: self._handel_excel_imports("import_member_dist_loads"))
        menubar.add_cascade(label="Import", menu=import_menu)

        add_menu = tk.Menu(menubar, tearoff=False)
        add_menu.add_command(label='Open Add Tabels', command=self._open_table_window)
        menubar.add_cascade(label="Add", menu=add_menu)

        analysis_menu = tk.Menu(menubar, tearoff=False)
        analysis_menu.add_command(label='Linear Analysis', command=self._linear_analysis)
        analysis_menu.add_command(label='Global Optimization', command=self._optimization_window)
        menubar.add_cascade(label="Analysis", menu=analysis_menu)

        return menubar

    def on_show(self) -> None:
        """
        Shows the top menu when this frame is displayed.
        """

        self.root_window.config(menu=self.menubar)

    def _exit(self) -> None:
        """
        Closes program.
        """

        self.root_window.quit()
        self.root_window.destroy()

    def _save_as(self) -> None:
        """
        Saves files under specified name and location.
        """

        self.file_path = _get_new_file_path(".structframe", [("Struct Frame files", "*.structframe")])

        if self.file_path is not None:
            self.file_path = self.frame_3d.save(self.file_path)

    def _save(self) -> None:
        """
        Saves files in current directory. If no current directory is specified, save under specified name and location.
        """

        if self.file_path is None:
            self._save_as()
        else:
            self.frame_3d.save(self.file_path)

    def _export_results(self) -> None:
        """
        Exports the analysis results to an Excel file.
        """

        self.file_path = _get_new_file_path(".xlsx", [("Excel Files", "*.xlsx")])

        if self.file_path is not None:
            self.frame_3d.export_results(self.file_path)

    def _handel_excel_imports(self, import_method_name: str) -> None:
        """
        Handels when a button to import info from an Excel file is chosen.

        :param import_method_name: What method is to be used. Must be a function within frame_3d.
        :type import_method_name: str
        """

        file_path = _get_existing_file_path([("Excel Files", "*.xlsx")])
        if file_path is not None:
            getattr(self.frame_3d, import_method_name)(file_path)

    def _open_table_window(self) -> None:
        """
        Opens the window of tables that show and alow the changing of the structure.
        """

        AddingTablesWindow(self.root_window, self)

    def _linear_analysis(self) -> None:
        """
        Runs the linear analysis tool.
        """

        self.frame_3d.linear_analysis()

    def _optimization_window(self) -> None:
        """
        Opens the window to input the data needed to run an optimization.
        """

        OptimizationPopUp(self.root_window, self)

def _get_new_file_path(file_type: str, file_type_name: list[tuple[str, str]] ) ->  str | None:
    """
    Gets the file path to save the frame to.

    :param file_type: String containing the file type aloued to be opend
    :type file_type: str
    :param file_type_name: list containing the names of the file types alowed and there file extension
    :type file_type_name: list[ tuple[ str ] ]
    :return: The file path to the save file or None if no file path is selected.
    :rtype: str or None
    """

    file_path = filedialog.asksaveasfilename(defaultextension=file_type,filetypes=file_type_name)
    if file_path: return file_path
    return None

def _get_existing_file_path(file_types: list[tuple[str, str]] ) ->  str | None:
    """
    Opens a file explorer dialogue using Tkinter and returns the selected file path.

    :param file_types:  list containing the names of the file types alowed and there file extension
    :type file_types: list[ tuple[ str ] ]
    :return: The selected file path or None if no file path is selected
    :rtype: str | None
    """

    file_path = filedialog.askopenfilename(title="Select a file", filetypes=file_types)
    if file_path:
        return file_path
    else:
        return None

"""
Handels the starting frame of the program when it is first lonched.
"""

from __future__ import annotations
import tkinter as tk
from tkinter import filedialog
from pathlib import Path
from typing import TYPE_CHECKING

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

class FrameMain(tk.Frame):
    """
    Object containing the starting frame of the program.
    """

    def __init__(self, parent: tk.Tk, controller: __main__.Structural_Solver):
        """
        Constructor: initializes the starting frame of the program.

        :param parent: Object containing the Tkinter for the program.
        :type parent: tkinter.Tk
        :param controller: Object containing the main program.
        :type controller: __main__.Structural_Solver
        """

        super().__init__(parent)
        self.controller = controller
        self.root_window = self.winfo_toplevel()

        self.menubar = self._create_top_menu()

    def _create_top_menu(self) -> tk.Menu:
        """
        Creates the top menu of the window.

        return The top menu.
        :rtype: tkinter.Menu
        """

        menubar = tk.Menu(self.root_window)

        file_menu = tk.Menu(menubar, tearoff=False)
        file_menu.add_command(label='Open', command= self._open)
        menubar.add_cascade(label='File', menu=file_menu)

        new_menu = tk.Menu(menubar, tearoff=False)
        new_menu.add_command(label='Structure', command= self._new_structure)
        menubar.add_cascade(label="New", menu=new_menu)

        return menubar

    def on_show(self) -> None:
        """
        Shows the top menu when this frame is displayed.
        """

        self.root_window.config(menu=self.menubar)

    def _open(self) -> None:
        """
        Lets the user select a save file and opens it.
        """

        file_types = [("Struct Frame files", "*.structframe")]
        file_path = select_file_gui(file_types)

        if file_path is not None:
            if Path(file_path).suffix == ".structframe":
                self.root_window.destroy()

                root_widget = tk.Tk()
                root_widget.title('3D Frame')
                self.controller.switch_to_frame_3d_frame(file_path)

    def _new_structure(self) -> None:
        """
        When a new structure is selected, create a pop-Up to get the info on what kind of structure is wanted.
        Then create that kind of structure and change the frame to it.
        """

        NewStructurePopUp(self.root_window, self)

def select_file_gui(file_types: list[tuple[str, str]] )  ->  str | None:
    """
    Lets the user select a file of the type: file_types.

    :param file_types: list containing the names of the file types alowed and there file extension.
    :type file_types: list[ tuple[ str ] ]
    :return: The selected file path or None if no file path is selected
    :rtype: str | None
    """

    file_path = filedialog.askopenfilename(title="Select a file", filetypes= file_types)
    if file_path: return file_path
    else: return None

class NewStructurePopUp(tk.Toplevel):
    """
    Object holding the custom popup window for the creation of a new structure.
    """

    def __init__(self, parent: tk.Tk, controller: FrameMain):
        """
        Constructor: initializes the custom popup window for the creation of a new structure.

        :param parent: Object containing the Tkinter for the program.
        :type parent: tkinter.Tk
        :param controller: Object containing the main program.
        :type controller: FrameMain
        """

        # Initialize the Toplevel window with its parent
        super().__init__(parent)
        self.controller = controller

        self.parent = parent
        WIDTH = 500
        HEIGHT = 170
        self._center_window(WIDTH, HEIGHT)
        self.title("Set Up New Structure")  # Set the title
        self.resizable(False, False)

        title = tk.Label(self, text="Set Up New Structure", font=('Helvetica', 20))
        w_title = title.winfo_reqwidth()
        h_title = title.winfo_reqheight()
        title.place(x=WIDTH / 2 - w_title / 2, y=0)

        dim_label = tk.Label(self, text="Dimentions: ", font=('Helvetica', 14))
        dim_label.place(x=0, y=h_title + 10)
        w_dim_label = dim_label.winfo_reqwidth()
        h_dim_label = dim_label.winfo_reqheight()

        DIM_OPTIONS = ["3D"]
        self._dim = tk.StringVar(self)
        self._dim.set(DIM_OPTIONS[0])
        dim_menu = tk.OptionMenu(self, self._dim, *DIM_OPTIONS)
        dim_menu.place(x=w_dim_label + 5, y=h_title + 10)

        type_label = tk.Label(self, text="Type: ", font=('Helvetica', 14))
        type_label.place(x=0, y=h_title + h_dim_label + 20)
        w_type_label = type_label.winfo_reqwidth()

        TYPE_OPTIONS = ["Frame"]
        self._type = tk.StringVar(self)
        self._type.set(TYPE_OPTIONS[0])
        type_menu = tk.OptionMenu(self, self._type, *TYPE_OPTIONS)
        type_menu.place(x=w_type_label + 5, y=h_title + h_dim_label + 20)

        tk.Button(self, text="Okay", command=self._end_pop_up).place(x=420, y=140)
        tk.Button(self, text="Cancel", command=self._cancel_pop_up).place(x=340, y=140)

        self.grab_set()
        self.parent.wait_window(self)

    def _center_window(self, width: int, height: int) -> None:
        """
        Centers the window on the screen and its dimensions.

        :param width: Width to make the window.
        :type width: int
        :param height: Height to make the window.
        :type height: int
        """

        main_window_x = self.parent.winfo_x()
        main_window_y = self.parent.winfo_y()
        main_window_width = self.parent.winfo_width()
        main_window_height = self.parent.winfo_height()

        x = (main_window_width // 2) + main_window_x - (width // 2)
        y = (main_window_height // 2) + main_window_y - (height // 2)

        self.geometry(f"{width}x{height}+{x}+{y}")

    def _end_pop_up(self) -> None:
        """
        Ends the pop-Up Window and creates the new structure and changes the main window to the structures window.
        """

        self.destroy()

        if self._dim.get() == "3D" and self._type.get() == "Frame":
            self.controller.controller.switch_to_frame_3d_frame()

    def _cancel_pop_up(self) -> None:
        """
        Cancels the creation of a new structure and closes the pop-Up Window.
        """

        self.destroy()

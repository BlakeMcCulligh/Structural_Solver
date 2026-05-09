"""
Holds everything to do with the pop-up to get info about the new structure being created.
"""

import tkinter as tk

import frame_3D_gui

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

class NewStructurePopUp:
    """
    Object for the new structure pop-up.
    """

    def __init__(self, _main_window_root):
        """
        Initializes the new structure pop-up.

        :param _main_window_root: The root of the main window creating the pop-up.
        """

        self._main_window_root = _main_window_root
        self._new_struct_popup_root = tk.Toplevel(_main_window_root)
        WIDTH = 500
        HEIGHT = 170
        self._center_window(WIDTH, HEIGHT)
        self._new_struct_popup_root.title("Set Up New Structure")  # Set the title
        self._new_struct_popup_root.resizable(False, False)

        title = tk.Label(self._new_struct_popup_root, text="Set Up New Structure", font=('Helvetica', 20))
        w_title = title.winfo_reqwidth()
        h_title = title.winfo_reqheight()
        title.place(x=WIDTH/2-w_title/2, y=0)

        dim_label = tk.Label(self._new_struct_popup_root, text="_dim: ", font=('Helvetica', 14))
        dim_label.place(x=0, y=h_title+10)
        w_dim_label = dim_label.winfo_reqwidth()
        h_dim_label = dim_label.winfo_reqheight()

        DIM_OPTIONS = ["2D", "3D"]
        self._dim = tk.StringVar(self._new_struct_popup_root)
        self._dim.set(DIM_OPTIONS[0])
        dim_menu = tk.OptionMenu(self._new_struct_popup_root, self._dim, *DIM_OPTIONS)
        dim_menu.place(x=w_dim_label+5, y=h_title+10)

        type_label = tk.Label(self._new_struct_popup_root, text="_type: ", font=('Helvetica', 14))
        type_label.place(x=0, y=h_title + h_dim_label + 20)
        w_type_label = type_label.winfo_reqwidth()

        TYPE_OPTIONS = ["Truss", "Frame"]
        self._type = tk.StringVar(self._new_struct_popup_root)
        self._type.set(TYPE_OPTIONS[0])
        type_menu = tk.OptionMenu(self._new_struct_popup_root, self._type, *TYPE_OPTIONS)
        type_menu.place(x=w_type_label+5, y=h_title + h_dim_label + 20)

        tk.Button(self._new_struct_popup_root, text="Okay", command=self._end_pop_up).place(x=420, y=140)
        tk.Button(self._new_struct_popup_root, text="Cancel", command=self._cancel_pop_up).place(x=340, y=140)

        self._new_struct_popup_root.grab_set()
        _main_window_root.wait_window(self._new_struct_popup_root)

    def _center_window(self, width, height):
        """
        Centers the window on the screen and its dimentions.

        :param width: Width to make the window.
        :param height: Height to make the window.
        """

        main_window_x = self._main_window_root.winfo_x()
        main_window_y = self._main_window_root.winfo_y()
        main_window_width = self._main_window_root.winfo_width()
        main_window_height = self._main_window_root.winfo_height()

        x = (main_window_width // 2) + main_window_x - (width // 2)
        y = (main_window_height // 2) + main_window_y - (height // 2)

        self._new_struct_popup_root.geometry(f"{width}x{height}+{x}+{y}")

    def _end_pop_up(self):
        """
        Ends the pop-Up Window and creates the new structure and changes the main window to the structures window.
        """

        self._new_struct_popup_root.destroy()

        if self._dim.get() == "3D" and self._type.get() == "Frame":
            self._new_struct_popup_root.destroy()
            self._main_window_root.destroy()

            root_widget = tk.Tk()
            root_widget.title('3D Frame')
            frame_3D_gui.window.MainWindow(root_widget)
            root_widget.mainloop()

    def _cancel_pop_up(self):
        """
        Cancels the creation of a new structure and closes the pop-Up Window.
        """

        self._new_struct_popup_root.destroy()

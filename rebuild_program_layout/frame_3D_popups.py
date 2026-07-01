"""
Handels top windows for the 3D frame.
"""

from __future__ import annotations
import tkinter as tk
from tkinter import ttk
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from rebuild_program_layout.frame_3d_frame import Frame3DFrame

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

from rebuild_program_layout.text_validate import validate_float, validate_index, validate_bool

class AddingTablesWindow(tk.Toplevel):
    """
    Object containing the top window that alows for the editing and changing of the 3D frame structure.
    """

    def __init__(self, parent: tk.Tk, controller: Frame3DFrame):
        """
        Constructor for the adding tables window.

        :param parent: Object containing the Tkinter for the program.
        :type parent: tkinter.Tk
        :param controller: Object containing the display frame for the 3D frame structure.
        :type controller: Frame3DFrame
        """

        super().__init__(parent)
        self.title("Sheets")
        self.geometry("800x700")
        self.resizable(False, False)

        self._controller = controller
        self._root_window = self.winfo_toplevel()

        self._val_float = (self._root_window.register(validate_float), '%P')
        self._val_index = (self._root_window.register(validate_index), '%P')
        self._val_bool = (self._root_window.register(validate_bool), '#P')

        self._tabs_notebook = None
        self._frames = []
        self._tables = []
        self._boxes = []
        self._buttons = []

        self._setup_window()

    def _setup_window(self) -> None:
        """
        Sets up all the tabs and tables of the window.
        """

        self._tabs_notebook = ttk.Notebook(self)
        self._tabs_notebook.pack(expand=True, fill="both")

        for i in range(8): self._frames.append(ttk.Frame(self._tabs_notebook))

        # Nodes
        self._tabs_notebook.add(self._frames[0], text="Nodes")
        self._tables.append(self._create_table(self._frames[0], ("Index", "X", "Y", "Z")))
        list_validate = [self._val_float,self._val_float,self._val_float]
        list_label = ["X:", "Y:", "Z:"]
        list_label_location = [[325,297],[325,347],[325,397]]
        list_box_location = [[350,300],[350,350],[350,400]]
        self._create_boxes(0, list_validate, list_label, list_label_location, list_box_location)

        # Materials
        self._tabs_notebook.add(self._frames[1], text="Materials")
        self._tables.append(self._create_table(self._frames[1], ("Index", "E", "G", "nu", "rho", "fy")))
        list_validate = [self._val_float, self._val_float, self._val_float, self._val_float, self._val_float]
        list_label = ["E:", "G:", "Poisons Ratio:", "Density:", "Yield Strength:"]
        list_label_location = [[125, 297], [325, 297], [40, 347], [70,397], [30,447]]
        list_box_location = [[150,300], [350,300], [150,350], [150,400], [150,450]]
        self._create_boxes(1, list_validate, list_label, list_label_location, list_box_location)

        # Members
        self._tabs_notebook.add(self._frames[2], text="Members")
        self._tables.append(self._create_table(self._frames[2], ("Index", "i Node", "j Node", "Material Id",
                                                                 "Set C.S.", "A", "Iy", "Iz", "J")))
        list_validate = [self._val_index, self._val_index, self._val_index, None, self._val_float, self._val_float,
                         self._val_float, self._val_float]
        list_label = ["i Node:", "j Node:", "Material Id:","Set C.S.:", "A:", "Iy:", "Iz:", "J:"]
        list_label_location = [[90,297], [290,297], [60,347], [275,347], [125,397], [325,397], [125,447], [325,447]]
        list_box_location = [[150,300], [350,300], [150,350], [350,350], [150,400], [350,400], [150,450], [350,450]]
        self._create_boxes(2, list_validate, list_label, list_label_location, list_box_location)

        # Supports
        self._tabs_notebook.add(self._frames[3], text="Supports")
        self._tables.append(self._create_table(self._frames[3], ("Index", "i Node", "D X", "D Y", "D Z.", "R X",
                                                                  "R Y", "R Z")))
        list_validate = [self._val_index, None, None, None, None, None, None]
        list_label = ["i Node:", "D X:", "D Y:", "D Z.:", "R X:", "R Y:", "R Z:"]
        list_label_location = [[90, 297], [312, 297], [115, 347], [315, 347], [112, 397], [315, 397], [115, 447]]
        list_box_location = [[150, 300], [350, 300], [150, 350], [350, 350], [150, 400], [350, 400], [150, 450]]
        self._create_boxes(3, list_validate, list_label, list_label_location, list_box_location)

        # Releases
        self._tabs_notebook.add(self._frames[4], text="Releases")
        self._tables.append(self._create_table(self._frames[4], ("Index", "i Member", "i DX", "i DY", "i DZ",
                                                                  "i RX", "i RY", "i RZ", "j DX", "j DY", "j DZ",
                                                                  "j RX",
                                                                  "j RY", "j RZ")))
        list_validate = [self._val_index, None, None, None, None, None, None, None, None, None, None, None, None]
        list_label = ["i Member:", "i DX:", "i DY:", "i DZ:", "i RX:", "i RY:", "i RZ:", "j DX:", "j DY:", "j DZ:",
                      "j RX:", "j RY:", "j RZ:"]
        list_label_location = [[70, 297], [310, 297], [110, 347], [310, 347], [110, 397], [310, 397], [110, 447],
                               [310, 447], [110, 497], [310, 497], [110, 547], [310, 547], [110, 597]]
        list_box_location = [[150, 300], [350, 300], [150, 350], [350, 350], [150, 400], [350, 400], [150, 450],
                             [350, 450], [150, 500], [350, 500], [150, 550], [350, 550], [150, 600]]
        self._create_boxes(4, list_validate, list_label, list_label_location, list_box_location)

        # Nodal Loads
        self._tabs_notebook.add(self._frames[5], text="Node Loads")
        self._tables.append(self._create_table(self._frames[5], ("Index", "i Node", "P X", "P Y", "P Z.",
                                                                    "M X", "M Y", "M Z", "Cases")))
        list_validate = [self._val_index, self._val_float, self._val_float, self._val_float, self._val_float,
                         self._val_float, self._val_float, self._val_index]
        list_label = ["i Node:", "P X:", "P Y:", "P Z.:", "M X:", "M Y:", "M Z:", "Cases:"]
        list_label_location = [[90, 297], [310, 297], [110, 347], [310, 347], [110, 397], [310, 397], [110, 447],
                               [290, 447]]
        list_box_location = [[150, 300], [350, 300], [150, 350], [350, 350], [150, 400], [350, 400], [150, 450],
                             [350, 450]]
        self._create_boxes(5, list_validate, list_label, list_label_location, list_box_location)

        # Member Point Loads
        self._tabs_notebook.add(self._frames[6], text="Member Point Loads")
        self._tables.append(self._create_table(self._frames[6], ("Index", "i Member", "x", "P X", "P Y",
                                                                            "P Z.", "M X", "M Y", "M Z", "Cases")))
        list_validate = [self._val_index, self._val_float, self._val_float, self._val_float, self._val_float,
                         self._val_float, self._val_float, self._val_float, self._val_index]
        list_label = ["i Member:", "x:", "P X:", "P Y:", "P Z:", "M X:", "M Y:", "M Z:", "Cases:"]
        list_label_location = [[70, 297], [325, 297], [110, 347], [310, 347], [110, 397], [310, 397], [110, 447],
                               [310, 447], [70, 497]]
        list_box_location = [[150, 300], [350, 300], [150, 350], [350, 350], [150, 400], [350, 400], [150, 450],
                             [350, 450], [150, 500]]
        self._create_boxes(6, list_validate, list_label, list_label_location, list_box_location)


        # Member Distributed Loads
        self._tabs_notebook.add(self._frames[7], text="Member Distributed Loads")
        self._tables.append(self._create_table(self._frames[7], ("Index", "i Member", "x 1", "x 2", "wx 1",
                                                                           "wx 2", "wy 1", "wy 2", "wz 1", "wz 2",
                                                                           "Case")))
        list_validate = [self._val_index, self._val_float, self._val_float, self._val_float, self._val_float,
                         self._val_float, self._val_float, self._val_float, self._val_float, self._val_index]
        list_label = ["i Member:", "x 1:", "x 2:", "wx 1:", "wx 2:", "wy 1:", "wy 2:", "wz 1:", "wz 2:", "Case:"]
        list_label_location = [[70, 297], [310, 297], [110, 347], [300, 347], [100, 397], [300, 397], [100, 447],
                               [300, 447], [100, 497], [290, 497]]
        list_box_location = [[150, 300], [350, 300], [150, 350], [350, 350], [150, 400], [350, 400], [150, 450],
                             [350, 450], [150, 500], [350, 500]]
        self._create_boxes(7, list_validate, list_label, list_label_location, list_box_location)

    def _create_table(self, tab: ttk.Frame, headings: tuple) -> ttk.Treeview:
        """
        Creates an input table.

        :param tab: Root of tab to put the table on.
        :type tab: tkinter.tkk.Frame
        :param headings: Headings to put on the table.
        :type headings: tuple
        :return: table
        :rtype: tkinter.ttk.Treeview
        """

        table_frame = ttk.Frame(tab)
        table_frame.pack(pady=10)

        table = ttk.Treeview(table_frame, columns=headings, show='headings', height=12)

        # adding scroll bar to table
        scrollbar = ttk.Scrollbar(table_frame, orient="vertical", command=table.yview)
        scrollbar.pack(side="right", fill="y")
        table.configure(yscrollcommand=scrollbar.set)

        # adding headings to table
        for i in range(len(headings)):
            table.heading(headings[i], text=headings[i])
            table.column(headings[i], width=50)

        # Updates the input boxes to be the current values in the selected table row
        table.bind("<<TreeviewSelect>>", self._get_selected_row)

        table.pack()

        return table

    def _create_boxes(self, tab_index: int, list_validate: list[tuple[str, str] | None], list_label: list[str],
                      list_label_location: list[list[int]], list_box_location: list[list[int]]) -> None:
        """
        Creates all the input boxes, lables, and buttons for the tab.

        :param tab_index: Index of the tab to add to.
        :type tab_index: int
        :param list_validate: list of what each value is to varidated as. If None, create a True or False drop down.
        :type list_validate: list[tuple[str, str] | None]
        :param list_label: List of lables for each input box.
        :type list_label: list[str]
        :param list_label_location: list of locations for each label.
        :type list_label_location: list[list[int]]
        :param list_box_location: list of locations for each input box.
        :type list_box_location: list[list[int]]
        """

        self._boxes.append([])
        for i in range(len(list_validate)):
            if list_validate[i] is None:
                self._boxes[tab_index].append(ttk.Combobox(self._frames[tab_index], values=("True", "False"),
                                                           state="readonly"))
            else:
                self._boxes[tab_index].append(tk.Entry(self._frames[tab_index], validate='key',
                                                       validatecommand=list_validate[i]))
            label = tk.Label(self._frames[tab_index], text=list_label[i], font=('Helvetica', 12))
            label.place(x=list_label_location[i][0], y=list_label_location[i][1])
            self._boxes[tab_index][i].place(x=list_box_location[i][0], y=list_box_location[i][1])
        self._buttons.append([])
        self._buttons[tab_index].append(tk.Button(self._frames[tab_index], text="Edit", command=self._edit_values))
        self._buttons[tab_index].append(tk.Button(self._frames[tab_index], text="Add", command=self._add_values))
        self._buttons[tab_index][0].place(x=550, y=300)
        self._buttons[tab_index][1].place(x=550, y=350)

    def _get_selected_row(self, event) -> None:
        """
        Updates the input boxes to be the current values in the selected table row.
        """

        num_col = [3,5,8,7,13,8,9,10]
        t_id = self._tabs_notebook.index("current")

        row_id = self._tables[t_id].focus()
        row_info = self._tables[t_id].item(row_id).get('values')
        for i in range(num_col[t_id]):
            self._boxes[t_id][i].delete(0, 'end')
            self._boxes[t_id][i].insert(0, row_info[i + 1])

    def _edit_values(self) -> None:
        """
        Changes the values in the selected row of the table to be the values in the input boxes.
        """

        num_col = [3, 5, 8, 7, 13, 8, 9, 10]
        t_id = self._tabs_notebook.index("current")

        row_id = self._tables[t_id].focus()
        row_info = self._tables[t_id].item(row_id).get('values')
        new_val = [row_info[0]]
        for i in range(num_col[t_id]): new_val.append(self._boxes[t_id][i].get()) # getting values from text boxes

        data_new = convert_str_to_float_or_bool(new_val, t_id, num_col)

        self._update_structure(t_id, data_new)

    def _add_values(self) -> None:
        """
        Adds a new row to the table with the values in the input boxes.
        """

        num_col = [3, 5, 8, 7, 13, 8, 9, 10]
        t_id = self._tabs_notebook.index("current")
        num_row = len(self._tables[t_id].get_children())
        new_val = [num_row]
        for i in range(num_col[t_id]): new_val.append(self._boxes[t_id][i].get()) # getting values from text boxes

        data_new = convert_str_to_float_or_bool(new_val, t_id, num_col)

        self._update_structure(t_id, data_new)

    def _update_structure(self, table, data_new):
        print("update_structure", data_new)
        pass
        #TODO update structue and table

def convert_str_to_float_or_bool(new_val: list[str | int], t_id: int, num_col: list[int]) -> list[float | bool]:
    """
    Coverts a list of strings and intagers into a list of a mixture of floats and booleans.

    :param new_val: list of strings to be converted.
    :type new_val: list[str or int]
    :param t_id: Index of the current tab.
    :type t_id: int
    :param num_col: Number of columns in each tab.
    :type num_col: list[int]
    :return: List of converted values
    :rtype: list[float | bool]
    """

    data_new = []
    for i in range(num_col[t_id]):
        if new_val[i + 1] == "":
            data_new.append(0.0)
        elif new_val[i + 1] == "True":
            data_new.append(True)
        elif new_val[i + 1] == "False":
            data_new.append(False)
        else:
            data_new.append(float(new_val[i + 1]))

    return data_new

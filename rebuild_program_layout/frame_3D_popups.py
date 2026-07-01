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

class OptimizationPopUp(tk.Toplevel):
    """
    Object containing everything to do with the pop-up when a cross-section optimization is selected.
    """
    def __init__(self, parent: tk.Tk, controller: Frame3DFrame):
        """
        Initializes the optimization pop-up object.

        :param parent: Object containing the Tkinter for the program.
        :type parent: tkinter.Tk
        :param controller: Object containing the display frame for the 3D frame structure.
        :type controller: Frame3DFrame
        """

        super().__init__(parent)

        self._controller = controller
        self._root_window = parent.winfo_toplevel()

        self._center_window(680, 500)
        self.title("Set Up New Structure")  # Set the title
        self.resizable(False, False)

        title = tk.Label(self, text="Set Up Optimization", font=('Helvetica', 20))
        w_title = title.winfo_reqwidth()
        h_title = title.winfo_reqheight()
        title.place(x=680 / 2 - w_title / 2, y=0)

        self._val_float = (self._root_window.register(validate_float), '%P')
        self._val_index = (self._root_window.register(validate_index), '%P')
        self._val_bool = (self._root_window.register(validate_bool), '#P')

        self._member_group_input_frame = None
        self._member_group_input_table = None
        self._set_member_group_input(h_title)

        self._group_settings_frame = None
        self._group_settings_table = None
        self._group_settings_input(h_title)

        self._cost_function_string_box = None
        self._dropdown_weight = None
        self._dropdown_reactions = None
        self._dropdown_internal_forces = None
        self._cost_function_input(h_title)

        ttk.Button(self, text="Okay", command=self._end_pop_up).place(x=580, y=470)
        ttk.Button(self, text="Cancel", command=self._cancel_pop_up).place(x=500, y=470)


    def _center_window(self, width: int, height: int) -> None:
        """
        Centers the pop-up window on the screen and sets its size.

        :param width: Width to set the pop-Up window to.
        :type width: int
        :param height: Height to set the pop-Up window to.
        :type height: int
        """

        main_window_x = self._root_window.winfo_x()
        main_window_y = self._root_window.winfo_y()
        main_window_width = self._root_window.winfo_width()
        main_window_height = self._root_window.winfo_height()

        x = (main_window_width // 2) + main_window_x - (width // 2)
        y = (main_window_height // 2) + main_window_y - (height // 2)

        self.geometry(f"{width}x{height}+{x}+{y}")

    def _set_member_group_input(self, h_title: int) -> None:
        """
        Initializes the tables that allows the input of values used to define what member group is applied to each
        non-cross-section set member.

        :param h_title: Highest of the title of the window.
        :type h_title: int
        """

        num_not_set_members, indices_not_set_members = self.get_non_set_members_data()

        d_label = tk.Label(self, text="Member Groups: ", font=('Helvetica', 14))
        d_label.place(x=10, y=h_title + 10)
        h_d_label = d_label.winfo_reqheight()

        self._member_group_input_frame = ttk.Frame(self, width=50, height=300)
        self._member_group_input_frame.place(x=10, y =h_title + h_d_label + 10)

        self._member_group_input_table = ttk.Treeview(self._member_group_input_frame,
                                                      columns=("Member Index", "Group Index"),
                                                      show='headings', height=10)

        scrollbar = ttk.Scrollbar(self._member_group_input_frame, orient="vertical",
                                  command=self._member_group_input_table.yview)
        scrollbar.pack(side="right", fill="y")

        self._member_group_input_table.configure(yscrollcommand=scrollbar.set)

        self._member_group_input_table.heading("Member Index", text="Member Index")
        self._member_group_input_table.column("Member Index", width=90)
        self._member_group_input_table.heading("Group Index", text="Group Index")
        self._member_group_input_table.column("Group Index", width=75)

        for i in range(len(indices_not_set_members)):
            self._member_group_input_table.insert("", "end", values=(indices_not_set_members[i], i%4))
            #TODO ^ used for testing perpuses ^

        self._member_group_input_table.bind('<Double-1>', self._edit_cell_member_group_input_table)

        self._member_group_input_table.pack()

    def _edit_cell_member_group_input_table(self, event) -> None:
        """
        Handles the cell editing of the member assignment to groups table.

        :param event: Event object from the double click of a cell.
        """

        # Identify which row and column were clicked
        region = self._member_group_input_table.identify_region(event.x, event.y)
        if region != "cell": return

        column = self._member_group_input_table.identify_column(event.x)
        column_index = int(column[1:]) - 1  # Convert '#1' to 0
        if column_index == 0: return
        iid = self._member_group_input_table.identify_row(event.y)

        # Get cell coordinates and dimensions
        x, y, width, height = self._member_group_input_table.bbox(iid, column)

        # Create the Entry widget overlay
        entry = ttk.Entry(self._member_group_input_frame, validate='key', validatecommand=self._val_index)
        entry.insert(0, self._member_group_input_table.item(iid)['values'][column_index])
        entry.place(x=x, y=y, width=width, height=height)
        entry.focus_set()

        # noinspection PyShadowingNames
        def save_edit(event) -> None:
            """
            Updates Table value, group setting table, and remove Entry

            :param event: Event object from the double click of a cell.
            """

            groups = list(sorted(set(self._get_data_from_member_group_input())))

            # Update Treeview value and remove Entry
            new_values = list(self._member_group_input_table.item(iid)['values'])
            newVal = entry.get()
            new_values[column_index] = newVal
            self._member_group_input_table.item(iid, values=new_values)
            entry.destroy()

            # updating group settings
            if int(newVal) not in groups:
                added = False
                for i in range(len(groups)-1):
                    if groups[i + 1] > int(newVal) > groups[i]:
                        self._group_settings_table.insert("", i + 1,
                                                          values=(newVal, "TubeHSS", 1, 5, "", "", 0.05, 0.5))
                        added = True
                if not added:
                    self._group_settings_table.insert("", "end", values=(newVal, "TubeHSS", 1, 5, "", "", 0.05, 0.5))

        entry.bind('<Return>', save_edit)
        entry.bind('<FocusOut>', lambda e: entry.destroy())

    def _group_settings_input(self, h_title: int) -> None:
        """
        Initializes the tables that allows the input of values used to define the cross-section groups.

        :param h_title: Highest of the title of the window.
        :type h_title: int
        """

        groups = sorted(set(self._get_data_from_member_group_input()))

        d_label = tk.Label(self, text="Groups: ", font=('Helvetica', 14))
        d_label.place(x=200, y=h_title + 10)
        h_d_label = d_label.winfo_reqheight()

        self._group_settings_frame = ttk.Frame(self, width=50, height=300)
        self._group_settings_frame.place(x=200, y=h_title + h_d_label + 10)

        self._group_settings_table = ttk.Treeview(self._group_settings_frame, columns=("Group Index", "Group _type",
                                                                                       "Min d","Max d", "Min b",
                                                                                       "Max b", "Min t", "Max t"),
                                                  show='headings', height=10)

        scrollbar = ttk.Scrollbar(self._group_settings_frame, orient="vertical",
                                  command=self._group_settings_table.yview)
        scrollbar.pack(side="right", fill="y")

        self._group_settings_table.configure(yscrollcommand=scrollbar.set)

        self._group_settings_table.heading("Group Index", text="Group Index")
        self._group_settings_table.column("Group Index", width=75)
        self._group_settings_table.heading("Group _type", text="Group _type")
        self._group_settings_table.column("Group _type", width=75)
        self._group_settings_table.heading("Min d", text="Min d")
        self._group_settings_table.column("Min d", width=50)
        self._group_settings_table.heading("Max d", text="Max d")
        self._group_settings_table.column("Max d", width=50)
        self._group_settings_table.heading("Min b", text="Min b")
        self._group_settings_table.column("Min b", width=50)
        self._group_settings_table.heading("Max b", text="Max b")
        self._group_settings_table.column("Max b", width=50)
        self._group_settings_table.heading("Min t", text="Min t")
        self._group_settings_table.column("Min t", width=50)
        self._group_settings_table.heading("Max t", text="Max t")
        self._group_settings_table.column("Max t", width=50)

        for i in groups:
            self._group_settings_table.insert("", "end", values=(i, "TubeHSS", 1, 5, "", "", 0.05, 0.5))

        self._group_settings_table.bind('<Double-1>', self._edit_cell_group_settings_table)

        self._group_settings_table.pack()

    def _edit_cell_group_settings_table(self, event) -> None:
        """
        Handles the cell editing of the member group settings table.

        :param event: Event object from the duble click of a cell.
        """

        # Identify which row and column were clicked
        region = self._group_settings_table.identify_region(event.x, event.y)
        if region != "cell": return

        column = self._group_settings_table.identify_column(event.x)
        column_index = int(column[1:]) - 1  # Convert '#1' to 0

        if column_index == 0: return

        iid = self._group_settings_table.identify_row(event.y)

        # Get cell coordinates and dimensions
        x, y, width, height = self._group_settings_table.bbox(iid, column)

        if column_index == 1:
            # Create the dropdown widget overlay
            dropdown = ttk.Combobox(self._group_settings_frame, values=("Angle", "RectHSS", "SquareHSS", "TubeHSS"),
                                    state="readonly")
            dropdown.set(self._group_settings_table.item(iid)['values'][column_index])
            dropdown.place(x=x, y=y, width=width, height=height)

            # noinspection PyShadowingNames
            def save_edit(event) -> None:
                """
                Updates Table value and remove Entry

                :param event: Event object from the double click of a cell.
                """

                new_values = list(self._group_settings_table.item(iid)['values'])
                new_values[column_index] = dropdown.get()
                self._group_settings_table.item(iid, values=new_values)
                dropdown.destroy()

            dropdown.bind('<<ComboboxSelected>>', save_edit)
            dropdown.bind('<FocusOut>', lambda e: dropdown.destroy())

        else:
            # Create the Entry widget overlay
            entry = ttk.Entry(self._group_settings_frame, validate='key', validatecommand=self._val_float)
            entry.insert(0, self._group_settings_table.item(iid)['values'][column_index])
            entry.place(x=x, y=y, width=width, height=height)
            entry.focus_set()

            # noinspection PyShadowingNames
            def save_edit(event):
                """
                Updates Table value and remove Entry

                :param event: Event object from the double click of a cell.
                """

                new_values = list(self._group_settings_table.item(iid)['values'])
                new_values[column_index] = entry.get()
                self._group_settings_table.item(iid, values=new_values)
                entry.destroy()

            entry.bind('<Return>', save_edit)
            entry.bind('<FocusOut>', lambda e: entry.destroy())

    def _get_data_from_member_group_input(self) -> list:
        """
        Gets data from the member group definition table.

        :return: Data from the member group definition table.
        :rtype: list
        """

        return [self._member_group_input_table.item(item)['values'][1]
                for item in self._member_group_input_table.get_children()]

    def _cost_function_input(self, hTitle: int) -> None:
        """
        Initializes and runs everything related to the input of the cost function.

        :param hTitle: Highest of the title of the window.
        :type hTitle: int
        """

        d_label = tk.Label(self, text="Cost Function: ", font=('Helvetica', 11))
        d_label.place(x=10, y=hTitle + 300)

        def validate_input_box(p: str) -> bool:
            """
            Restricts what can be input into the cost function text box. Ran every time a character is changed.

            :param p: Updated string of what is typed in the cross-section text box.
            :type p: str
            :return: Whether the new string is allowed or not (True or False).
            :rtype: bool
            """

            if p == "":
                return True

            keyword = None
            for i in range(len(p)):
                if (p[i].isdigit() or p[i] == "[" or p[i] == "]" or p[i] == "(" or p[i] == ")" or p[i] == " " or
                        p[i] == "+" or p[i] == "-" or p[i] == "*" or p[i] == "/" or p[i] == "."):
                    keyword = None
                elif (p[i] == "D" or p[i] == "R" or p[i] == "W" or p[i] == "I" or p[i] == "m" or p[i] == "a"  or
                      p[i] == "s") and keyword is None:
                    keyword = p[i]
                elif (p[i] == "X" or p[i] == "Y" or p[i] == "Z") and (keyword == "D" or keyword == "R"):
                    keyword = keyword + p[i]
                elif ((p[i] == "e" and keyword == "W") or (p[i] == "i" and keyword == "We") or
                      (p[i] == "g" and keyword == "Wei") or (p[i] == "h" and keyword == "Weig") or
                      (p[i] == "t" and keyword == "Weigh")):
                    keyword = keyword + p[i]
                elif ((p[i] == "e" and keyword == "R") or (p[i] == "a" and keyword == "Re") or
                      (p[i] == "c" and keyword == "Rea") or (p[i] == "t" and keyword == "Reac") or
                      (p[i] == "i" and keyword == "React") or (p[i] == "o" and keyword == "Reacti") or
                      (p[i] == "n" and keyword == "Reactio") or (p[i] == "s" and keyword == "Reaction")):
                    keyword = keyword + p[i]
                elif ((p[i] == "n" and keyword == "I") or (p[i] == "t" and keyword == "In") or
                      (p[i] == "e" and keyword == "Int") or (p[i] == "r" and keyword == "Inte") or
                      (p[i] == "n" and keyword == "Inter") or (p[i] == "a" and keyword == "Intern") or
                      (p[i] == "a" and keyword == "Intern") or (p[i] == "l" and keyword == "Interna") or
                      (p[i] == "F" and keyword == "Internal") or (p[i] == "o" and keyword == "InternalF") or
                      (p[i] == "r" and keyword == "InternalFo") or (p[i] == "c" and keyword == "InternalFor") or
                      (p[i] == "e" and keyword == "InternalForc") or (p[i] == "s" and keyword == "InternalForce")):
                    keyword = keyword + p[i]
                elif ((p[i] == "a" and keyword == "m") or (p[i] == "x" and keyword == "ma")  or
                      (p[i] == "i" and keyword == "m") or (p[i] == "n" and keyword == "mi")  or
                      (p[i] == "b" and keyword == "a")  or (p[i] == "s" and keyword == "ab")):
                    keyword = keyword + p[i]
                elif (p[i] == "u" and keyword == "s") or (p[i] == "m" and keyword == "su"):
                    keyword = keyword + p[i]
                else:
                    return False
            return True

        vcmd = self.register(validate_input_box)

        self._cost_function_string_box = ttk.Entry(self, validate="key", validatecommand=(vcmd, '%P'))
        self._cost_function_string_box.place(x=10, y=hTitle + 325, width=350, height=20)


        ttk.Button(self, text="Information", command=self._info_pop_up).place(x=10, y=hTitle + 360)

        d_label = tk.Label(self, text="Weight Needed: ", font=('Helvetica', 11))
        d_label.place(x=400, y=hTitle + 300)
        h_d_label = d_label.winfo_reqheight()

        self._dropdown_weight = ttk.Combobox(self, values=("False", "True"), state="readonly")
        self._dropdown_weight.set("False")
        self._dropdown_weight.place(x=580, y=hTitle + 300, width=75, height=h_d_label)

        d_label = tk.Label(self, text="Reactions Needed: ", font=('Helvetica', 11))
        d_label.place(x=400, y=hTitle + 325)
        h_d_label = d_label.winfo_reqheight()

        self._dropdown_reactions = ttk.Combobox(self, values=("False", "True"), state="readonly")
        self._dropdown_reactions.set("False")
        self._dropdown_reactions.place(x=580, y=hTitle + 325, width=75, height=h_d_label)

        d_label = tk.Label(self, text="Internal Forces Needed: ", font=('Helvetica', 11))
        d_label.place(x=400, y=hTitle + 350)
        h_d_label = d_label.winfo_reqheight()

        self._dropdown_internal_forces = ttk.Combobox(self, values=("False", "True"), state="readonly")
        self._dropdown_internal_forces.set("False")
        self._dropdown_internal_forces.place(x=580, y=hTitle + 350, width=75, height=h_d_label)

    def _info_pop_up(self) -> None:
        """
        Creates a pop-up that displays info about the cost function and how to enter use it.
        """

        info_pop_up = tk.Toplevel(self)
        info_pop_up.title("Cost Function Info")
        WIDTH = 500
        HEIGHT = 500
        main_window_x = self._root_window.winfo_x()
        main_window_y = self._root_window.winfo_y()
        main_window_width = self._root_window.winfo_width()
        main_window_height = self._root_window.winfo_height()
        x = (main_window_width // 2) + main_window_x - (WIDTH // 2)
        y = (main_window_height // 2) + main_window_y - (HEIGHT // 2)
        info_pop_up.geometry(f"{WIDTH}x{HEIGHT}+{x}+{y}")
        info_pop_up.resizable(False, False)

        text_area = tk.Text(info_pop_up, height=HEIGHT, width=WIDTH, font=('Helvetica', 11))
        text_area.place(x=0,y=0 )
        text_area.pack()

        string = (" Cost Function Info\n ----- Variables -----\n 'DX': Node Deflection in the x direction.\n "
                  "'DY': Node Deflection in the y direction.\n 'DZ': Node Deflection in the z direction.\n "
                  "'RX': Node Rotation Deflection in the x direction.\n "
                  "'RY': Node Rotation Deflection in the y direction.\n "
                  "'RZ': Node Rotation Deflection in the z direction.\n  "
                  "All Deflections: Index 1: Load Cases, Index 2: Node.\n \n 'Weight': Weight of the structure.\n "
                  "'Reactions': Index 1: Load Cases,\n"
                  "                     Index 2: Node,\n"
                  "                     Index 3: Direction (X, Y, Z, RX, RY, RZ).\n"
                  " 'InternalForces': Index 1: Direction Type (Forces, Moments),\n"
                  "                            Index 2: Direction (X, Y, Z),\n"
                  "                            Index 3: Info (Magnitude, Governing Case).\n"
                  " To be able to use 'Weight','Reactions','InternalForces' variables they must\n"
                  " be selected to run.\n \n \n  ----- Functions -----\n ' + ': addition\n ' - ': subtraction\n"
                  " ' * ': multiplication\n ' / ': division\n ' ** ': exponent\n"
                  " ' max(list) ': finds the maximum value of list.\n ' min(list) ': finds the minimum value of list.\n"
                  " ' abs(float) ': finds the absolute value of float.")

        text_area.insert(tk.END, string)
        text_area.config(state="disabled")

    def _end_pop_up(self) -> None:
        """
        Gets input table, closes the pop-up, and starts the cross-section optimization.
        """

        group_assignments = [self._member_group_input_table.item(item)['values'][1]
                            for item in self._member_group_input_table.get_children()]
        group_types = [self._group_settings_table.item(item)['values'][1]
                      for item in self._group_settings_table.get_children()]
        d_low_bounds = [self._group_settings_table.item(item)['values'][2]
                        for item in self._group_settings_table.get_children()]
        d_high_bounds = [self._group_settings_table.item(item)['values'][3]
                       for item in self._group_settings_table.get_children()]
        b_low_bounds = [self._group_settings_table.item(item)['values'][4]
                        for item in self._group_settings_table.get_children()]
        b_high_bounds = [self._group_settings_table.item(item)['values'][5]
                       for item in self._group_settings_table.get_children()]
        t_low_bounds = [self._group_settings_table.item(item)['values'][6]
                        for item in self._group_settings_table.get_children()]
        t_high_bounds = [self._group_settings_table.item(item)['values'][7]
                       for item in self._group_settings_table.get_children()]

        cost_function = self._cost_function_string_box.get()
        weight_run = self._dropdown_weight.get()
        reaction_run = self._dropdown_reactions.get()
        internal_forces_run = self._dropdown_internal_forces.get()

        lower_bounds = []
        upper_bounds = []
        for i in range(len(group_types)):
            lower_bounds.append(float(d_low_bounds[i]))
            upper_bounds.append(float(d_high_bounds[i]))
            if group_types[i] == "Angle" or group_types[i] == "RectHSS":
                lower_bounds.append(float(b_low_bounds[i]))
                upper_bounds.append(float(b_high_bounds[i]))
            lower_bounds.append(float(t_low_bounds[i]))
            upper_bounds.append(float(t_high_bounds[i]))

        self.destroy()

        (self._controller.frame_3d.GlobalOptimization(group_assignments, group_types, lower_bounds, upper_bounds,
                                                      cost_function, weight_run, reaction_run, internal_forces_run))

    def _cancel_pop_up(self) -> None:
        """
        Removes the pop-up window.
        """

        self.destroy()

    def get_non_set_members_data(self) -> list[int | list[int]]:
        """
        Gets the number and indexes of what member cross-sections are not set and need to be optimized.

        :return: Number of non-set members. Indices of what member are not set.
        :rtype: list[int | list[int]]
        """

        #TODO

        num = 0
        non_indices = []

        #self._controller.frame_3d

        return [num, non_indices]
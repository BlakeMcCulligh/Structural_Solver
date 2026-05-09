"""
Holds the object that handels everything to do with the optimization pop-up.
"""

import tkinter as tk
from tkinter import ttk

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

class OptimizationPopUp:
    """
    Object containing everything to do with the pop-up when a cross-section optimization is selected.
    """

    def __init__(self, MainWindow, MainWindowRoot, Data):
        """
        Initalizes the optimization pop-up object.

        :param MainWindow: Object storing the main window.
        :param MainWindowRoot: Root of the main window.
        :param Data: Data needed to define the frame.
        """

        self._main_window = MainWindow
        self._main_window_root = MainWindowRoot
        self._data = Data
        self._root_pop_up = tk.Toplevel(MainWindowRoot)
        WIDTH = 680
        HEIGHT = 500
        self._center_window(WIDTH, HEIGHT)
        self._root_pop_up.title("Set Up New Structure")  # Set the title
        self._root_pop_up.resizable(False, False)

        title = tk.Label(self._root_pop_up, text="Set Up Optimization", font=('Helvetica', 20))
        w_title = title.winfo_reqwidth()
        h_title = title.winfo_reqheight()
        title.place(x=WIDTH/2-w_title/2, y=0)

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


        ttk.Button(self._root_pop_up, text="Okay", command=self._end_pop_up).place(x=580, y=470)
        ttk.Button(self._root_pop_up, text="Cancel", command=self._cancel_pop_up).place(x=500, y=470)

    def _center_window(self, width, height):
        """
        Centers the pop-up window on the screen and sets its size.

        :param width: Width to set the pop-Up window to.
        :param height: Height to set the pop-Up window to.
        """

        main_window_x = self._main_window_root.winfo_x()
        main_window_y = self._main_window_root.winfo_y()
        main_window_width = self._main_window_root.winfo_width()
        main_window_height = self._main_window_root.winfo_height()

        x = (main_window_width // 2) + main_window_x - (width // 2)
        y = (main_window_height // 2) + main_window_y - (height // 2)

        self._root_pop_up.geometry(f"{width}x{height}+{x}+{y}")

    def _set_member_group_input(self, h_title):
        """
        Initalizes the tables that alows the input of values used to define what member group is applied to each
        non-cross-section set member.

        :param h_title:  Highet of the title of the window.
        """

        num_not_set_members, indices_not_set_members = get_non_set_members_data(self._data)

        d_label = tk.Label(self._root_pop_up, text="Member Groups: ", font=('Helvetica', 14))
        d_label.place(x=10, y=h_title + 10)
        h_d_label = d_label.winfo_reqheight()

        self._member_group_input_frame = ttk.Frame(self._root_pop_up, width=50, height=300)
        self._member_group_input_frame.place(x=10, y =h_title + h_d_label + 10)

        self._member_group_input_table = ttk.Treeview(self._member_group_input_frame, columns=("Member Index", "Group Index"),
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
            self._member_group_input_table.insert("", "end", values=(indices_not_set_members[i], i))

        self._member_group_input_table.bind('<Double-1>', self._edit_cell_member_group_input_table)

        self._member_group_input_table.pack()

    def _edit_cell_member_group_input_table(self, event):
        """
        Handels the cell editing of the member assignment to groups table.

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
        entry = ttk.Entry(self._member_group_input_frame)
        entry.insert(0, self._member_group_input_table.item(iid)['values'][column_index])
        entry.place(x=x, y=y, width=width, height=height)
        entry.focus_set()

        # noinspection PyShadowingNames
        def save_edit(event):
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
                        self._group_settings_table.insert("", i + 1, values=(newVal, "TubeHSS", 1, 5, "", "", 0.05, 0.5))
                        added = True
                if not added:
                    self._group_settings_table.insert("", "end", values=(newVal, "TubeHSS", 1, 5, "", "", 0.05, 0.5))

        entry.bind('<Return>', save_edit)
        entry.bind('<FocusOut>', lambda e: entry.destroy())

    def _group_settings_input(self, h_title):
        """
        Initalizes the tables that alows the input of values used to define the cross-section groups.

        :param h_title: Highet of the title of the window.
        """

        groups = sorted(set(self._get_data_from_member_group_input()))

        d_label = tk.Label(self._root_pop_up, text="Groups: ", font=('Helvetica', 14))
        d_label.place(x=200, y=h_title + 10)
        h_d_label = d_label.winfo_reqheight()

        self._group_settings_frame = ttk.Frame(self._root_pop_up, width=50, height=300)
        self._group_settings_frame.place(x=200, y=h_title + h_d_label + 10)

        self._group_settings_table = ttk.Treeview(self._group_settings_frame, columns=("Group Index", "Group _type", "Min d",
                                                                                 "Max d", "Min b", "Max b", "Min t",
                                                                                 "Max t"),
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

        # for i in range(20):
        #     self._group_settings_table.insert("", "end", values=(i, 1))

        self._group_settings_table.bind('<Double-1>', self._edit_cell_group_settings_table)

        self._group_settings_table.pack()

    def _edit_cell_group_settings_table(self, event):
        """
        Handels the cell editing of the member group settings table.

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
            def save_edit(event):
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
            entry = ttk.Entry(self._group_settings_frame)
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

    def _get_data_from_member_group_input(self):
        """
        Gets data from the member group definition table.

        :return: Data from the member group definition table.
        """

        return [self._member_group_input_table.item(item)['values'][1]
                for item in self._member_group_input_table.get_children()]

    def _cost_function_input(self, hTitle):
        """
        Initalizes and runs everything related to the input of the cost function.

        :param hTitle: Highet of the title of the window.
        """

        d_label = tk.Label(self._root_pop_up, text="Cost Function: ", font=('Helvetica', 11))
        d_label.place(x=10, y=hTitle + 300)

        def validate_input_box(p):
            """
            Restricts what can be input into the cost function text box. Ran every time a character is changed.

            :param p: Updated string of what is typed in the cross-section text box.
            :return: Whether the new string is alowed or not (True or False).
            """

            if p == "":
                return True

            keyword = None
            for i in range(len(p)):
                if p[i].isdigit() or p[i] == "[" or p[i] == "]" or p[i] == "(" or p[i] == ")" or p[i] == " " or p[i] == "+" or p[i] == "-" or p[i] == "*" or p[i] == "/" or p[i] == ".":
                    keyword = None
                elif (p[i] == "D" or p[i] == "R" or p[i] == "W" or p[i] == "I" or p[i] == "m" or p[i] == "a") and keyword is None:
                    keyword = p[i]
                elif (p[i] == "X" or p[i] == "Y" or p[i] == "Z") and (keyword == "D" or keyword == "R"):
                    keyword = keyword + p[i]
                elif (p[i] == "e" and keyword == "W") or (p[i] == "i" and keyword == "We") or (p[i] == "g" and keyword == "Wei") or (p[i] == "h" and keyword == "Weig") or (p[i] == "t" and keyword == "Weigh"):
                    keyword = keyword + p[i]
                elif (p[i] == "e" and keyword == "R") or (p[i] == "a" and keyword == "Re") or (p[i] == "c" and keyword == "Rea") or (p[i] == "t" and keyword == "Reac") or (p[i] == "i" and keyword == "React") or (p[i] == "o" and keyword == "Reacti") or (p[i] == "n" and keyword == "Reactio") or (p[i] == "s" and keyword == "Reaction"):
                    keyword = keyword + p[i]
                elif (p[i] == "n" and keyword == "I") or (p[i] == "t" and keyword == "In") or (p[i] == "e" and keyword == "Int") or (p[i] == "r" and keyword == "Inte") or (p[i] == "n" and keyword == "Inter") or (p[i] == "a" and keyword == "Intern") or (p[i] == "a" and keyword == "Intern") or (p[i] == "l" and keyword == "Interna") or (p[i] == "F" and keyword == "Internal") or (p[i] == "o" and keyword == "InternalF") or (p[i] == "r" and keyword == "InternalFo") or (p[i] == "c" and keyword == "InternalFor") or (p[i] == "e" and keyword == "InternalForc") or (p[i] == "s" and keyword == "InternalForce"):
                    keyword = keyword + p[i]
                elif (p[i] == "a" and keyword == "m") or (p[i] == "x" and keyword == "ma")  or (p[i] == "i" and keyword == "m") or (p[i] == "n" and keyword == "mi")  or (p[i] == "b" and keyword == "a")  or (p[i] == "s" and keyword == "ab"):
                    keyword = keyword + p[i]
                else:
                    return False
            return True

        vcmd = self._root_pop_up.register(validate_input_box)

        self._cost_function_string_box = ttk.Entry(self._root_pop_up, validate="key", validatecommand=(vcmd, '%P'))
        self._cost_function_string_box.place(x=10, y=hTitle + 325, width=350, height=20)


        ttk.Button(self._root_pop_up, text="Information", command=self._info_pop_up).place(x=10, y=hTitle + 360)

        d_label = tk.Label(self._root_pop_up, text="Weight Needed: ", font=('Helvetica', 11))
        d_label.place(x=400, y=hTitle + 300)
        h_d_label = d_label.winfo_reqheight()

        self._dropdown_weight = ttk.Combobox(self._root_pop_up, values=("False", "True"), state="readonly")
        self._dropdown_weight.set("False")
        self._dropdown_weight.place(x=580, y=hTitle + 300, width=75, height=h_d_label)

        d_label = tk.Label(self._root_pop_up, text="Reactions Needed: ", font=('Helvetica', 11))
        d_label.place(x=400, y=hTitle + 325)
        h_d_label = d_label.winfo_reqheight()

        self._dropdown_reactions = ttk.Combobox(self._root_pop_up, values=("False", "True"), state="readonly")
        self._dropdown_reactions.set("False")
        self._dropdown_reactions.place(x=580, y=hTitle + 325, width=75, height=h_d_label)

        d_label = tk.Label(self._root_pop_up, text="Internal Forces Needed: ", font=('Helvetica', 11))
        d_label.place(x=400, y=hTitle + 350)
        h_d_label = d_label.winfo_reqheight()

        self._dropdown_internal_forces = ttk.Combobox(self._root_pop_up, values=("False", "True"), state="readonly")
        self._dropdown_internal_forces.set("False")
        self._dropdown_internal_forces.place(x=580, y=hTitle + 350, width=75, height=h_d_label)

    def _info_pop_up(self):
        """
        Creates a pop-up that displays info about the cost function and how to enter use it.
        """

        info_pop_up = tk.Toplevel(self._root_pop_up)
        info_pop_up.title("Cost Function Info")
        WIDTH = 500
        HEIGHT = 500
        main_window_x = self._main_window_root.winfo_x()
        main_window_y = self._main_window_root.winfo_y()
        main_window_width = self._main_window_root.winfo_width()
        main_window_height = self._main_window_root.winfo_height()
        x = (main_window_width // 2) + main_window_x - (WIDTH // 2)
        y = (main_window_height // 2) + main_window_y - (HEIGHT // 2)
        info_pop_up.geometry(f"{WIDTH}x{HEIGHT}+{x}+{y}")
        info_pop_up.resizable(False, False)

        text_area = tk.Text(info_pop_up, height=HEIGHT, width=WIDTH, font=('Helvetica', 11))
        text_area.place(x=0,y=0 )
        text_area.pack()

        string = (" Cost Function Info\n ----- Vaiables -----\n 'DX': Node Deflection in the x direction.\n "
                  "'DY': Node Deflection in the y direction.\n 'DZ': Node Deflection in the z direction.\n "
                  "'RX': Node Rotation Deflection in the x direction.\n "
                  "'RY': Node Rotation Deflection in the y direction.\n "
                  "'RZ': Node Rotation Deflection in the z direction.\n  "
                  "All Deflections: Index 1: Load Casses, Index 2: Node.\n \n 'Weight': Weight of the structure.\n "
                  "'Reactions': Index 1: Load Casses,\n"
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

    def _end_pop_up(self):
        """
        Gets input table, closes the pop-up, and starts the ccross-section optimization.
        """

        self._root_pop_up.destroy()

        group_assignments = [self._member_group_input_table.item(item)['values'][1]
                            for item in self._member_group_input_table.get_children()]
        group_types = [self._group_settings_table.item(item)['values'][1]
                      for item in self._group_settings_table.get_children()]
        d_high_bounds = [self._group_settings_table.item(item)['values'][2]
                        for item in self._group_settings_table.get_children()]
        d_low_bounds = [self._group_settings_table.item(item)['values'][3]
                       for item in self._group_settings_table.get_children()]
        b_high_bounds = [self._group_settings_table.item(item)['values'][4]
                        for item in self._group_settings_table.get_children()]
        b_low_bounds = [self._group_settings_table.item(item)['values'][5]
                       for item in self._group_settings_table.get_children()]
        t_high_bounds = [self._group_settings_table.item(item)['values'][6]
                        for item in self._group_settings_table.get_children()]
        t_low_bounds = [self._group_settings_table.item(item)['values'][7]
                       for item in self._group_settings_table.get_children()]

        cost_function = self._cost_function_string_box.get()
        weight_run = self._dropdown_weight.get()
        reaction_run = self._dropdown_reactions.get()
        internal_forces_run = self._dropdown_internal_forces.get()

        lower_bounds = []
        upper_bounds = []
        for i in range(len(group_types)):
            lower_bounds.append(d_low_bounds)
            upper_bounds.append(d_high_bounds)
            if group_types[i] == "Angle" or group_types[i] == "RectHSS":
                lower_bounds.append(b_low_bounds)
                upper_bounds.append(b_high_bounds)
            lower_bounds.append(t_low_bounds)
            upper_bounds.append(t_high_bounds)

        self._main_window.GlobalOptimization(self, group_assignments, group_types, lower_bounds, upper_bounds, cost_function,
                                             weight_run, reaction_run, internal_forces_run)

    def _cancel_pop_up(self):
        """
        Removes the pop-up window.
        """

        self._root_pop_up.destroy()

def get_non_set_members_data(d):
    """
    Gets the number and indeces of what member cross-sections are not set and need to be optimized.

    :param d: Info that defines the frame.
    :return: Number of non-set members. Indeces of what member are not set.
    """

    num = 0
    non_indices = []
    for i in range(len(d.Members[0])):
        if not d.Members[3]:
            num += 1
            non_indices.append(i)

    return num, non_indices
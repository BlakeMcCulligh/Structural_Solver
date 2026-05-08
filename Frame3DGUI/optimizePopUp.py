import tkinter as tk
from tkinter import ttk

class OptimizationPopUp:
    """
    Object containing everything to do with the pop-up when a cross-section optimization is selected.
    """

    def __init__(self, mainwindow, root, data):
        """
        Initalizes the optimization pop-up object.

        :param mainwindow: Object storing the main window.
        :param root: Root of the main window.
        :param data: Data needed to define the frame.
        """

        self.mainwindow = mainwindow
        self.root = root
        self.data = data
        self.top = tk.Toplevel(root)
        width = 680
        height = 500
        self.centerWindow(width, height)
        self.top.title("Set Up New Structure")  # Set the title
        self.top.resizable(False, False)

        title = tk.Label(self.top, text="Set Up Optimization", font=('Helvetica', 20))
        wTitle = title.winfo_reqwidth()
        hTitle = title.winfo_reqheight()
        title.place(x=width/2-wTitle/2, y=0)

        self.MemberGroupInputFrame = None
        self.MemberGroupInputTable = None
        self.setMemberGroupInput(hTitle)

        self.GroupSettingsFrame = None
        self.GroupSettingsTable = None
        self.GroupSettingsInput(hTitle)

        self.costFunctionStringBox = None
        self.dropdownWeight = None
        self.dropdownReactions = None
        self.dropdownInternalForces = None
        self.costFunctionInput(hTitle)


        ttk.Button(self.top, text="Okay", command=self.endPopUp).place(x=580, y=470)
        ttk.Button(self.top, text="Cancel", command=self.cancelPopUp).place(x=500, y=470)

    def centerWindow(self, width, height):
        """
        Centers the pop-up window on the screen and sets its size.

        :param width: Width to set the pop-up window to.
        :param height: Height to set the pop-up window to.
        """

        mainwindow_x = self.root.winfo_x()
        mainwindow_y = self.root.winfo_y()
        mainwindow_width = self.root.winfo_width()
        mainwindow_height = self.root.winfo_height()

        x = (mainwindow_width // 2) + mainwindow_x - (width // 2)
        y = (mainwindow_height // 2) + mainwindow_y - (height // 2)

        self.top.geometry(f"{width}x{height}+{x}+{y}")

    def setMemberGroupInput(self, hTitle):
        """
        Initalizes the tables that alows the input of values used to define what member group is applied to each
        non-cross-section set member.

        :param hTitle:  Highet of the title of the window.
        """

        numNotSetMembers, indicesNotSetMembers = getNonSetMembersData(self.data)

        DLabel = tk.Label(self.top, text="Member Groups: ", font=('Helvetica', 14))
        DLabel.place(x=10, y=hTitle + 10)
        hDLabel = DLabel.winfo_reqheight()

        self.MemberGroupInputFrame = ttk.Frame(self.top, width=50, height=300)
        self.MemberGroupInputFrame.place(x=10, y = hTitle + hDLabel + 10)

        self.MemberGroupInputTable = ttk.Treeview(self.MemberGroupInputFrame, columns=("Member Index","Group Index"),
                                                  show='headings', height=10)

        scrollbar = ttk.Scrollbar(self.MemberGroupInputFrame, orient="vertical",
                                  command=self.MemberGroupInputTable.yview)
        scrollbar.pack(side="right", fill="y")

        self.MemberGroupInputTable.configure(yscrollcommand=scrollbar.set)

        self.MemberGroupInputTable.heading("Member Index", text="Member Index")
        self.MemberGroupInputTable.column("Member Index", width=90)
        self.MemberGroupInputTable.heading("Group Index", text="Group Index")
        self.MemberGroupInputTable.column("Group Index", width=75)

        for i in range(len(indicesNotSetMembers)):
            self.MemberGroupInputTable.insert("", "end", values=(indicesNotSetMembers[i], i))

        self.MemberGroupInputTable.bind('<Double-1>', self.edit_cell_MemberGroupInputTable)

        self.MemberGroupInputTable.pack()

    def edit_cell_MemberGroupInputTable(self, event):
        """
        Handels the cell editing of the member assignment to groups table.

        :param event: Event object from the double click of a cell.
        """

        # Identify which row and column were clicked
        region = self.MemberGroupInputTable.identify_region(event.x, event.y)
        if region != "cell": return

        column = self.MemberGroupInputTable.identify_column(event.x)
        column_index = int(column[1:]) - 1  # Convert '#1' to 0
        if column_index == 0: return
        iid = self.MemberGroupInputTable.identify_row(event.y)

        # Get cell coordinates and dimensions
        x, y, width, height = self.MemberGroupInputTable.bbox(iid, column)

        # Create the Entry widget overlay
        entry = ttk.Entry(self.MemberGroupInputFrame)
        entry.insert(0, self.MemberGroupInputTable.item(iid)['values'][column_index])
        entry.place(x=x, y=y, width=width, height=height)
        entry.focus_set()

        # noinspection PyShadowingNames
        def save_edit(event):
            """
            Updates Table value, group setting table, and remove Entry

            :param event: Event object from the double click of a cell.
            """

            groups = list(sorted(set(self.getDataFromMemberGroupInput())))

            # Update Treeview value and remove Entry
            new_values = list(self.MemberGroupInputTable.item(iid)['values'])
            newVal = entry.get()
            new_values[column_index] = newVal
            self.MemberGroupInputTable.item(iid, values=new_values)
            entry.destroy()

            # updating group settings
            if int(newVal) not in groups:
                added = False
                for i in range(len(groups)-1):
                    if groups[i + 1] > int(newVal) > groups[i]:
                        self.GroupSettingsTable.insert("", i+1, values=(newVal, "TubeHSS", 1, 5, "", "", 0.05, 0.5))
                        added = True
                if not added:
                    self.GroupSettingsTable.insert("", "end", values=(newVal, "TubeHSS", 1, 5, "", "", 0.05, 0.5))

        entry.bind('<Return>', save_edit)
        entry.bind('<FocusOut>', lambda e: entry.destroy())

    def GroupSettingsInput(self, hTitle):
        """
        Initalizes the tables that alows the input of values used to define the cross-section groups.

        :param hTitle: Highet of the title of the window.
        """

        groups = sorted(set(self.getDataFromMemberGroupInput()))

        DLabel = tk.Label(self.top, text="Groups: ", font=('Helvetica', 14))
        DLabel.place(x=200, y=hTitle + 10)
        hDLabel = DLabel.winfo_reqheight()

        self.GroupSettingsFrame = ttk.Frame(self.top, width=50, height=300)
        self.GroupSettingsFrame.place(x=200, y=hTitle + hDLabel + 10)

        self.GroupSettingsTable = ttk.Treeview(self.GroupSettingsFrame, columns=("Group Index", "Group Type", "Min d",
                                                                                 "Max d", "Min b", "Max b", "Min t",
                                                                                 "Max t"),
                                                  show='headings', height=10)

        scrollbar = ttk.Scrollbar(self.GroupSettingsFrame, orient="vertical",
                                  command=self.GroupSettingsTable.yview)
        scrollbar.pack(side="right", fill="y")

        self.GroupSettingsTable.configure(yscrollcommand=scrollbar.set)

        self.GroupSettingsTable.heading("Group Index", text="Group Index")
        self.GroupSettingsTable.column("Group Index", width=75)
        self.GroupSettingsTable.heading("Group Type", text="Group Type")
        self.GroupSettingsTable.column("Group Type", width=75)
        self.GroupSettingsTable.heading("Min d", text="Min d")
        self.GroupSettingsTable.column("Min d", width=50)
        self.GroupSettingsTable.heading("Max d", text="Max d")
        self.GroupSettingsTable.column("Max d", width=50)
        self.GroupSettingsTable.heading("Min b", text="Min b")
        self.GroupSettingsTable.column("Min b", width=50)
        self.GroupSettingsTable.heading("Max b", text="Max b")
        self.GroupSettingsTable.column("Max b", width=50)
        self.GroupSettingsTable.heading("Min t", text="Min t")
        self.GroupSettingsTable.column("Min t", width=50)
        self.GroupSettingsTable.heading("Max t", text="Max t")
        self.GroupSettingsTable.column("Max t", width=50)


        for i in groups:
            self.GroupSettingsTable.insert("", "end", values=(i, "TubeHSS", 1, 5, "", "", 0.05, 0.5))

        # for i in range(20):
        #     self.GroupSettingsTable.insert("", "end", values=(i, 1))

        self.GroupSettingsTable.bind('<Double-1>', self.edit_cell_GroupSettingsTable)

        self.GroupSettingsTable.pack()

    def edit_cell_GroupSettingsTable(self, event):
        """
        Handels the cell editing of the member group settings table.

        :param event: Event object from the duble click of a cell.
        """

        # Identify which row and column were clicked
        region = self.GroupSettingsTable.identify_region(event.x, event.y)
        if region != "cell": return

        column = self.GroupSettingsTable.identify_column(event.x)
        column_index = int(column[1:]) - 1  # Convert '#1' to 0

        if column_index == 0: return

        iid = self.GroupSettingsTable.identify_row(event.y)

        # Get cell coordinates and dimensions
        x, y, width, height = self.GroupSettingsTable.bbox(iid, column)

        if column_index == 1:
            # Create the dropdown widget overlay
            dropdown = ttk.Combobox(self.GroupSettingsFrame, values=("Angle", "RectHSS", "SquareHSS", "TubeHSS"),
                                    state="readonly")
            dropdown.set(self.GroupSettingsTable.item(iid)['values'][column_index])
            dropdown.place(x=x, y=y, width=width, height=height)

            # noinspection PyShadowingNames
            def save_edit(event):
                """
                Updates Table value and remove Entry

                :param event: Event object from the double click of a cell.
                """

                new_values = list(self.GroupSettingsTable.item(iid)['values'])
                new_values[column_index] = dropdown.get()
                self.GroupSettingsTable.item(iid, values=new_values)
                dropdown.destroy()

            dropdown.bind('<<ComboboxSelected>>', save_edit)
            dropdown.bind('<FocusOut>', lambda e: dropdown.destroy())

        else:
            # Create the Entry widget overlay
            entry = ttk.Entry(self.GroupSettingsFrame)
            entry.insert(0, self.GroupSettingsTable.item(iid)['values'][column_index])
            entry.place(x=x, y=y, width=width, height=height)
            entry.focus_set()

            # noinspection PyShadowingNames
            def save_edit(event):
                """
                Updates Table value and remove Entry

                :param event: Event object from the double click of a cell.
                """

                new_values = list(self.GroupSettingsTable.item(iid)['values'])
                new_values[column_index] = entry.get()
                self.GroupSettingsTable.item(iid, values=new_values)
                entry.destroy()

            entry.bind('<Return>', save_edit)
            entry.bind('<FocusOut>', lambda e: entry.destroy())

    def getDataFromMemberGroupInput(self):
        """
        Gets data from the member group definition table.

        :return: Data from the member group definition table.
        """

        return [self.MemberGroupInputTable.item(item)['values'][1]
                for item in self.MemberGroupInputTable.get_children()]

    def costFunctionInput(self, hTitle):
        """
        Initalizes and runs everything related to the input of the cost function.

        :param hTitle: Highet of the title of the window.
        """

        DLabel = tk.Label(self.top, text="Cost Function: ", font=('Helvetica', 11))
        DLabel.place(x=10, y=hTitle + 300)

        def validate_inputBox(p):
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

        vcmd = self.top.register(validate_inputBox)

        self.costFunctionStringBox = ttk.Entry(self.top, validate="key", validatecommand=(vcmd, '%P'))
        self.costFunctionStringBox.place(x=10, y=hTitle + 325, width=350, height=20)


        ttk.Button(self.top, text="Information", command=self.infoPopUp).place(x=10, y=hTitle + 360)

        DLabel = tk.Label(self.top, text="Weight Needed: ", font=('Helvetica', 11))
        DLabel.place(x=400, y=hTitle + 300)
        hDLabel = DLabel.winfo_reqheight()

        self.dropdownWeight = ttk.Combobox(self.top, values=("False", "True"), state="readonly")
        self.dropdownWeight.set("False")
        self.dropdownWeight.place(x=580, y=hTitle + 300, width=75, height=hDLabel)

        DLabel = tk.Label(self.top, text="Reactions Needed: ", font=('Helvetica', 11))
        DLabel.place(x=400, y=hTitle + 325)
        hDLabel = DLabel.winfo_reqheight()

        self.dropdownReactions = ttk.Combobox(self.top, values=("False", "True"), state="readonly")
        self.dropdownReactions.set("False")
        self.dropdownReactions.place(x=580, y=hTitle + 325, width=75, height=hDLabel)

        DLabel = tk.Label(self.top, text="Internal Forces Needed: ", font=('Helvetica', 11))
        DLabel.place(x=400, y=hTitle + 350)
        hDLabel = DLabel.winfo_reqheight()

        self.dropdownInternalForces = ttk.Combobox(self.top, values=("False", "True"), state="readonly")
        self.dropdownInternalForces.set("False")
        self.dropdownInternalForces.place(x=580, y=hTitle + 350, width=75, height=hDLabel)

    def infoPopUp(self):
        """
        Creates a pop-up that displays info about the cost function and how to enter use it.
        """

        infoPopUp = tk.Toplevel(self.top)
        infoPopUp.title("Cost Function Info")
        width = 500
        height = 500
        mainwindow_x = self.root.winfo_x()
        mainwindow_y = self.root.winfo_y()
        mainwindow_width = self.root.winfo_width()
        mainwindow_height = self.root.winfo_height()
        x = (mainwindow_width // 2) + mainwindow_x - (width // 2)
        y = (mainwindow_height // 2) + mainwindow_y - (height // 2)
        infoPopUp.geometry(f"{width}x{height}+{x}+{y}")
        infoPopUp.resizable(False, False)

        text_area = tk.Text(infoPopUp, height=height, width=width, font=('Helvetica', 11))
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

    def endPopUp(self):
        """
        Gets input table, closes the pop-up, and starts the ccross-section optimization.
        """

        self.top.destroy()

        GroupAssignments = [self.MemberGroupInputTable.item(item)['values'][1]
                            for item in self.MemberGroupInputTable.get_children()]
        GroupTypes = [self.GroupSettingsTable.item(item)['values'][1]
                      for item in self.GroupSettingsTable.get_children()]
        d_highBounds = [self.GroupSettingsTable.item(item)['values'][2]
                        for item in self.GroupSettingsTable.get_children()]
        d_LowBounds = [self.GroupSettingsTable.item(item)['values'][3]
                       for item in self.GroupSettingsTable.get_children()]
        b_highBounds = [self.GroupSettingsTable.item(item)['values'][4]
                        for item in self.GroupSettingsTable.get_children()]
        b_LowBounds = [self.GroupSettingsTable.item(item)['values'][5]
                       for item in self.GroupSettingsTable.get_children()]
        t_highBounds = [self.GroupSettingsTable.item(item)['values'][6]
                        for item in self.GroupSettingsTable.get_children()]
        t_LowBounds = [self.GroupSettingsTable.item(item)['values'][7]
                       for item in self.GroupSettingsTable.get_children()]

        costFunction = self.costFunctionStringBox.get()
        weightRun = self.dropdownWeight.get()
        reactionRun = self.dropdownReactions.get()
        internalForcesRun = self.dropdownInternalForces.get()

        lowerBounds = []
        upperBounds = []
        for i in range(len(GroupTypes)):
            lowerBounds.append(d_LowBounds)
            upperBounds.append(d_highBounds)
            if GroupTypes[i] == "Angle" or GroupTypes[i] == "RectHSS":
                lowerBounds.append(b_LowBounds)
                upperBounds.append(b_highBounds)
            lowerBounds.append(t_LowBounds)
            upperBounds.append(t_highBounds)

        self.mainwindow.GlobalOptimization(self, GroupAssignments, GroupTypes, lowerBounds, upperBounds, costFunction,
                                           weightRun, reactionRun, internalForcesRun)

    def cancelPopUp(self):
        """
        Removes the pop-up window.
        """

        self.top.destroy()

def getNonSetMembersData(d):
    """
    Gets the number and indeces of what member cross-sections are not set and need to be optimized.

    :param d: Info that defines the frame.
    :return: Number of non-set members. Indeces of what member are not set.
    """

    num = 0
    nonIndices = []
    for i in range(len(d.members[0])):
        if not d.members[3]:
            num += 1
            nonIndices.append(i)

    return num, nonIndices
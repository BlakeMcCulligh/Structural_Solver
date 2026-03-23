import tkinter as tk
from tkinter import ttk


class NewStructurePopUp:
    def __init__(self, root):
        self.root = root
        self.top = tk.Toplevel(root)
        width = 500
        height = 170
        self.centerWindow(width, height)
        self.top.title("Set Up New Structure")  # Set the title
        self.top.resizable(False, False)

        # Add widgets to the Toplevel window

        title = tk.Label(self.top, text="Set Up New Structure", font=('Helvetica', 20))
        wTitle = title.winfo_reqwidth()
        hTitle = title.winfo_reqheight()
        title.place(x=width/2-wTitle/2, y=0)


        DLabel = tk.Label(self.top, text="Dimentions: ", font=('Helvetica', 14))
        DLabel.place(x=0, y=hTitle+10)
        wDLabel = DLabel.winfo_reqwidth()
        hDLabel = DLabel.winfo_reqheight()

        DimentionsOPTIONS = ["2D", "3D"]
        self.Dimentions = tk.StringVar(self.top)
        self.Dimentions.set(DimentionsOPTIONS[0])
        Dimentions_menu = tk.OptionMenu(self.top, self.Dimentions, *DimentionsOPTIONS)
        Dimentions_menu.place(x=wDLabel+5, y=hTitle+10)

        TLabel = tk.Label(self.top, text="Type: ", font=('Helvetica', 14))
        TLabel.place(x=0, y=hTitle + hDLabel + 20)
        wTLabel = TLabel.winfo_reqwidth()

        TypeOPTIONS = ["Truss", "Frame"]
        self.Type = tk.StringVar(self.top)
        self.Type.set(TypeOPTIONS[0])
        Type_menu = tk.OptionMenu(self.top, self.Type, *TypeOPTIONS)
        Type_menu.place(x=wTLabel+5, y=hTitle + hDLabel + 20)

        tk.Button(self.top, text="Okay", command=self.endPopUp).place(x=420, y=140)
        tk.Button(self.top, text="Cancel", command=self.cancelPopUp).place(x=340, y=140)

        # Optional: ensure the pop-up stays on top of the main window
        self.top.grab_set()  # Prevents interaction with the main window while the pop-up is open
        root.wait_window(self.top)  # Waits for the 'top' window to be destroyed before proceeding in mainloop

    def centerWindow(self, width, height):

        mainwindow_x = self.root.winfo_x()
        mainwindow_y = self.root.winfo_y()
        mainwindow_width = self.root.winfo_width()
        mainwindow_height = self.root.winfo_height()

        x = (mainwindow_width // 2) + mainwindow_x - (width // 2)
        y = (mainwindow_height // 2) + mainwindow_y - (height // 2)

        self.top.geometry(f"{width}x{height}+{x}+{y}")

    def endPopUp(self):
        self.top.destroy()
        print(self.Dimentions.get())
        print(self.Type.get())
        #todo

    def cancelPopUp(self):
        self.top.destroy()

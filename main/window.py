
import tkinter as tk

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
        menubar.add_cascade(label="New", menu=new_menu)


        open_menu = tk.Menu(menubar, tearoff="off")
        open_menu.add_command(label='Browse Files', command=self.browseFilesToOpen)
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


    def browseFilesToOpen(self):
        print("Long Term Goals")
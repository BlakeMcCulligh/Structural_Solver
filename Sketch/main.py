import tkinter as tk
from Sketch.constraints.constraints import Constraints
from Sketch.geometry import Geometry
from Sketch.gui.gui import GUI
from Sketch.solver.solver import Solver

solver = None
gui = None

def geometry_changed_by_GUI(active_point):
    global solver
    solver.solve(active_point)

def constraints_changed_by_GUI():
    global solver
    solver.solve(None)

def geometry_changed_by_solver():
    global gui
    global solver
    gui.degrees_of_freedom = solver.degrees_of_freedom
    gui.redraw_geometry()

def sketch():
    global solver
    global gui

    geometry = Geometry()
    constraints = Constraints()

    solver = Solver(geometry, geometry_changed_by_solver, constraints)

    root_widget = tk.Tk()
    root_widget.title('2D Geometric Constraint Solver')

    gui = GUI(root_widget, geometry, geometry_changed_by_GUI, constraints, constraints_changed_by_GUI, mainWindow = None, objective = None)
    gui.pack(fill="both", expand=True)

    #root_widget.mainloop()

    return geometry

def startSketch(root_widget, mainWindow, objective, nextFunction):
    global solver
    global gui

    geometry = Geometry()
    constraints = Constraints()

    solver = Solver(geometry, geometry_changed_by_solver, constraints)

    gui = GUI(root_widget, geometry, geometry_changed_by_GUI, constraints, constraints_changed_by_GUI, mainWindow, objective, nextFunction)
    gui.pack(fill="both", expand=True)

    return geometry

#geometry = sketch()

import tkinter as tk

from Sketch import main
from main.window import MainWindow

if __name__ == '__main__':

    root_widget = tk.Tk()
    root_widget.title('Structural Solver')

    MainWindow(root_widget, main.sketch)

    root_widget.mainloop()
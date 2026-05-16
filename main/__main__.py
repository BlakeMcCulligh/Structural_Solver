"""
Starts the program.
"""
#test
import tkinter as tk

from main.window import MainWindow

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

if __name__ == '__main__':

    root_widget = tk.Tk()
    root_widget.title('Structural Solver')

    MainWindow(root_widget)

    root_widget.mainloop()
"""
Main of the program.
Handels starting the program and holds the main class and window.
"""

import tkinter as tk
from concurrent.futures import ThreadPoolExecutor

from rebuild_program_layout.frame_3d_frame import Frame3DFrame
from rebuild_program_layout.frame_main import FrameMain

__author__ = "Blake McCulligh"
__copyright__ = ""
__credits__ = ["Blake McCulligh"]

__license__ = ""
__version__ = ""
__maintainer__ = "Blake McCulligh"
__email__ = "bmcculli@uwaterloo.ca"
__status__ = ""

class Structural_Solver:
    """
    Object containing the program.
    """

    def __init__(self, root: tk.Tk, executor: ThreadPoolExecutor):
        """
        Constructor: Starts the program.

        :param root: Root of the main window.
        :type root: tkinter.Tk
        :param executor: Executor for running calculations in a seprate thread.
        :type executor:ThreadPoolExecutor
        """

        self.root = root
        self.root.geometry("400x300")

        self.main_frame = FrameMain(parent=self.root, controller=self)
        self.frame_3d_f = Frame3DFrame(parent=self.root, controller=self)

        self.current_frame = self.main_frame
        self.current_frame.pack(fill="both", expand=True)

        self._center_window()

        self.current_frame.on_show()

        self.executor = executor

    def _center_window(self) -> None:
        """
        Centers the window on the screen and sets its dimensions.
        """

        WIDTH = 1000
        HEIGHT = 800

        screen_width = self.root.winfo_screenwidth()
        screen_height = self.root.winfo_screenheight()

        if WIDTH > screen_width or HEIGHT > screen_height:
            WIDTH = screen_width * 0.8
            HEIGHT = screen_height * 0.8
        x = (screen_width // 2) - (WIDTH // 2)
        y = (screen_height // 2) - (HEIGHT // 2) - 50

        self.root.geometry(f"{WIDTH}x{HEIGHT}+{x}+{y}")

    def switch_to_frame_3d_frame(self, file_path: str | None = None) -> None:
        """
        Switches from the opening display frame, to the dispaly frame for a 3D frame structure.

        :param file_path: File path to the 3D frame structure being opened. Optional.
        :type file_path: str or None
        """

        self.current_frame.pack_forget()
        self.current_frame = self.frame_3d_f
        self.current_frame.pack(fill="both", expand=True)
        self.current_frame.on_show()

        if file_path is not None:
            self.current_frame.frame_3d.open_frame(file_path)
            self.current_frame.frame_3d.open_results(file_path.split('.')[0])

if __name__ == "__main__":
    r = tk.Tk()

    e = ThreadPoolExecutor(max_workers=1)

    Structural_Solver(r, e)

    r.mainloop()
import tkinter as tk
from tkinter import ttk
import multiprocessing
import time


# 1. Define the task that runs in an independent process
def heavy_calculation_process(process_id, data_queue):
    """This runs entirely in its own process memory space."""
    time.sleep(3)  # Simulate a 3-second heavy computational task
    result = f"Process {process_id} finished successfully!"
    data_queue.put(result)  # Safely send the result back to the main process


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Multiprocessing with Tkinter")
        self.geometry("400x250")

        # Create a thread/process safe queue
        self.result_queue = multiprocessing.Queue()
        self.process_counter = 0

        # Setup simple UI elements
        self.status_label = ttk.Label(self, text="Click below to start independent processes.")
        self.status_label.pack(pady=20)

        self.start_btn = ttk.Button(self, text="Launch New Process", command=self.start_new_process)
        self.start_btn.pack(pady=10)

        self.listbox = tk.Listbox(self)
        self.listbox.pack(pady=10, fill="both", expand=True)

        # 2. Begin the polling loop immediately
        self.check_queue_loop()

    def start_new_process(self):
        self.process_counter += 1
        self.listbox.insert(tk.END, f"Started Process #{self.process_counter}...")

        # 3. Create and initialize the separate process
        new_process = multiprocessing.Process(
            target=heavy_calculation_process,
            args=(self.process_counter, self.result_queue)
        )
        new_process.daemon = True  # Ensures process dies if main GUI window is closed
        new_process.start()

    def check_queue_loop(self):
        """Checks the queue for new messages without freezing the GUI."""
        try:
            # Look for available data without blocking
            while True:
                result_message = self.result_queue.get_nowait()
                self.listbox.insert(tk.END, result_message)
        except multiprocessing.queues.Empty:
            # The queue is currently empty; do nothing
            pass

        # 4. Schedule this function to run again in 100 milliseconds
        self.after(100, self.check_queue_loop)


# 5. Guard statement is MANDATORY for Python multiprocessing to prevent recursive loops
if __name__ == "__main__":
    app = App()
    app.mainloop()
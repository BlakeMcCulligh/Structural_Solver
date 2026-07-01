import tkinter as tk
from concurrent.futures import ThreadPoolExecutor
from scipy.optimize import differential_evolution


def objective_function(x):
    # Your expensive computations here
    return x[0] ** 2 + x[1] ** 2


def run_optimizer_in_thread():
    bounds = [(-5, 5), (-5, 5)]

    # Run SciPy optimizer in the background (uses all available cores)
    result = differential_evolution(objective_function, bounds, workers=-1)

    # Safely update your Tkinter UI with the results on the main thread
    root.after(0, update_ui, result.x, result.fun)


def start_optimization():
    # Submit the optimizer to a background thread to prevent UI lockup
    executor.submit(run_optimizer_in_thread)


def update_ui(optimal_x, optimal_val):
    result_label.config(text=f"Success! Best X: {optimal_x}, Value: {optimal_val:.4f}")


# Initialize Tkinter GUI
if __name__ == "__main__":
    root = tk.Tk()
    tk.Button(root, text="Start Optimization", command=start_optimization).pack(pady=20)
    result_label = tk.Label(root, text="Idle")
    result_label.pack(pady=10)

    # Start the Thread Pool Executor
    executor = ThreadPoolExecutor(max_workers=1)

    root.mainloop()
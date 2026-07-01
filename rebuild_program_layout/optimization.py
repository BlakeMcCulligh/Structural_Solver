"""
Handels the optimization of the structure.
"""

from scipy import optimize

# TODO

def objective_function(x):

    return x[0] ** 2 + x[1] ** 2

def run_optimizer(window_root, bounds):

    # Run SciPy optimizer in the background (uses all available cores)
    result = optimize.differential_evolution(objective_function, bounds, updating='deferred', workers=-1)

    # update your Tkinter UI with the results on the main thread
    window_root.after(0, return_optimization_results, result)

def start_optimization(executor, window_root, bounds):
    # run optimizer in a background thread
    executor.submit(lambda: run_optimizer(window_root, bounds))

def return_optimization_results(result):
    print("resalts outside: ", result)
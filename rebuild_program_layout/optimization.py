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
    print("resalts outside 2: ", result)

def global_optimzation_results(results):
    print("Global Optimization Results: ", results)

def start_global_optimization(executor, window_root, frame, group_assignments, group_types, lower_bounds, upper_bounds,
                              cost_function, weight_run, reaction_run, internal_forces_run):
    print("Call Frame Optimzer")
    future = executor.submit(lambda: frame.optimize(window_root.after, group_assignments, group_types, lower_bounds, upper_bounds,
                 cost_function, weight_run, reaction_run, internal_forces_run))

    # Wait for the thread to finish processing
    future.done()

    error = future.exception()
    if error:
        print(f"The thread failed with error type: {type(error).__name__}")
        print(f"Error message: {error}")

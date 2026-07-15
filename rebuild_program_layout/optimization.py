"""
Handels the optimization of the structure.
"""

# TODO

def start_global_optimization(executor, frame3D, solver_frame, group_assignments, group_types, lower_bounds,
                              upper_bounds,cost_function, weight_run, reaction_run, internal_forces_run):

    global_opt = executor.submit(lambda: solver_frame.optimize(group_assignments, group_types, lower_bounds, upper_bounds,
                 cost_function, weight_run, reaction_run, internal_forces_run))

    # Wait for the thread to finish processing
    global_opt.done()
    frame3D.handel_global_optimization_results(global_opt.result())

    error = global_opt.exception()
    if error:
        print(f"The thread failed with error type: {type(error).__name__}")
        print(f"Error message: {error}")

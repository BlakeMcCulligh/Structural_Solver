import cvxpy as cp
import numpy as np

# Data
volume_bound = 5
l = np.array([5, 6, 7, 3, 5], dtype=float)
m = len(l)

d = 2 * np.sum([1, 0, 0, 0, 1])

B = np.array([
    [1, 2],
    [2, 3],
    [3, 4],
    [4, 5],
    [5, 1]
], dtype=float)

E = np.diag([5, 5, 5, 5, 5])

# Variables
w = cp.Variable()
a = cp.Variable(m)

# f must be length 2 to match B'
f = cp.Variable(2)

# Diagonal term
D = cp.diag(cp.multiply(l, a))

# Lower-right block
K = B.T @ E @ D @ B   # (2×2)

# SDP block matrix
M = cp.bmat([
    [cp.reshape(w, (1, 1)), cp.reshape(f, (1, 2))],
    [cp.reshape(f, (2, 1)), K]
])

# Constraints
constraints = [
    M >> 0,
    l @ a <= volume_bound,
    a >= 0
]

# Problem
prob = cp.Problem(cp.Minimize(w), constraints)
prob.solve(solver=cp.SCS)

print("Status:", prob.status)
print("w:", w.value)
print("a:", a.value)

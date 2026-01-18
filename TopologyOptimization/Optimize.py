import cvxpy as cp
import numpy as np


class TrussTopologyOptimization:
    def __init__(self, volumeBound, members, nodes, memberLengths, E):
        self.volumeBound = volumeBound
        self.members = members
        self.nodes = nodes
        self.l = memberLengths
        self.m = len(self.l)
        self.Ediagenal = np.diag([E] * self.m)

        self.B = self.build_B()
        self.w = None
        self.A = None
        self.prob = None

    def build_B(self):
        n_nodes = len(self.nodes)
        n_members = len(self.members)
        d = 2 * n_nodes

        B = np.zeros((n_members, d))

        for k, (i, j) in enumerate(self.members):
            xi, yi = self.nodes[i]
            xj, yj = self.nodes[j]

            dx = xj - xi
            dy = yj - yi
            L = np.sqrt(dx ** 2 + dy ** 2)

            c = dx / L
            s = dy / L

            B[k, 2 * i] = -c / L
            B[k, 2 * i + 1] = -s / L
            B[k, 2 * j] = c / L
            B[k, 2 * j + 1] = s / L

        return B

    def optimize(self):
        self.w = cp.Variable()
        self.A = cp.Variable(self.m)

        # f must be length 2 to match B'
        d = self.B.shape[1]
        f = cp.Variable(d)

        # Diagonal term
        D = cp.diag(cp.multiply(self.l, self.A))

        # Lower-right block
        K = self.B.T @ self.Ediagenal @ D @ self.B  # (2×2)

        # SDP block matrix
        M = cp.bmat([
            [cp.reshape(self.w, (1, 1), order="F"), cp.reshape(f, (1, d), order="F")],
            [cp.reshape(f, (d, 1), order="F"), K]
        ])

        # Constraints
        constraints = [
            M >> 0,
            self.l @ self.A <= self.volumeBound,
            self.A >= 0
        ]

        # Problem
        self.prob = cp.Problem(cp.Minimize(self.w), constraints)
        self.prob.solve(solver=cp.SCS)

    def printResults(self):
        print("Status:", self.prob.status)
        print("w:", self.w.value)
        print("a:", self.A.value)



# volume_bound = 5
# B = np.array([
#     [1, 2],
#     [2, 3],
#     [3, 4],
#     [4, 5],
#     [5, 1]
# ], dtype=float)
# l = np.array([5, 6, 7, 3, 5], dtype=float)
# E = 5
#
# opt = TrussTopologyOptimization(volume_bound, B, [1], l, E)
# opt.optimize()
# opt.printResults()


import cvxpy as cp
import numpy as np


class TrussTopologyOptimization:
    def __init__(self, volumeBound, members, nodes, memberLengths, E, supports, loads):
        self.volumeBound = volumeBound
        self.members = members
        self.nodes = nodes
        self.l = memberLengths
        self.m = len(self.l)
        self.E = E

        self.w = None
        self.A = None
        self.prob = None

        # node_id : (fix_x, fix_y)
        self.supports = supports

        # node_id : (Fx, Fy)
        self.loads = loads

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

    def get_dof_sets(self):
        fixed = []
        for node, (fix_x, fix_y) in self.supports.items():
            if fix_x:
                fixed.append(2 * node)
            if fix_y:
                fixed.append(2 * node + 1)

        all_dofs = set(range(2 * len(self.nodes)))
        free = sorted(all_dofs - set(fixed))

        return free, sorted(fixed)

    def build_force_vector(self):
        F = np.zeros(2 * len(self.nodes))
        for node, (Fx, Fy) in self.loads.items():
            F[2 * node] = Fx
            F[2 * node + 1] = Fy
        return F

    def optimize(self):
        B = self.build_B()

        free_dofs, fixed_dofs = self.get_dof_sets()

        F_full = self.build_force_vector()
        F_free = F_full[free_dofs]

        B_free = B[:, free_dofs]

        d = B_free.shape[1]

        self.A = cp.Variable(len(self.members), nonneg=True)
        self.w = cp.Variable()
        f = cp.Variable(d)

        D = cp.diag(cp.multiply(self.E * self.l, self.A))
        K = B_free.T @ D @ B_free

        # SDP block matrix
        M = cp.bmat([
            [cp.reshape(self.w, (1, 1), order="F"), cp.reshape(f, (1, d), order="F")],
            [cp.reshape(f, (d, 1), order="F"), K]
        ])

        # Constraints
        constraints = [
            M >> 0,
            self.l @ self.A <= self.volumeBound,
            self.A >= 0,
            f == F_free
        ]

        # Problem
        self.prob = cp.Problem(cp.Minimize(self.w), constraints)
        self.prob.solve(solver=cp.SCS,verbose=True)

    def printResults(self):
        print("Status:", self.prob.status)
        print("w:", self.w.value)
        print("a:", self.A.value)
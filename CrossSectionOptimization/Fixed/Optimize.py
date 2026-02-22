import numpy as np
from scipy.optimize import minimize

from CrossSectionOptimization.Fixed.CrossSectionCatalog import Catalog, CrossSection
#from CrossSectionOptimization.Fixed.SolveFrame import FrameSolver
from CrossSectionOptimization.Fixed.Structure import structure



# 2D frame
def beam_stiffness(E, A, I, L):
    k = np.zeros((6, 6))

    k[0, 0] = k[3, 3] = E * A / L
    k[0, 3] = k[3, 0] = -E * A / L

    k[1, 1] = k[4, 4] = 12 * E * I / L**3
    k[1, 4] = k[4, 1] = -12 * E * I / L**3

    k[1, 2] = k[2, 1] = 6 * E * I / L**2
    k[1, 5] = k[5, 1] = 6 * E * I / L**2
    k[4, 2] = k[2, 4] = -6 * E * I / L**2
    k[4, 5] = k[5, 4] = -6 * E * I / L**2

    k[2, 2] = k[5, 5] = 4 * E * I / L
    k[2, 5] = k[5, 2] = 2 * E * I / L

    return k

def rotation_matrix(theta):
    c = np.cos(theta)
    s = np.sin(theta)

    R = np.zeros((6, 6))
    R[:3, :3] = [[c, s, 0], [-s, c, 0], [0, 0, 1]]
    R[3:, 3:] = R[:3, :3]
    return R

class FrameSolver:
    def __init__(self, model, catalog, E=210e9):
        self.model = model
        self.catalog = catalog
        self.E = E

    def assemble(self, x, p):
        n_dof = 3 * len(self.model.nodes)
        K = np.zeros((n_dof, n_dof))

        for j, member in enumerate(self.model.members):
            n1, n2 = member.n1, member.n2
            L = member.length(self.model.nodes)
            theta = member.angle(self.model.nodes)

            K0 = np.eye(6) * 1e-9

            Ke = K0.copy()
            for i, sec in enumerate(self.catalog.sections):
                w = DMO.simp(x[j, i], p)
                k_i = beam_stiffness(self.E, sec.A, sec.I, L)
                Ke += w * (k_i - K0)

            R = rotation_matrix(theta)
            Ke = R.T @ Ke @ R

            dofs = [
                3*n1, 3*n1+1, 3*n1+2,
                3*n2, 3*n2+1, 3*n2+2
            ]

            for a in range(6):
                for b in range(6):
                    K[dofs[a], dofs[b]] += Ke[a, b]

        return K

    def solve(self, x, p):
        K = self.assemble(x, p)
        n_dof = K.shape[0]
        F = np.zeros(n_dof)

        for node, load in self.model.loads.items():
            F[3*node:3*node+3] += load

        fixed = []
        for i, node in enumerate(self.model.nodes):
            for d, f in enumerate(node.fixed):
                if f:
                    fixed.append(3*i + d)

        free = np.setdiff1d(np.arange(n_dof), fixed)

        U = np.zeros(n_dof)
        U[free] = np.linalg.solve(K[np.ix_(free, free)], F[free])

        return K, U, F












class DMO:
    @staticmethod
    def simp(x, p):
        return x**p

    @staticmethod
    def simp_deriv(x, p):
        return p * x**(p - 1)

    @staticmethod
    def ramp(x, r):
        return x / (1 + r * (1 - x))

def mass(x, model, catalog, seq):
    m = 0.0
    for j, member in enumerate(model.members):
        L = member.length(model.nodes)
        for i, sec in enumerate(catalog.sections):
            w = DMO.ramp(x[j, i], seq)
            m += w * sec.rho * sec.A * L
    return m

class FrameOptimizer:
    def __init__(self, model, catalog, alpha=0.3):
        self.model = model
        self.catalog = catalog
        self.alpha = alpha
        self.fem = FrameSolver(model, catalog)

    def objective(self, x_flat):
        x = x_flat.reshape(len(self.model.members), self.catalog.n)
        K, U, F = self.fem.solve(x, p=2)
        C = F @ U
        m = mass(x, self.model, self.catalog, seq=2)
        return self.alpha * C + (1 - self.alpha) * m

    def optimize(self):
        n_m = len(self.model.members)
        n_c = self.catalog.n

        x0 = np.ones(n_m * n_c) / n_c

        cons = {
            "type": "eq",
            "fun": lambda x: np.sum(
                x.reshape(n_m, n_c), axis=1
            ) - 1.0
        }

        res = minimize(
            self.objective,
            x0,
            method="SLSQP",
            constraints=[cons],
            options={"ftol": 1e-6, "disp": True}
        )

        return res

# testing
model = structure()

n0 = model.add_node(0, 0, fixed=(True, True, True))
n1 = model.add_node(4, 0)
n2 = model.add_node(4, 3)

model.add_member(n0, n1)
model.add_member(n1, n2)

model.add_load(n2, fy=-10_000)

catalog = Catalog([
    CrossSection(A=0.002, I=8e-6, rho=7850, sigma_y=235e6),
    CrossSection(A=0.004, I=3e-5, rho=7850, sigma_y=235e6),
])

opt = FrameOptimizer(model, catalog, alpha=0.4)
result = opt.optimize()

print(result.x.reshape(len(model.members), catalog.n))

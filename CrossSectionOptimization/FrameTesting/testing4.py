import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import minimize


def run_large_scale_2D_frame_DMO(
    nodes,
    members,
    candidates,
    load_cases,
    fixed_dofs,
    load_weights=None,
    E=210e9,
    rho=7850,
    sigma_y=235e6,
    alpha=0.3,
    lambda_min_required=1.0,
    continuation_steps=5
):
    """
    Large-scale 2D FrameTesting DMO Optimizer
    Efficient for 100–300+ members
    """

    n_nodes = len(nodes)
    n_members = len(members)
    n_cand = len(candidates)
    total_dof = 3 * n_nodes

    if load_weights is None:
        load_weights = np.ones(len(load_cases)) / len(load_cases)

    free_dofs = np.setdiff1d(np.arange(total_dof), fixed_dofs)

    epsilon = 1e-8

    # ==========================================================
    # PRECOMPUTE ELEMENT MATRICES
    # ==========================================================

    elem_dofs = []
    elem_lengths = []
    elem_k = []

    def beam_stiff(E, A, I, x1, y1, x2, y2):
        L = np.hypot(x2 - x1, y2 - y1)
        c = (x2 - x1) / L
        s = (y2 - y1) / L

        k_local = np.array([
            [A*E/L,0,0,-A*E/L,0,0],
            [0,12*E*I/L**3,6*E*I/L**2,0,-12*E*I/L**3,6*E*I/L**2],
            [0,6*E*I/L**2,4*E*I/L,0,-6*E*I/L**2,2*E*I/L],
            [-A*E/L,0,0,A*E/L,0,0],
            [0,-12*E*I/L**3,-6*E*I/L**2,0,12*E*I/L**3,-6*E*I/L**2],
            [0,6*E*I/L**2,2*E*I/L,0,-6*E*I/L**2,4*E*I/L]
        ])

        T = np.array([
            [c,s,0,0,0,0],
            [-s,c,0,0,0,0],
            [0,0,1,0,0,0],
            [0,0,0,c,s,0],
            [0,0,0,-s,c,0],
            [0,0,0,0,0,1]
        ])

        return T.T @ k_local @ T, L

    for (n1, n2) in members:
        x1, y1 = nodes[n1]
        x2, y2 = nodes[n2]

        dofs = np.array([
            3*n1, 3*n1+1, 3*n1+2,
            3*n2, 3*n2+1, 3*n2+2
        ])

        elem_dofs.append(dofs)

        k_list = []
        for cand in candidates:
            k, L = beam_stiff(E, cand["A"], cand["I"], x1, y1, x2, y2)
            k_list.append(k)

        elem_k.append(k_list)
        elem_lengths.append(L)

    elem_lengths = np.array(elem_lengths)

    # ==========================================================
    # INTERPOLATION
    # ==========================================================

    def wc(x, p): return x**p
    def wd(x, s): return x / (1 + s*(1-x))

    # ==========================================================
    # ASSEMBLY (SPARSE)
    # ==========================================================

    def assemble_K(x, p):

        x = x.reshape(n_members, n_cand)

        rows = []
        cols = []
        data = []

        for e in range(n_members):

            dofs = elem_dofs[e]
            member_sum = np.sum(x[e])

            Ke = np.zeros((6,6))

            for i in range(n_cand):
                Ke += wc(x[e,i], p) * elem_k[e][i]

            Ke += epsilon * (1 - member_sum) * np.eye(6)

            for i in range(6):
                for j in range(6):
                    rows.append(dofs[i])
                    cols.append(dofs[j])
                    data.append(Ke[i,j])

        K = sp.coo_matrix((data,(rows,cols)),
                          shape=(total_dof,total_dof)).tocsr()

        return K

    # ==========================================================
    # SOLVE ALL LOAD CASES
    # ==========================================================

    def solve_all(x, p):

        K = assemble_K(x, p)
        Kff = K[free_dofs][:, free_dofs]

        displacements = []

        for F in load_cases:
            uf = spla.spsolve(Kff, F[free_dofs])
            u = np.zeros(total_dof)
            u[free_dofs] = uf
            displacements.append(u)

        return displacements, K

    # ==========================================================
    # STABILITY (SMALLEST EIGENVALUE)
    # ==========================================================

    def stability_constraint(x, p):

        u_all, K = solve_all(x, p)
        Kff = K[free_dofs][:, free_dofs]

        worst_margin = 1e9

        for u in u_all:

            Kg = sp.lil_matrix((total_dof,total_dof))

            for e in range(n_members):
                dofs = elem_dofs[e]
                ue = u[dofs]
                L = elem_lengths[e]

                axial = (ue[3] - ue[0]) / L
                N = E * axial

                k_geo = (N/(30*L))*np.array([
                    [0,0,0,0,0,0],
                    [0,36,3*L,0,-36,3*L],
                    [0,3*L,4*L**2,0,-3*L,-L**2],
                    [0,0,0,0,0,0],
                    [0,-36,-3*L,0,36,-3*L],
                    [0,3*L,-L**2,0,-3*L,4*L**2]
                ])

                for i in range(6):
                    for j in range(6):
                        Kg[dofs[i],dofs[j]] += k_geo[i,j]

            Kg = Kg.tocsr()
            Kgff = Kg[free_dofs][:,free_dofs]

            # smallest eigenvalue only
            eig = spla.eigs(Kff, k=1, M=Kgff, which='SM',
                            return_eigenvectors=False)
            lam = np.real(eig[0])

            worst_margin = min(worst_margin, lam - lambda_min_required)

        return worst_margin

    # ==========================================================
    # OBJECTIVE
    # ==========================================================

    def objective(x, p, s):

        u_all, _ = solve_all(x, p)

        compliance = 0
        for w, u, F in zip(load_weights, u_all, load_cases):
            compliance += w * (F @ u)

        xmat = x.reshape(n_members, n_cand)

        mass = np.sum(
            wd(xmat, s) *
            np.array([c["A"] for c in candidates]) *
            rho *
            elem_lengths[:,None]
        )

        return alpha * compliance + (1-alpha) * mass

    # ==========================================================
    # OPTIMIZATION LOOP
    # ==========================================================

    x0 = np.ones((n_members,n_cand)) / n_cand
    x_flat = x0.flatten()
    bounds = [(0,1)] * len(x_flat)

    for cont in range(continuation_steps):

        p = 1 + cont
        s = cont

        constraints = []

        for e in range(n_members):
            idx = slice(e*n_cand,(e+1)*n_cand)
            constraints.append({
                "type":"ineq",
                "fun":lambda x,idx=idx: 1 - np.sum(x[idx])
            })

        constraints.append({
            "type":"ineq",
            "fun":lambda x,p=p: stability_constraint(x,p)
        })

        res = minimize(objective, x_flat,
                       args=(p,s),
                       method="SLSQP",
                       bounds=bounds,
                       constraints=constraints,
                       options={"maxiter":200})

        x_flat = res.x

    return {
        "weights": x_flat.reshape(n_members,n_cand),
        "selected": np.argmax(x_flat.reshape(n_members,n_cand),axis=1)
    }
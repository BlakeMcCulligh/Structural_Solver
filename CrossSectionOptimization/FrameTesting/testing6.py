import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import minimize
import matplotlib.pyplot as plt


# ==========================================================
# OPTIMIZER (Efficient 2D DMO)
# ==========================================================

def run_large_scale_2D_frame_DMO(
    nodes,
    members,
    candidates,
    load_cases,
    fixed_dofs,
    E=210e9,
    rho=7850,
    alpha=0.3,
    continuation_steps=3
):

    n_nodes = len(nodes)
    n_members = len(members)
    n_cand = len(candidates)
    total_dof = 3*n_nodes
    free = np.setdiff1d(np.arange(total_dof), fixed_dofs)

    epsilon = 1e-8

    elem_dofs = []
    elem_lengths = []
    elem_k = []

    def beam_stiff(A,I,x1,y1,x2,y2):
        L = np.hypot(x2-x1,y2-y1)
        c = (x2-x1)/L
        s = (y2-y1)/L

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

    for (n1,n2) in members:
        x1,y1 = nodes[n1]
        x2,y2 = nodes[n2]
        dofs = np.array([3*n1,3*n1+1,3*n1+2,
                         3*n2,3*n2+1,3*n2+2])
        elem_dofs.append(dofs)

        k_list = []
        for cand in candidates:
            k,L = beam_stiff(cand["A"],cand["I"],x1,y1,x2,y2)
            k_list.append(k)

        elem_k.append(k_list)
        elem_lengths.append(L)

    elem_lengths = np.array(elem_lengths)

    def wc(x,p): return x**p
    def wd(x,s): return x/(1+s*(1-x))

    def assemble_K(x,p):
        x = x.reshape(n_members,n_cand)
        K = sp.lil_matrix((total_dof,total_dof))

        for e in range(n_members):
            Ke = np.zeros((6,6))
            for i in range(n_cand):
                Ke += wc(x[e,i],p)*elem_k[e][i]
            dofs = elem_dofs[e]
            for i in range(6):
                for j in range(6):
                    K[dofs[i],dofs[j]] += Ke[i,j]

        return K.tocsr()

    def solve(x,p):
        K = assemble_K(x,p)
        Kff = K[free][:,free]
        F = load_cases[0]
        uf = spla.spsolve(Kff,F[free])
        u = np.zeros(total_dof)
        u[free] = uf
        return u

    def objective(x,p,s):
        u = solve(x,p)
        compliance = load_cases[0] @ u
        xmat = x.reshape(n_members,n_cand)
        A_vec = np.array([c["A"] for c in candidates])
        mass = np.sum(wd(xmat,s)*A_vec*rho*elem_lengths[:,None])
        return alpha*compliance + (1-alpha)*mass

    x0 = np.ones((n_members,n_cand))/n_cand
    x_flat = x0.flatten()
    bounds = [(0,1)]*len(x_flat)

    for cont in range(continuation_steps):
        p = 1+cont
        s = cont
        res = minimize(objective,x_flat,
                       args=(p,s),
                       method="SLSQP",
                       bounds=bounds,
                       options={"maxiter":100})
        x_flat = res.x

    x_opt = x_flat.reshape(n_members,n_cand)
    selected = np.argmax(x_opt,axis=1)
    u_final = solve(x_flat,p)

    return x_opt, selected, u_final


# ==========================================================
# PLOTTING FUNCTION
# ==========================================================

def plot_structure(nodes, members, selected, u, scale=1.0, title=""):

    plt.figure(figsize=(8,6))

    nodes_def = nodes.copy()
    nodes_def[:,0] = nodes_def[:,0] +  u[0::3]*scale
    nodes_def[:,1] = nodes_def[:,1] +  u[1::3]*scale

    colors = ['blue','red','green','purple','orange']

    for i,(n1,n2) in enumerate(members):
        c = colors[selected[i] % len(colors)]

        # undeformed
        plt.plot([nodes[n1,0],nodes[n2,0]],
                 [nodes[n1,1],nodes[n2,1]],
                 '--',color='gray',alpha=0.5)

        # deformed
        plt.plot([nodes_def[n1,0],nodes_def[n2,0]],
                 [nodes_def[n1,1],nodes_def[n2,1]],
                 color=c,linewidth=2)

    plt.scatter(nodes[:,0],nodes[:,1],color='black')
    plt.axis('equal')
    plt.title(title)
    plt.show()


# ==========================================================
# EXAMPLE PROBLEMS
# ==========================================================

def run_examples():

    candidates = [
        {"A":5e-4,"I":1e-6},
        {"A":1e-3,"I":4e-6}
    ]

    # Example 1: Cantilever
    nodes = np.array([[0,0],[3,0]])
    members=[(0,1)]
    F = np.zeros(6)
    F[4]=-10000

    x,sel,u = run_large_scale_2D_frame_DMO(
        nodes,members,candidates,[F],[0,1,2])

    plot_structure(nodes,members,sel,u,scale=100,
                   title="Cantilever Beam")

    # Example 2: Triangle
    nodes = np.array([[0,0],[4,0],[2,3]])
    members=[(0,1),(1,2),(0,2)]
    F = np.zeros(9)
    F[7]=-15000

    x,sel,u = run_large_scale_2D_frame_DMO(
        nodes,members,candidates,[F],[0,1,2])

    plot_structure(nodes,members,sel,u,scale=50,
                   title="Triangle FrameTesting")

    # Example 3: Portal FrameTesting
    nodes = np.array([[0,0],[4,0],[0,3],[4,3]])
    members=[(0,1),(2,3),(0,2),(1,3)]
    F = np.zeros(12)
    F[10]=-20000

    x,sel,u = run_large_scale_2D_frame_DMO(
        nodes,members,candidates,[F],[0,1,2])

    plot_structure(nodes,members,sel,u,scale=40,
                   title="Portal FrameTesting")


# ==========================================================
# RUN DEMOS
# ==========================================================

run_examples()
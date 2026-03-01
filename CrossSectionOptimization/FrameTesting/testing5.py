import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import minimize


def run_demo_problems():

    # ==========================================================
    # CORE OPTIMIZER
    # ==========================================================

    def run_large_scale_2D_frame_DMO(
        nodes,
        members,
        candidates,
        load_cases,
        fixed_dofs,
        load_weights=None,
        E=210e9,
        rho=7850,
        alpha=0.3,
        lambda_min_required=1.0,
        continuation_steps=4
    ):

        n_nodes = len(nodes)
        n_members = len(members)
        n_cand = len(candidates)
        total_dof = 3*n_nodes

        if load_weights is None:
            load_weights = np.ones(len(load_cases))/len(load_cases)

        free = np.setdiff1d(np.arange(total_dof), fixed_dofs)

        epsilon = 1e-8

        # ------------------------------------------------------
        # PRECOMPUTE ELEMENT STIFFNESS
        # ------------------------------------------------------

        elem_dofs = []
        elem_lengths = []
        elem_k = []

        def beam_stiff(A,I,x1,y1,x2,y2):

            L = np.hypot(x2-x1,y2-y1)
            if L < 1e-12:
                raise ValueError("Zero length element")

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

            dofs = np.array([
                3*n1,3*n1+1,3*n1+2,
                3*n2,3*n2+1,3*n2+2
            ])

            elem_dofs.append(dofs)

            k_list = []
            for cand in candidates:
                k,L = beam_stiff(cand["A"],cand["I"],x1,y1,x2,y2)
                k_list.append(k)

            elem_k.append(k_list)
            elem_lengths.append(L)

        elem_lengths = np.array(elem_lengths)

        # ------------------------------------------------------
        # INTERPOLATION
        # ------------------------------------------------------

        def wc(x,p): return x**p
        def wd(x,s): return x/(1+s*(1-x))

        # ------------------------------------------------------
        # ASSEMBLY
        # ------------------------------------------------------

        def assemble_K(x,p):

            x = x.reshape(n_members,n_cand)

            rows=[]
            cols=[]
            data=[]

            for e in range(n_members):

                dofs = elem_dofs[e]
                Ke = np.zeros((6,6))
                member_sum = np.sum(x[e])

                for i in range(n_cand):
                    Ke += wc(x[e,i],p)*elem_k[e][i]

                Ke += epsilon*(1-member_sum)*np.eye(6)

                for i in range(6):
                    for j in range(6):
                        rows.append(dofs[i])
                        cols.append(dofs[j])
                        data.append(Ke[i,j])

            return sp.coo_matrix((data,(rows,cols)),
                                 shape=(total_dof,total_dof)).tocsr()

        # ------------------------------------------------------
        # SOLVE
        # ------------------------------------------------------

        def solve_all(x,p):

            K = assemble_K(x,p)
            Kff = K[free][:,free]

            displacements=[]

            for F in load_cases:
                try:
                    uf = spla.spsolve(Kff,F[free])
                except:
                    return None,None
                u = np.zeros(total_dof)
                u[free]=uf
                displacements.append(u)

            return displacements,K

        # ------------------------------------------------------
        # STABILITY CHECK
        # ------------------------------------------------------

        def stability_constraint(x,p):

            u_all,K = solve_all(x,p)
            if u_all is None:
                return -1e3

            Kff = K[free][:,free]

            worst = 1e9

            for u in u_all:

                Kg = sp.lil_matrix((total_dof,total_dof))

                for e in range(n_members):

                    dofs = elem_dofs[e]
                    ue = u[dofs]
                    L = elem_lengths[e]

                    axial = (ue[3]-ue[0])/L
                    N = E*axial

                    if N >= 0:
                        continue  # no compression → no buckling

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
                Kgff = Kg[free][:,free]

                if Kgff.nnz == 0:
                    continue

                try:
                    eig = spla.eigs(Kff,k=1,M=Kgff,
                                    which='SM',
                                    return_eigenvectors=False)
                    lam = np.real(eig[0])
                except:
                    return -1e3

                worst = min(worst, lam-lambda_min_required)

            return worst

        # ------------------------------------------------------
        # OBJECTIVE
        # ------------------------------------------------------

        def objective(x,p,s):

            u_all,_ = solve_all(x,p)
            if u_all is None:
                return 1e12

            compliance=0
            for w,u,F in zip(load_weights,u_all,load_cases):
                compliance+=w*(F@u)

            xmat=x.reshape(n_members,n_cand)
            A_vec=np.array([c["A"] for c in candidates])

            mass=np.sum(
                wd(xmat,s)*A_vec*rho*elem_lengths[:,None]
            )

            return alpha*compliance+(1-alpha)*mass

        # ------------------------------------------------------
        # OPTIMIZATION
        # ------------------------------------------------------

        x0=np.ones((n_members,n_cand))/n_cand
        x_flat=x0.flatten()
        bounds=[(0,1)]*len(x_flat)

        for cont in range(continuation_steps):

            p=1+cont
            s=cont

            cons=[]

            for e in range(n_members):
                idx=slice(e*n_cand,(e+1)*n_cand)
                cons.append({
                    "type":"ineq",
                    "fun":lambda x,idx=idx:1-np.sum(x[idx])
                })

            cons.append({
                "type":"ineq",
                "fun":lambda x,p=p:stability_constraint(x,p)
            })

            res=minimize(objective,x_flat,
                         args=(p,s),
                         method="SLSQP",
                         bounds=bounds,
                         constraints=cons,
                         options={"maxiter":150})

            x_flat=res.x

        return x_flat.reshape(n_members,n_cand)

    # ==========================================================
    # EXAMPLE 1: Cantilever Beam
    # ==========================================================

    nodes=np.array([[0,0],[2,0]])
    members=[(0,1)]
    candidates=[
        {"A":5e-4,"I":1e-6},
        {"A":1e-3,"I":4e-6},
        {"A": 1e-2, "I": 8e-6},
        {"A": 1e-1, "I": 12e-6},
        {"A": 1, "I": 24e-6},

    ]
    F=np.zeros(6)
    F[4]=-10000

    result1=run_large_scale_2D_frame_DMO(
        nodes,members,candidates,
        [F],[0,1,2]
    )

    print("Example 1 (Cantilever) Selected:",
          np.argmax(result1,axis=1))
    print("Expected: largest section selected\n")

    # ==========================================================
    # EXAMPLE 2: Triangle Truss
    # ==========================================================

    nodes=np.array([[0,0],[3,0],[1.5,2]])
    members=[(0,1),(1,2),(0,2)]
    F=np.zeros(9)
    F[7]=-10000

    result2=run_large_scale_2D_frame_DMO(
        nodes,members,candidates,
        [F],[0,1,2]
    )

    print("Example 2 (Triangle) Selected:",
          np.argmax(result2,axis=1))
    print("Expected: compression member larger\n")

    # ==========================================================
    # EXAMPLE 3: 5-Bay FrameTesting (~10 members)
    # ==========================================================

    nodes=np.array([
        [0,0],[2,0],[4,0],
        [0,3],[2,3],[4,3]
    ])

    members=[
        (0,1),(1,2),
        (3,4),(4,5),
        (0,3),(1,4),(2,5),
        (1,3),(2,4)
    ]

    F=np.zeros(18)
    F[16]=-20000

    result3=run_large_scale_2D_frame_DMO(
        nodes,members,candidates,
        [F],[0,1,2]
    )

    print("Example 3 (FrameTesting) Selected:",
          np.argmax(result3,axis=1))
    print("Expected: columns larger than beams\n")


# ==============================================================
# RUN ALL DEMOS
# ==============================================================

run_demo_problems()
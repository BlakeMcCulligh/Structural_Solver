import numpy as np
from scipy.optimize import minimize


def run_full_dmo_frame_optimizer():

    # ==========================================================
    # PROBLEM SETUP
    # ==========================================================

    E = 210e9
    rho = 7850
    sigma_y = 235e6
    alpha = 0.3

    nodes = np.array([
        [0.0, 0.0],
        [3.0, 0.0],
        [3.0, 3.0]
    ])

    members = [
        (0, 1),
        (1, 2),
        (0, 2)
    ]

    candidates = [
        {"A": 5e-4, "I": 8e-6,  "c":0.05},
        {"A": 8e-4, "I": 2e-5,  "c":0.06},
        {"A": 1.2e-3,"I": 4e-5, "c":0.08}
    ]

    n_nodes = len(nodes)
    n_members = len(members)
    n_cand = len(candidates)

    total_dof = 3*n_nodes
    fixed_dofs = [0,1,2]
    free_dofs = list(set(range(total_dof))-set(fixed_dofs))

    F = np.zeros(total_dof)
    F[7] = -10000

    # ==========================================================
    # FEM
    # ==========================================================

    def beam_stiffness(E,A,I,x1,y1,x2,y2):
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

        return T.T @ k_local @ T, L, c, s

    K_elem = []
    lengths = []
    transforms = []

    for (n1,n2) in members:
        x1,y1 = nodes[n1]
        x2,y2 = nodes[n2]
        k_list = []
        k_base, L, c, s = beam_stiffness(E,1,1,x1,y1,x2,y2)
        for cand in candidates:
            k,_ ,_,_ = beam_stiffness(E,cand["A"],cand["I"],x1,y1,x2,y2)
            k_list.append(k)
        K_elem.append(k_list)
        lengths.append(L)
        transforms.append((c,s))

    # Weak stiffness for topology
    epsilon = 1e-6

    # ==========================================================
    # INTERPOLATION
    # ==========================================================

    def wc(x,p): return x**p
    def dwc(x,p): return p*x**(p-1)

    def ws(x,q): return x**q
    def dws(x,q): return q*x**(q-1)

    def wd(x,s): return x/(1+s*(1-x))
    def dwd(x,s): return (1+s)/(1+s*(1-x))**2

    # ==========================================================
    # FEM SOLVER
    # ==========================================================

    def solve_structure(x,p):

        x = x.reshape(n_members,n_cand)
        K = np.zeros((total_dof,total_dof))

        for j,(n1,n2) in enumerate(members):
            dofs=[3*n1,3*n1+1,3*n1+2,3*n2,3*n2+1,3*n2+2]
            member_sum = np.sum(x[j])
            for i in range(n_cand):
                K[np.ix_(dofs,dofs)] += wc(x[j,i],p)*K_elem[j][i]
            K[np.ix_(dofs,dofs)] += epsilon*(1-member_sum)*np.eye(6)

        Kff = K[np.ix_(free_dofs,free_dofs)]
        uf = np.linalg.solve(Kff,F[free_dofs])
        u = np.zeros(total_dof)
        u[free_dofs]=uf
        return u,K

    # ==========================================================
    # OBJECTIVE
    # ==========================================================

    def objective(x,p,s):
        u,K = solve_structure(x,p)
        compliance = F@u
        xmat = x.reshape(n_members,n_cand)
        mass=0
        for j in range(n_members):
            for i in range(n_cand):
                mass += wd(xmat[j,i],s)*rho*candidates[i]["A"]*lengths[j]
        return alpha*compliance+(1-alpha)*mass

    # ==========================================================
    # STRESS CONSTRAINTS
    # ==========================================================

    def stress_constraints(x,p,q):

        xmat=x.reshape(n_members,n_cand)
        u,_=solve_structure(x,p)

        constraints=[]

        for j,(n1,n2) in enumerate(members):
            dofs=[3*n1,3*n1+1,3*n1+2,3*n2,3*n2+1,3*n2+2]
            ue=u[dofs]
            L=lengths[j]
            c,s=transforms[j]

            axial_strain=(ue[3]-ue[0])/L
            curvature=(ue[5]-ue[2])/L

            for i,cand in enumerate(candidates):

                sigma_ax=E*axial_strain
                sigma_b=E*cand["c"]*curvature
                sigma_vm=np.sqrt((sigma_ax+sigma_b)**2)

                g = sigma_vm/sigma_y - ws(xmat[j,i],q)
                constraints.append(g)

        return np.array(constraints)

    # ==========================================================
    # GRADIENT
    # ==========================================================

    def gradient(x,p,s):
        xmat=x.reshape(n_members,n_cand)
        u,K=solve_structure(x,p)
        grad=np.zeros_like(xmat)

        for j,(n1,n2) in enumerate(members):
            dofs=[3*n1,3*n1+1,3*n1+2,3*n2,3*n2+1,3*n2+2]
            ue=u[dofs]
            for i in range(n_cand):
                dK = dwc(xmat[j,i],p)*K_elem[j][i]
                dC = -ue@dK@ue
                dM = dwd(xmat[j,i],s)*rho*candidates[i]["A"]*lengths[j]
                grad[j,i]=alpha*dC+(1-alpha)*dM

        return grad.flatten()

    # ==========================================================
    # CONTINUATION
    # ==========================================================

    x0=np.ones((n_members,n_cand))/n_cand
    x_flat=x0.flatten()
    bounds=[(0,1)]*len(x_flat)

    for cont in range(5):

        p=1+cont
        q=max(1,p-1)
        s=cont

        cons=[]

        # topology inequality
        for j in range(n_members):
            idx=slice(j*n_cand,(j+1)*n_cand)
            cons.append({
                "type":"ineq",
                "fun":lambda x,idx=idx: 1-np.sum(x[idx])
            })

        # stress constraints
        cons.append({
            "type":"ineq",
            "fun":lambda x,p=p,q=q: -stress_constraints(x,p,q)
        })

        res=minimize(objective,x_flat,args=(p,s),
                     jac=gradient,method="SLSQP",
                     bounds=bounds,constraints=cons,
                     options={"ftol":1e-6,"maxiter":200})

        x_flat=res.x

        xmat=x_flat.reshape(n_members,n_cand)
        if np.all(np.max(xmat,axis=1)>0.95):
            break

    return {
        "weights":xmat,
        "selected":np.argmax(xmat,axis=1),
        "active_members":np.sum(xmat,axis=1)>0.01,
        "objective":objective(x_flat,p,s)
    }

result = run_full_dmo_frame_optimizer()
print(result["weights"])
print("Selected:", result["selected"])

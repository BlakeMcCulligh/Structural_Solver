import numpy as np
from scipy.optimize import minimize
from scipy.linalg import eigvals


def run_full_dmo_frame_optimizer_with_stability():

    # ==========================================================
    # MATERIAL
    # ==========================================================

    E = 210e9
    rho = 7850
    sigma_y = 235e6
    alpha = 0.3
    lambda_min_required = 1.0

    # ==========================================================
    # GEOMETRY
    # ==========================================================

    nodes = np.array([
        [0.0, 0.0],
        [3.0, 0.0],
        [3.0, 3.0]
    ])

    members = [(0,1),(1,2),(0,2)]

    candidates = [
        {"A":5e-4,"I":8e-6,"c":0.05},
        {"A":8e-4,"I":2e-5,"c":0.06},
        {"A":1.2e-3,"I":4e-5,"c":0.08}
    ]

    n_nodes=len(nodes)
    n_members=len(members)
    n_cand=len(candidates)
    total_dof=3*n_nodes

    fixed=[0,1,2]
    free=list(set(range(total_dof))-set(fixed))

    F=np.zeros(total_dof)
    F[7]=-10000

    epsilon=1e-6

    # ==========================================================
    # BEAM STIFFNESS
    # ==========================================================

    def beam_stiff(E,A,I,x1,y1,x2,y2):

        L=np.hypot(x2-x1,y2-y1)
        c=(x2-x1)/L
        s=(y2-y1)/L

        k_local=np.array([
            [A*E/L,0,0,-A*E/L,0,0],
            [0,12*E*I/L**3,6*E*I/L**2,0,-12*E*I/L**3,6*E*I/L**2],
            [0,6*E*I/L**2,4*E*I/L,0,-6*E*I/L**2,2*E*I/L],
            [-A*E/L,0,0,A*E/L,0,0],
            [0,-12*E*I/L**3,-6*E*I/L**2,0,12*E*I/L**3,-6*E*I/L**2],
            [0,6*E*I/L**2,2*E*I/L,0,-6*E*I/L**2,4*E*I/L]
        ])

        T=np.array([
            [c,s,0,0,0,0],
            [-s,c,0,0,0,0],
            [0,0,1,0,0,0],
            [0,0,0,c,s,0],
            [0,0,0,-s,c,0],
            [0,0,0,0,0,1]
        ])

        return T.T@k_local@T,L,c,s

    # ==========================================================
    # PRECOMPUTE
    # ==========================================================

    K_elem=[]
    lengths=[]
    transforms=[]

    for (n1,n2) in members:
        x1,y1=nodes[n1]
        x2,y2=nodes[n2]
        k_list=[]
        for cand in candidates:
            k,L,c,s=beam_stiff(E,cand["A"],cand["I"],x1,y1,x2,y2)
            k_list.append(k)
        K_elem.append(k_list)
        lengths.append(L)
        transforms.append((c,s))

    # ==========================================================
    # INTERPOLATION
    # ==========================================================

    def wc(x,p): return x**p
    def dwc(x,p): return p*x**(p-1)

    def ws(x,q): return x**q
    def wd(x,s): return x/(1+s*(1-x))

    # ==========================================================
    # SOLVER
    # ==========================================================

    def solve(x,p):

        x=x.reshape(n_members,n_cand)
        K=np.zeros((total_dof,total_dof))

        for j,(n1,n2) in enumerate(members):
            dofs=[3*n1,3*n1+1,3*n1+2,3*n2,3*n2+1,3*n2+2]
            member_sum=np.sum(x[j])
            for i in range(n_cand):
                K[np.ix_(dofs,dofs)] += wc(x[j,i],p)*K_elem[j][i]
            K[np.ix_(dofs,dofs)] += epsilon*(1-member_sum)*np.eye(6)

        Kff=K[np.ix_(free,free)]
        uf=np.linalg.solve(Kff,F[free])
        u=np.zeros(total_dof)
        u[free]=uf
        return u,K

    # ==========================================================
    # GEOMETRIC STIFFNESS
    # ==========================================================

    def geometric_stiffness(x,u,p):

        x=x.reshape(n_members,n_cand)
        Kg=np.zeros((total_dof,total_dof))

        for j,(n1,n2) in enumerate(members):

            dofs=[3*n1,3*n1+1,3*n1+2,3*n2,3*n2+1,3*n2+2]
            ue=u[dofs]
            L=lengths[j]

            axial=(ue[3]-ue[0])/L
            N=E*axial*np.sum(wc(x[j],p))

            k_geo=(N/(30*L))*np.array([
                [0,0,0,0,0,0],
                [0,36,3*L,0,-36,3*L],
                [0,3*L,4*L**2,0,-3*L,-L**2],
                [0,0,0,0,0,0],
                [0,-36,-3*L,0,36,-3*L],
                [0,3*L,-L**2,0,-3*L,4*L**2]
            ])

            Kg[np.ix_(dofs,dofs)] += k_geo

        return Kg

    # ==========================================================
    # BUCKLING CONSTRAINT
    # ==========================================================

    def stability_constraint(x,p):

        u,K=solve(x,p)
        Kg=geometric_stiffness(x,u,p)

        Kff=K[np.ix_(free,free)]
        Kgff=Kg[np.ix_(free,free)]

        eigs=eigvals(Kff,Kgff)
        eigs=np.real(eigs[np.isreal(eigs)])
        lam=np.min(eigs)

        return lam - lambda_min_required

    # ==========================================================
    # OBJECTIVE
    # ==========================================================

    def objective(x,p,s):
        u,_=solve(x,p)
        compliance=F@u
        xmat=x.reshape(n_members,n_cand)
        mass=0
        for j in range(n_members):
            for i in range(n_cand):
                mass+=wd(xmat[j,i],s)*rho*candidates[i]["A"]*lengths[j]
        return alpha*compliance+(1-alpha)*mass

    # ==========================================================
    # CONTINUATION LOOP
    # ==========================================================

    x0=np.ones((n_members,n_cand))/n_cand
    x_flat=x0.flatten()
    bounds=[(0,1)]*len(x_flat)

    for cont in range(5):

        p=1+cont
        q=max(1,p-1)
        s=cont

        cons=[]

        for j in range(n_members):
            idx=slice(j*n_cand,(j+1)*n_cand)
            cons.append({
                "type":"ineq",
                "fun":lambda x,idx=idx: 1-np.sum(x[idx])
            })

        cons.append({
            "type":"ineq",
            "fun":lambda x,p=p: stability_constraint(x,p)
        })

        res=minimize(objective,x_flat,args=(p,s),
                     method="SLSQP",
                     bounds=bounds,
                     constraints=cons,
                     options={"maxiter":200})

        x_flat=res.x

    return {
        "weights":x_flat.reshape(n_members,n_cand),
        "selected":np.argmax(x_flat.reshape(n_members,n_cand),axis=1)
    }

result = run_full_dmo_frame_optimizer_with_stability()
print(result["weights"])
print("Selected:", result["selected"])
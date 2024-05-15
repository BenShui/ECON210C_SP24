import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

# Parameters
gam = 1
phi = 1
chi = 1
beta = 0.99
rho = 0.99

T = 500

# Calibrating theta
nu_arr = [0.25, 0.5, 0.999, 2, 4]
theta_arr = []
for nu in nu_arr:
    theta_arr.append((1 - beta) / ( (1 - beta) + ( (-(1 - beta) + np.sqrt(beta**2 - 2*beta + 5))/2 )**nu))

# Useful Matrices
I = sp.sparse.eye(T)
Ip1 = sp.sparse.diags([np.ones(T-1)], [1], (T, T))
Im1 = sp.sparse.diags([np.ones(T-1)], [-1], (T, T))
O = sp.sparse.csr_matrix((T, T))

# m
m = np.zeros((T, 1))
m[0] = 1
for t in range(1, T):
    m[t] = rho * m[t-1]

# Iteration
## Record c, p, q IRFs
C, P, Q = [], [], []
for i in range(len(nu_arr)):
    nu = nu_arr[i]
    theta = theta_arr[i]

    # Market Clearing Condition
    Phigmc = I
    Phigmy = -I
    Phigmq = O
    Phigmx = O
    Phigmw = O

    # Euler Equation
    Phieulc = nu*I-nu*Ip1
    Phieulx = (gam-nu)*I-(gam-nu)*Ip1
    Phieulq = I
    Phieuly = O
    Phieulw = O

    dHdY = sp.sparse.bmat([[Phigmc, Phigmq, Phigmx, Phigmy, Phigmw],
                       [Phieulc, Phieulq, Phieulx, Phieuly, Phieulw]])

    # dHdU
    Phigmn = O
    Phieuln = O
    Phigmp = O
    Phieulp = I-Ip1
    dHdU = sp.sparse.bmat([[Phigmn, Phigmp],
                       [Phieuln, Phieulp]])

    # dYdU and dYdZ
    Phiyn = I
    Phiyp = O
    Phiyep = O
    Phiwp = I
    Phiwn = O
    Phiwep = O
    dYFdU = sp.sparse.bmat([[Phiyn, Phiyp],
                        [Phiwn, Phiwp]])
    dYFdZ = sp.sparse.bmat([[Phiyep],
                        [Phiwep]])
    omega = (theta**(1/nu)*(1-beta)**(1-1/nu))/((1-theta)**(1/nu)+theta**(1/nu)*(1-beta)**(1-1/nu))
    eta = beta/(nu*(1-beta))
    Phicn = -eta*phi/(gam*eta-(gam-nu)*omega*eta)*I
    Phicp = (gam-nu)*omega*eta/(gam*eta-(gam-nu)*omega*eta)*I
    Phiqn = -phi/(gam*eta-(gam-nu)*omega*eta)*I
    Phiqp = gam/(gam*eta-(gam-nu)*omega*eta)*I
    Phicm = -(gam-nu)*omega*eta/(gam*eta-(gam-nu)*omega*eta)*I
    Phiqm = -gam/(gam*eta-(gam-nu)*omega*eta)*I
    Phixn = (1-omega)*Phicn
    Phixp = (1-omega)*Phicp - omega*I
    Phixm = (1-omega)*Phicm + omega*I
    dYHdU = sp.sparse.bmat([[Phicn, Phicp],
                        [Phiqn, Phiqp],
                        [Phixn, Phixp]])
    dYHdZ = sp.sparse.bmat([[Phicm],
                        [Phiqm],
                        [Phixm]])
    dYdU = sp.sparse.bmat([[dYHdU],
                           [dYFdU]])
    dYdZ = sp.sparse.bmat([[dYHdZ],
                       [dYFdZ]])

    # dHdU and dHdZ
    # Note that here dHdU has additional direct effect term
    dHdU = dHdY @ dYdU + dHdU
    dHdZ = dHdY @ dYdZ

    # Jacobian
    dUdZ = - sp.sparse.linalg.spsolve(dHdU, dHdZ)
    dYdZ = dYdU @ dUdZ + dYdZ
    dXdZ = sp.sparse.bmat([[dUdZ],
                           [dYdZ]])
    X = dXdZ @ m
    # unpack X into its components n,q,c,q,x,y,w
    n = X[0:T]
    p = X[T:2*T]
    c = X[2*T:3*T]
    q = X[3*T:4*T]
    x = X[4*T:5*T]
    y = X[5*T:6*T]
    w = X[6*T:7*T]

    # Record
    C.append(c)
    P.append(p)
    Q.append(q)

plt.figure(figsize=(10,7))
for i in range(len(C)):
    plt.plot(C[i], label='nu='+str(nu_arr[i]))
plt.xlabel('Time')
plt.ylabel('Consumption')
plt.title('Consumption Change')
plt.legend()
plt.show()


plt.figure(figsize=(10,7))
for i in range(len(P)):
    plt.plot(P[i], label='nu='+str(nu_arr[i]))
plt.xlabel('Time')
plt.ylabel('Price')
plt.title('Price Change')
plt.legend()
plt.show()

plt.figure(figsize=(10,7))
for i in range(len(Q)):
    plt.plot(Q[i], label='nu='+str(nu_arr[i]))
plt.xlabel('Time')
plt.ylabel('Nominal Interest Rate')
plt.title('Nominal Interest Rate Change')
plt.legend()
plt.show()
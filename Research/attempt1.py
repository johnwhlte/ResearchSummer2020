# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 09:29:11 2020

@author: johnw
"""

import numpy as np
import matplotlib.pyplot as plt

u_ns = []
xs = []
ts = []

def y(x):
    return 3*x**2 + 2*x

def y2(x,t):
    return 2*x**2*t

def solver_FE_simple(I, a, f, L, dt, F, T):
    """
    Simplest expression of the computational algorithm
    using the Forward Euler method and explicit Python loops.
    For this method F <= 0.5 for stability.
    """
    import time;  t0 = time.process_time()  # For measuring the CPU time

    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)       # Mesh points in space
    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]

    u   = np.zeros(Nx+1)
    u_n = np.zeros(Nx+1)

    # Set initial condition u(x,0) = I(x)
    for i in range(0, Nx+1):
        u_n[i] = I(x[i])

    for n in range(0, Nt):
        # Compute u at inner mesh points
        for i in range(1, Nx):
            u[i] = u_n[i] + F*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + \
                   dt*f(x[i], t[n])

        # Insert boundary conditions
        u[0] = 0;  u[Nx] = 0

        # Switch variables before next step
        #u_n[:] = u  # safe, but slow
        u_n, u = u, u_n

    t1 = time.process_time()
    return u_n, x, t, t1-t0  # u_n holds latest u

h = solver_FE_simple(y,1,y2,10,0.1,0.5,10)

for i in range(0,22):
    u_ns.append(h[0][i])
for i in range(0,22):
    xs.append(h[1][i])
for i in range(0,99):
    ts.append(h[2][i])
    
plt.plot(xs,u_ns,'r-')
plt.show()



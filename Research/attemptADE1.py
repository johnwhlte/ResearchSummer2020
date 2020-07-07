# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 12:45:30 2020

@author: johnw
"""

import numpy as np
import matplotlib.pyplot as plt
import time

#I is the initial condition
#a is the diffusion coefficient
#f is the source function
#L is the length
#dt is the time step
#D is the coefficient on the second derivative
#U is the advecitve coefficient
#T is the length of time (amount of time steps)

u_ns = []
xs = []
ts = []

def initial_value_function(x):
    return 2

def source_or_sink(x,t):
    return 0.001*x**2*t

def solver_ADE_Simple(dx,dt,L,T,a,D,U,I,f):

    t0 = time.process_time()

    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    dx = np.sqrt(a*dt/D)
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
            u[i] = u_n[i] - (U*dt/(2*dx))*(u_n[i+1] - u_n[i-1]) + ((D*dt)/dx**2)*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + \
                   dt*f(x[i], t[n])

        # Insert boundary conditions
        u[0] = 1;  u[Nx] = 0

        # Switch variables before next step
        #u_n[:] = u  # safe, but slow
        u_n, u = u, u_n

    t1 = time.process_time()
    return u_n, x, t, t1-t0  # u_n holds latest u

sol = solver_ADE_Simple(0.1,0.1,10,10,1,0.01,0.005,initial_value_function,source_or_sink)
for i in range(0,len(sol[0])):
    u_ns.append(sol[0][i])
for i in range(0,len(sol[1])):
    xs.append(sol[1][i])
for i in range(0,len(sol[2])):
    ts.append(sol[2][i])

plt.plot(xs,u_ns,'r-')

plt.show()

print(sol)

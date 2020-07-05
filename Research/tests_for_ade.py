import numpy as np
import matplotlib.pyplot as plt
import math

def Diffusion(x_step,y_step,x_length,y_length,u_n,D,time):
    diffusion = [[0 for i in range(0,math.floor(x_length/x_step))] for j in range(0,math.floor(y_length/y_step))]
    for i in range(1,math.floor(x_length/x_step)-1):
        for j in range(1, math.floor(y_length/y_step)-1):
            diffusion[i][j] = time*D*(((u_n[i+1][j] - 2*u_n[i][j] + u_n[i-1][j])/(x_step**2)) + ((u_n[i][j+1] - 2*u_n[i][j] + u_n[i][j-1])/(y_step**2)))
    return diffusion

u = [[0 for i in range(0,10)] for j in range(0,10)]

for i in range(0,10):
    u[0][i] = 1
for i in range(0,10):
    u[-1][i] = 0
for i in range(0,10):
    u[i][0] = 1
for i in range(0,10):
    u[i][-1] = 0

one = u[3][2] + 2

print(one)
for t in range(0,100):
    test1 = Diffusion(0.01,0.01,0.1,0.1,u,1,0.000045)

    for i in range(0,10):
        for j in range(0,10):
            u[i][j] = u[i][j] + test1[i][j]

for i in range(0,10):
    for j in range(0,10):
        if u[i][j] >=0.9:
            plt.plot(i,j,'ro')
        elif u[i][j] >=0.7:
            plt.plot(i,j,'o',color="orange")
        elif u[i][j] >=0.5:
            plt.plot(i,j,'yo')
        elif u[i][j] >= 0.3:
            plt.plot(i,j,'go')
        elif u[i][j] >= 0.1:
            plt.plot(i,j,'bo')
        else:
            plt.plot(i,j,'mo')

plt.show()


print(u)

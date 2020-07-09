import math
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

xlength = 1
ylength = 1
deltax = 0.01
deltay = 0.01
D = 1
deltat = (1/(2*D))*(((deltax*deltay)**2)/((deltax**2) + (deltay**2)))
time = 10*deltat

class Mesh_Object:
    def __init__(self,x,y,u_xy):
        self.x = x
        self.y = y
        self.u_xy = u_xy

x = np.arange(0,xlength,deltax)
y = np.arange(0,xlength,deltay)
xx, yy = np.meshgrid(x,y)

def Diffusion(x_step,y_step,x_length,y_length,u_n,D,x,y):
    diffusion = 0*x + 0*y
    for i in range(1,math.floor(y_length/y_step)-1):
        for j in range(1, math.floor(x_length/x_step)-1):
            diffusion[i][j] = D*(((u_n[i+1][j] - 2*u_n[i][j] + u_n[i-1][j])/(x_step**2)) + ((u_n[i][j+1] - 2*u_n[i][j] + u_n[i][j-1])/(y_step**2)))
    return diffusion

def Advection(xstars,ystars,u_xy,deltax,deltay,xlength,ylength,x,y):
    f = interpolate.interp2d(x,y,u_xy)
    advected = f(xstars,ystars)

    return advected

def Velocity_Field(x_step,y_step,x_length,y_length,x,y):
    velocities = [[[0,0] for i in range(0,math.floor(y_length/y_step))] for j in range(0,math.floor(x_length/x_step))]
    for i in range(0,math.floor(x_length/x_step)):
        for j in range(0,math.floor(y_length/y_step)):
            velocities[i][j] = [1,1]
    return velocities

def find_xstar_ystar(v,x,y,deltax,deltay):
    empty_2d = 0*x + 0*y
    xstars = []
    ystars = []
    for i in range(0,len(empty_2d[:][0])):
        for j in range(0,len(empty_2d[0][:])):
            xstars.append((i + np.sqrt(v[i][j][0]**2 + v[i][j][1]**2))*math.cos(v[i][j][0]/(np.sqrt(v[i][j][0]**2 + v[i][j][1]**2))))
            ystars.append((j + np.sqrt(v[i][j][0]**2 + v[i][j][1]**2))*math.sin(v[i][j][1]/(np.sqrt(v[i][j][0]**2 + v[i][j][1]**2))))

    return xstars,ystars

def initial_condition(x,y):
    return 0*x + 0*y

z = initial_condition(xx,yy)


Two_D_Mesh = Mesh_Object(x,y,z)


for i in range(0,math.floor(ylength/deltay)):
    Two_D_Mesh.u_xy[0][i] = 1

for i in range(0,math.floor(ylength/deltay)):
    Two_D_Mesh.u_xy[-1][i] = 0
for i in range(0,math.floor(xlength/deltax)):
    Two_D_Mesh.u_xy[i][0] = 0

for i in range(0,math.floor(xlength/deltax)):
    Two_D_Mesh.u_xy[i][-1] = 0

diffused = Diffusion(deltax,deltay,xlength,ylength,z,D,xx,yy)
velocity = Velocity_Field(deltax,deltay,xlength,ylength,xx,yy)
xstars,ystars = find_xstar_ystar(velocity,xx,yy,deltax,deltay)
#print(xstars,ystars)
for t in range(0,math.floor(time/deltat)):
    diffused = Diffusion(deltax,deltay,xlength,ylength,z,D,xx,yy)
    advected = Advection(xstars,ystars,Two_D_Mesh.u_xy,deltax,deltay,xlength,ylength,x,y)
    for i in range(0,math.floor(xlength/deltax)):
        for j in range(0,math.floor(ylength/deltay)):
            Two_D_Mesh.u_xy[i][j] = Two_D_Mesh.u_xy[i][j] + deltat*diffused[i][j] + advected[i][j]
            if Two_D_Mesh.u_xy[i][j] > 1:
                Two_D_Mesh.u_xy[i][j] = 1
            elif Two_D_Mesh.u_xy[i][j] <0:
                Two_D_Mesh.u_xy[i][j] = 0

plt.imshow(Two_D_Mesh.u_xy, cmap='viridis', interpolation='nearest')
plt.title('2-D Diffusion with ' + str(math.floor(time/deltat)) + ' time steps')
plt.show()

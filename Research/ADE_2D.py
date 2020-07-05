import math
import matplotlib.pyplot as plt
import numpy as np

#--------setting parameters------------------------
x_step_size = 0.01
y_step_size = 0.01
magnitude_spatial_x = 1
magnitude_spatial_y = 1
Diffusion_Coefficient = 1
time_step = (1/(2*Diffusion_Coefficient))*(((x_step_size*y_step_size)**2)/((x_step_size**2) + (y_step_size**2)))
length_of_time = 10000*time_step
print(time_step)
#-----------defining the mesh class------------
class Mesh_Object:

    def __init__(self,x,y,u_xy):
        self.x = x
        self.y = y
        self.u_xy = u_xy

#------math functions------------------------

def Diffusion(x_step,y_step,x_length,y_length,u_n,D,time):
    diffusion = [[0 for i in range(0,math.floor(x_length/x_step))] for j in range(0,math.floor(y_length/y_step))]
    for i in range(1,math.floor(x_length/x_step)-1):
        for j in range(1, math.floor(y_length/y_step)-1):
            diffusion[i][j] = time*D*(((u_n[i+1][j] - 2*u_n[i][j] + u_n[i-1][j])/(x_step**2)) + ((u_n[i][j+1] - 2*u_n[i][j] + u_n[i][j-1])/(y_step**2)))
    return diffusion
#----initialize the concentration field, define mesh domain ------------
null_u = [[0 for i in range(0,math.floor(magnitude_spatial_x/x_step_size))] for j in range(0,math.floor(magnitude_spatial_y/y_step_size))]

def initial_condition(u,x_step,y_step,x_length,y_length):
    for i in range(0,math.floor(x_length/x_step)):
        for j in range(0,math.floor(y_length/y_step)):
            u[i][j] = 0
    return u

def create_x_and_y_points(x_step,y_step,x_length,y_length):
    x = np.linspace(0,math.floor(x_length/x_step))
    y = np.linspace(0,math.floor(y_length/y_step))

    return x,y

init_u_xy = initial_condition(null_u,x_step_size,y_step_size,magnitude_spatial_x,magnitude_spatial_y)
x_field, y_field = create_x_and_y_points(x_step_size,y_step_size,magnitude_spatial_x,magnitude_spatial_y)

#----------creating the mesh object---------------
Two_D_Mesh = Mesh_Object(x_field,y_field,init_u_xy)
print(Two_D_Mesh.u_xy)

#-----set the boundary conditions-------------------

for i in range(0,math.floor(magnitude_spatial_y/y_step_size)):
    Two_D_Mesh.u_xy[0][i] = 1

for i in range(0,math.floor(magnitude_spatial_y/y_step_size)):
    Two_D_Mesh.u_xy[-1][i] = 0
for i in range(0,math.floor(magnitude_spatial_x/x_step_size)):
    Two_D_Mesh.u_xy[i][0] = 1

for i in range(0,math.floor(magnitude_spatial_x/x_step_size)):
    Two_D_Mesh.u_xy[i][-1] = 0

#-----attempt at time stepping---------------------

for t in range(0,math.floor(length_of_time/time_step)):
    diffused_profile = Diffusion(x_step_size,y_step_size,magnitude_spatial_x,magnitude_spatial_y,Two_D_Mesh.u_xy,Diffusion_Coefficient,time_step)
    for i in range(0,math.floor(magnitude_spatial_x/x_step_size)):
        for j in range(0,math.floor(magnitude_spatial_y/y_step_size)):
            Two_D_Mesh.u_xy[i][j] = Two_D_Mesh.u_xy[i][j] + diffused_profile[i][j]
            if Two_D_Mesh.u_xy[i][j] > 1:
                Two_D_Mesh.u_xy[i][j] = 1
            elif Two_D_Mesh.u_xy[i][j] <0:
                Two_D_Mesh.u_xy[i][j] = 0

print(Two_D_Mesh.u_xy)
plt.imshow(Two_D_Mesh.u_xy, cmap='viridis', interpolation='nearest')
plt.title('2-D Diffusion with ' + str(math.floor(length_of_time/time_step)) + ' time steps')
plt.show()

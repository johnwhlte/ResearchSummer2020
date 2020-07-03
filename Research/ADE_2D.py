import math
import matplotlib.pyplot as plt
import numpy as np

#--------setting parameters------------------------
x_step_size = 0.01
y_step_size = 0.01
magnitude_spatial_x = 1
magnitude_spatial_y = 1
null_u = []
#-----------defining the mesh class------------
class Mesh_Object:

    def __init__(self,x,y,u_xy):
        self.x = x
        self.y = y
        self.u_xy = u_xy

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

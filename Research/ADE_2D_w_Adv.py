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
length_of_time = 100*time_step
velocity = [0.01,3]
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
            diffusion[i][j] = D*(((u_n[i+1][j] - 2*u_n[i][j] + u_n[i-1][j])/(x_step**2)) + ((u_n[i][j+1] - 2*u_n[i][j] + u_n[i][j-1])/(y_step**2)))
    return diffusion

def MOC(v_mag,v_angle,v_direct,delta_t,delta_x,delta_y,x_length,y_length,u):
    u_pstar = [[0 for i in range(0,math.floor(x_length/delta_x))] for j in range(0,math.floor(y_length/delta_y))]
    #find the point to use, calculates the nearest corner
    x_stepper = 1
    y_stepper = -1
    distance_pstar = v_mag*delta_t
    d_xstar = v_mag*delta_t*math.cos(v_angle)
    d_ystar = v_mag*delta_t*math.sin(v_angle)
    num_of_xsteps = math.floor(d_xstar/delta_x)
    num_of_ysteps = math.floor(d_ystar/delta_y)
    #assign signs to the steps
    if v_direct == 'down' or v_direct == 'down_left' or v_direct == 'down-right':
        num_of_ysteps = - num_of_ysteps
        y_stepper = -1
        d_ystar = -d_ystar
    if v_direct == 'left' or v_direct == 'down_left' or v_direct == 'up-left':
        num_of_xsteps = - num_of_xsteps
        x_stepper = -1
        d_xstar = -d_xstar
    if v_direct == 'left' or 'right':
        y_stepper = 0
    if v_direct =='down' or 'up':
        x_stepper = 0
    for i in range(1,math.floor(x_length/delta_x)-1):
        for j in range(1,math.floor(x_length/delta_x)-1):
            store_info_surrounding = [[0 for k in range(0,4)] for l in range(0,3)]
            store_info_surrounding[0][0] = u[i + num_of_xsteps][j + num_of_ysteps]
            store_info_surrounding[0][1] = u[i + num_of_xsteps + x_stepper][j + num_of_ysteps]
            store_info_surrounding[0][2] = u[i + num_of_xsteps][j + num_of_ysteps + y_stepper]
            store_info_surrounding[0][3] = u[i + num_of_xsteps + x_stepper][j + num_of_ysteps + y_stepper]
            store_info_surrounding[1][0] = abs(((j+d_ystar) - (j + num_of_ysteps))/((i + d_xstar) - (i + num_of_xsteps)))
            store_info_surrounding[1][1] = abs(((j+d_ystar) - (j + num_of_ysteps))/((i + d_xstar) - (i + num_of_xsteps + x_stepper)))
            store_info_surrounding[1][2] = abs(((j+d_ystar) - (j + num_of_ysteps + y_stepper))/((i + d_xstar) - (i + num_of_xsteps)))
            store_info_surrounding[1][3] = abs(((j+d_ystar) - (j + num_of_ysteps + y_stepper))/((i + d_xstar) - (i + num_of_xsteps + x_stepper)))
            store_info_surrounding[2][0] = store_info_surrounding[1][0]/(store_info_surrounding[1][0] + store_info_surrounding[1][1] + store_info_surrounding[1][2] + store_info_surrounding[1][3])
            store_info_surrounding[2][1] = store_info_surrounding[1][1]/(store_info_surrounding[1][0] + store_info_surrounding[1][1] + store_info_surrounding[1][2] + store_info_surrounding[1][3])
            store_info_surrounding[2][2] = store_info_surrounding[1][2]/(store_info_surrounding[1][0] + store_info_surrounding[1][1] + store_info_surrounding[1][2] + store_info_surrounding[1][3])
            store_info_surrounding[2][3] = store_info_surrounding[1][3]/(store_info_surrounding[1][0] + store_info_surrounding[1][1] + store_info_surrounding[1][2] + store_info_surrounding[1][3])

            u_pstar[i][j] = store_info_surrounding[0][0]*store_info_surrounding[2][0] + store_info_surrounding[0][1]*store_info_surrounding[2][1] + store_info_surrounding[0][2]*store_info_surrounding[2][2] + store_info_surrounding[0][3]*store_info_surrounding[2][3]

    return u_pstar
#-----directional,magnitude,angle check for MOC------------------
def check_direction_and_magnitude(v):
    if v[0] == 0 and v[1] == 0:
        return 0
    elif v[0] == 0 and v[1] > 0:
        direction = 'up'
    elif v[0] == 0 and v[1] < 0:
        direction = 'down'
    elif v[0] < 0 and v[1] == 0:
        direction = 'left'
    elif v[0] > 0 and v[1] == 0:
        direction = 'right'
    elif v[0] <0 and v[1] < 0:
        direction = 'down-left'
    elif v[0] <0 and v[1] >0:
        direction = 'up-left'
    elif v[0] >0 and v[1] < 0:
        direction = 'down-left'
    elif v[0] > 0 and v[1] > 0:
        direction = 'up-right'
    else:
        print('error in velocity')

    magnitude = math.sqrt(v[0]**2 + v[1]**2)
    angle = math.acos((abs(v[0])/magnitude))

    return direction,magnitude,angle

vel_dir,vel_mag,vel_ang = check_direction_and_magnitude(velocity)

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
    advected_profile = MOC(vel_mag,vel_ang,vel_dir,time_step,x_step_size,y_step_size,magnitude_spatial_x,magnitude_spatial_y,Two_D_Mesh.u_xy)
    for i in range(0,math.floor(magnitude_spatial_x/x_step_size)):
        for j in range(0,math.floor(magnitude_spatial_y/y_step_size)):
            Two_D_Mesh.u_xy[i][j] = Two_D_Mesh.u_xy[i][j] + time_step*diffused_profile[i][j] + advected_profile[i][j]
            if Two_D_Mesh.u_xy[i][j] > 1:
                Two_D_Mesh.u_xy[i][j] = 1
            elif Two_D_Mesh.u_xy[i][j] <0:
                Two_D_Mesh.u_xy[i][j] = 0

print(Two_D_Mesh.u_xy)
plt.imshow(Two_D_Mesh.u_xy, cmap='viridis', interpolation='nearest')
plt.title('2-D Diffusion with ' + str(math.floor(length_of_time/time_step)) + ' time steps')
plt.show()

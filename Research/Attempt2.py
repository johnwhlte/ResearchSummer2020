#importing the necessary packages
import matplotlib.pyplot as plt
import numpy as np
import math

#tunable parameters for setting the constant values of the ADE
length_of_time = 100
time_step_size = 0.02
x_step_size = 0.5*time_step_size
Diffusivity_Coefficient = 0.0005
Advection_Coefficient = 0.05
Boundary1 = 1
Boundary2 = 0
x_initial = 0
length = 2
sigma = x_step_size

#Here we will create a field object, which will store the basic domain information, as well as the concentration
#profile initially and as we calculate over the time steps
class One_D_Field:

    def __init__(self,x,u_x): #initializer function for python that sets characteristic variables for objects that inherit from this class
        self.x = x #setting the characteristic variable "domain"
        self.u_x = u_x #setting the caracteristic variable "concentration profile"

#Below is the diffusion equation using Euler forward method of FDE. The inputs are as follows...
#x = Domain, u_x = concentration profile, BC1 = boundary condition at first edge, BC2 = boundary condition at the second edge
#dx = size of the mesh (spacing between the mesh points), dt = time step length, D = Diffusivity, SOS = Source or sink if
#applicable, t = time variable

def Diffusion_Function(x,u_x,BC1,BC2,dx,dt,D,SOS,t):
    for i in range(1,len(x)-2):
        u_x[i] = u_x[i] + ((D*dt)/dx**2)*(u_x[i+1] - 2*u_x[i] + u_x[i-1]) + SOS(x[i],t)

    u_x[0] = BC1
    u_x[len(x) - 1] = BC2

    return u_x #returns the diffused profile

#Below is the advection equation using the Euler forward method for FDE. All of the input variables are the same
#as the Diffusion_Function, except the D is replaced with U, the advection coefficient.

def Advection_Function(x,u_x,BC1,BC2,dx,dt,U,SOS,t):
    for i in range(1,len(x)-2):
        u_x[i] = u_x[i] - ((U*dt)/(2*dx))*(u_x[i+1]-u_x[i-1]) + SOS(x[i],t)

    u_x[0] = BC1
    u_x[len(x) - 1] = BC2

    return u_x #returns the advected profile

def Source_or_Sink(x,t): #defines the source or sink function
    return 0

def Initial_Value(x,sigma): #defines the initial value condition
    return np.exp((x**2)/(2*sigma)**2)

domain_x = np.linspace(0,length,math.floor((length/x_step_size)))#defines the domain
conc_profile = np.zeros(len(domain_x)) #intializes a concentration profile

One_Dimensional_Field = One_D_Field(domain_x, conc_profile) #constructs the mesh as our object, inheriting from the One_D_Field class.

for i in range(0,len(One_Dimensional_Field.x)): #loop that sets the initial concentration profile to the initial condition.
    One_Dimensional_Field.u_x[i] = Initial_Value(One_Dimensional_Field.x[i],sigma)
print(One_Dimensional_Field.u_x)
#Below is the loop which we apply multiple time steps. The idea here is that we are able to choose simply diffusion
#through the Diffusion_Function and simply advection through the Advection_Function or both by commenting out the
#one we do not need. Here the values are placed into the field object for holding ad subsequent calculations.
#Then each concentration profile is plotted one by one, so you can see the trend over many iterations
for i in range (0,math.floor(length_of_time/time_step_size)): #the loop counter is determined by the length of time over the time step
    calc_Diffusion = Diffusion_Function(One_Dimensional_Field.x,One_Dimensional_Field.u_x,Boundary1,Boundary2,x_step_size,time_step_size,Diffusivity_Coefficient,Source_or_Sink,(i*time_step_size))
    One_Dimensional_Field.u_x = calc_Diffusion #stores the field object's concentration profile as the diffused profile
    calc_Advection = Advection_Function(One_Dimensional_Field.x,One_Dimensional_Field.u_x,Boundary1,Boundary2,x_step_size,time_step_size,Advection_Coefficient,Source_or_Sink,(i*time_step_size))
    One_Dimensional_Field.u_x = calc_Advection #stores the field object's concentration profile as the advected profile
    #plt.plot(One_Dimensional_Field.x,One_Dimensional_Field.u_x,label=str(i)) #plots the concentration profile over the spatial variable

print("Pe = " + str((Advection_Coefficient*time_step_size/(x_step_size))/((Diffusivity_Coefficient*time_step_size)/(x_step_size**2))))
#Calls for the plot of the legend and to display the plot
plt.plot(One_Dimensional_Field.x,One_Dimensional_Field.u_x,label=str(i)) #plots the concentration profile over the spatial variable
plt.legend()
plt.show()

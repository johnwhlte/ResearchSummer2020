#importing the necessary packages
import matplotlib.pyplot as plt
import numpy as np
import math

#tunable parameters for setting the constant values of the ADE
only_diff1_adv2_both3=3  # only diffusion=1, only advection=2, both advection and diffusion=3
length = 1
length_of_time =0.3
Diffusivity_Coefficient = 1
constant_velocity = 10 # Pe
x_step_size =0.01
Boundary1 = 0
Boundary2 = 0
x_initial = 0
sigma = 4*x_step_size
#------------------------------------------------------------------------------
#Here we will create a field object, which will store the basic domain information, as well as the concentration
#profile initially and as we calculate over the time steps
class One_D_Field:

    def __init__(self,x,u_x,prev_u_x): #initializer function for python that sets characteristic variables for objects that inherit from this class
        self.x = x #setting the characteristic variable "domain"
        self.u_x = u_x #setting the caracteristic variable "concentration profile"
        self.prev_u_x = prev_u_x

#------------------------------------------------------------------------------
#choosing time step based on CFL condtion

def Set_Timestep(dx,D,V,diff1_adv2_both3):
    dt_diff = 0.45*dx**2/D #  D dt/dx^2 < 0.5
    dt_adv  = 0.9*dx/V
    if diff1_adv2_both3==1:
        dt=dt_diff
    elif diff1_adv2_both3==2:
            dt=dt_adv
    else:
            dt =min(dt_diff,dt_adv)
    return dt #returns the diffused profile
#------------------------------------------------------------------------------
#Below is the diffusion equation using Euler forward method of FDE. The inputs are as follows...
#x = Domain, u_x = concentration profile, BC1 = boundary condition at first edge, BC2 = boundary condition at the second edge
#dx = size of the mesh (spacing between the mesh points), dt = time step length, D = Diffusivity, SOS = Source or sink if
#applicable, t = time variable

def Diffusion_Function(x,u_n,dx,dt,D):
    diffusion = np.zeros(len(u_n))
    for i in range(1,len(x)-2):
        diffusion[i] = ((D)/dx**2)*(u_n[i+1] - 2*u_n[i] + u_n[i-1])
    return diffusion #returns the diffused profile

#------------------------------------------------------------------------------
#Below is the advection equation using the Euler forward method for FDE. All of the input variables are the same
#as the Diffusion_Function, except the D is replaced with U, the advection coefficient.

def Advection_Function(x,u_n,dx,dt,U):
    advection = np.zeros(len(u_n))
    for i in range(1,len(x)-2):
        advection[i] =  - ((U)/(2*dx))*(u_n[i+1]-u_n[i-1])
    return advection #returns the advected profile
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def Source_or_Sink(x,t): #defines the source or sink function
    return 0

#------------------------------------------------------------------------------
def Initial_Value(x,sigma): #defines the initial value condition
    return np.exp(-((x-0.5)**2)/(2*sigma**2))
    #return 0
#______________________________________________________________________________
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
domain_x = np.linspace(0,length,math.floor((length/x_step_size)))#defines the domain
time_step_size=Set_Timestep(x_step_size,Diffusivity_Coefficient,constant_velocity,only_diff1_adv2_both3)
conc_profile = np.zeros(len(domain_x)) #intializes a concentration profile

One_Dimensional_Field = One_D_Field(domain_x, conc_profile, conc_profile) #constructs the mesh as our object, inheriting from the One_D_Field class.

for i in range(1,len(One_Dimensional_Field.x)-1): #loop that sets the initial concentration profile to the initial condition.
    One_Dimensional_Field.u_x[i] = Initial_Value(One_Dimensional_Field.x[i],sigma)
One_Dimensional_Field.u_x[0]=Boundary1;
One_Dimensional_Field.u_x[-1]=Boundary2;
print(One_Dimensional_Field.u_x)

#Below is the loop which we apply multiple time steps. The idea here is that we are able to choose simply diffusion
#through the Diffusion_Function and simply advection through the Advection_Function or both by commenting out the
#one we do not need. Here the values are placed into the field object for holding ad subsequent calculations.
#Then each concentration profile is plotted one by one, so you can see the trend over many iterations

for i in range (0,math.floor(length_of_time/time_step_size)): #the loop counter is determined by the length of time over the time step
    One_Dimensional_Field.prev_u_x = One_Dimensional_Field.u_x
    if i%50 == 0:
        plt.plot(One_Dimensional_Field.x,One_Dimensional_Field.u_x,label=str(i/10)) #plots the concentration profile over the spatial variable
    calc_Diffusion = time_step_size*Diffusion_Function(One_Dimensional_Field.x,One_Dimensional_Field.u_x,x_step_size,time_step_size,Diffusivity_Coefficient)
    calc_Advection = time_step_size*Advection_Function(One_Dimensional_Field.x,One_Dimensional_Field.u_x,x_step_size,time_step_size,constant_velocity)
    One_Dimensional_Field.u_x = One_Dimensional_Field.u_x+ np.multiply(int(only_diff1_adv2_both3!=2),calc_Diffusion)+ np.multiply(int(only_diff1_adv2_both3!=1),calc_Advection)+ time_step_size*Source_or_Sink(One_Dimensional_Field.u_x,i*time_step_size) #stores the field object's concentration profile as the advected profile    One_Dimensional_Field.u_x[One_Dimensional_Field.u_x>1]=1
    One_Dimensional_Field.u_x[One_Dimensional_Field.u_x>1]=1
    One_Dimensional_Field.u_x[One_Dimensional_Field.u_x<0]=0
    One_Dimensional_Field.u_x[0]=Boundary1;
    One_Dimensional_Field.u_x[-1]=Boundary2;

print("Pe = " + str(constant_velocity))

#Calls for the plot of the legend and to display the plot
#plt.plot(One_Dimensional_Field.x,One_Dimensional_Field.u_x,label=str(i)) #plots the concentration profile over the spatial variable
plt.xlabel("X")
plt.ylabel("Concentration")
plt.title("Pe = " + str(constant_velocity)+" each line is profile after 100 time steps")
plt.show()
#------------------------------------------------------------------------------

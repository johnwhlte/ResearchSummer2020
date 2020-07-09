import math
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
#--------setting parameters------------------------
x_step_size = 0.01
y_step_size = 0.01
c = 0.1
Lx = 2
Ly = 1
One_D_Velocity = 1
Diffusion_Coefficient = 1
time_step = (1/(2*Diffusion_Coefficient))*(((x_step_size*y_step_size)**2)/((x_step_size**2) + (y_step_size**2)))
#time_step = c*x_step_size/One_D_Velocity
length_of_time = 150*time_step
sigma = 4*x_step_size
alpha = 1

#-----------defining the mesh class------------
class Mesh_Object_2D:
    def __init__(self,Lx,Ly,dx,dy):
        self.Lx=Lx
        self.Ly=Ly
        self.dx=dx
        self.dy=dy
        self.nx=math.floor(Lx/dx)
        self.ny=math.floor(Ly/dy)
        self.linespacex=np.linspace(0,1,self.nx)
        self.linespacey=np.linspace(0,1,self.ny)
        self.x ,self.y=np.meshgrid(self.linespacex,self.linespacey)
        self.u_xy =  [[0 for i in range(0,self.nx)] for j in range(0,self.ny)]

    def set_initial_condition(self):
        self.u_xy=  [[0 for i in range(0,self.nx)] for j in range(0,self.ny)]
    def set_boundary_conditions(self):
        self.set_leftNBC()
        self.set_rightNBC()
        self.set_topDBC()
        self.set_bottomDBC()

    def set_leftDBC(self):
        for i in range(0,self.ny):
            self.u_xy[0][i] = 1
    def set_rightDBC(self):
        for i in range(0,self.ny):
            self.u_xy[-1][i] = 1
    def set_topDBC(self):
        for i in range(0,self.ny):
            self.u_xy[i][0] = 1
    def set_bottomDBC(self):
        for i in range(0,self.ny):
            self.u_xy[i][-1] = 0

    def set_leftNBC(self):
        for i in range(0,self.ny):
            self.u_xy[0][i] = self.u_xy[1][i]
    def set_rightNBC(self):
        for i in range(0,self.ny):
            self.u_xy[-1][i] = self.u_xy[-2][i]
    def set_topNBC(self):
        for i in range(0,self.ny):
            self.u_xy[i][0] = self.u_xy[i][1]
    def set_bottomNBC(self):
        for i in range(0,self.ny):
            self.u_xy[i][-1] = self.u_xy[i][-2]

    def set_diffusion_coefficient(self,D):
        self.D=D

    def set_time(self,T):
        self.T=T
    def set_time_step(self,dt):
        self.dt=dt
        self.nt=math.floor(self.T/self.dt)
    def set_velocity(self,vx,vy):
        self.velocity= [[[vx,vy] for i in range(0,self.nx)] for j in range(0,self.ny)]

    def set_u_xy(self, u):
        self.u_xy=u
        for i in range(0,self.nx):
            for j in range(0,self.ny):
                if self.u_xy[i][j] > 1:
                    self.u_xy[i][j] = 1
                elif self.u_xy[i][j] <0:
                    self.u_xy[i][j] = 0
    def DiffusionStep(self):
        diffusion = [[0 for i in range(0,self.nx)] for j in range(0,self.ny)]
        u_n=self.u_xy
        for i in range(1,self.nx-1):
            for j in range(1, self.ny-1):
                diffusion[i][j] = self.D*(((u_n[i+1][j] - 2*u_n[i][j] + u_n[i-1][j])/(self.dx**2)) + ((u_n[i][j+1] - 2*u_n[i][j] + u_n[i][j-1])/(self.dy**2)))
        return diffusion
    def MOC(self):
        advective = [[0 for i in range(0,self.nx)] for j in range(0,self.ny)]
        #find the point to use, calculates the nearest corner

        f = interp2d(self.linespacex, self.linespacey, np.array(self.u_xy), kind='linear')
        for p in range(1, self.nx-1):
            for q in range(1, self.ny-1):
                x_star=self.x[p][q]-self.velocity[p][q][0]*self.dt
                y_star=self.x[p][q]-self.velocity[p][q][1]*self.dt
                advective[p][q] = np.double(f(x_star, y_star))
        return advective


    def advection_diffusion(self):
        for t in range(0,self.nt):
            advected_profile = self.MOC()
            self.set_u_xy(advected_profile)
            self.set_boundary_conditions()
            diffused_profile = self.DiffusionStep()
            self.set_u_xy(self.u_xy + np.multiply(self.dt,diffused_profile) )
            self.set_boundary_conditions()

class Mesh_Object_1D:

        def __init__(self,Lx,dx):
            self.Lx=Lx
            self.dx=dx
            self.nx=math.floor(Lx/dx)
            self.linspacex=np.linspace(0,Lx,self.nx)
            self.x = self.linspacex
            self.u_x =  [0 for i in range(0,self.nx)]

        def set_initial_condition(self,sigma):
            self.u_x =  [0 for i in range(0,self.nx)]
            for i in range(0,self.nx):
                self.u_x[i] = np.exp(-((self.x[i]-1)**2/(2*sigma**2)))


        def set_boundary_conditions(self):
            self.set_leftDBC()
            self.set_rightNBC()

        def set_leftDBC(self):
            self.u_x[0] = 0

        def set_rightNBC(self):
            self.u_x[-1] = self.u_x[-2]

        def set_diffusion_coefficient(self,D):
            self.D=D

        def set_time(self,T):
            self.T=T

        def set_time_step(self,dt):
            self.dt=dt
            self.nt=math.floor(self.T/self.dt)

        def set_velocity(self,vx):
            self.velocity = np.zeros(self.nx)
            for i in range(0,self.nx):
                self.velocity[i] = vx

        def set_u_x(self, u):
            self.u_x=u
            for i in range(0,self.nx):
                if self.u_x[i] > 1:
                    self.u_x[i] = 1
                elif self.u_x[i] <0:
                    self.u_x[i] = 0

        def DiffusionStep(self):
            diffusion = np.zeros(len(self.u_x))
            for i in range(1,self.nx-1):
                diffusion[i] = ((self.D)/self.dx**2)*(self.u_x[i+1] - 2*self.u_x[i] + self.u_x[i-1])
            return diffusion

        def MOC(self):
            advective = np.zeros(len(self.u_x))
            f = interpolate.interp1d(self.linspacex,np.array(self.u_x),kind="linear")
            for i in range(1,self.nx-1):
                xstar = self.x[i] - self.velocity[i]*self.dt
                advective[i] = np.double(f(xstar))
            return advective

        def advection_diffusion(self):
            for t in range(0,self.nt):
                advected_profile = self.MOC()
                self.set_u_x(advected_profile)
                self.set_boundary_conditions()
                diffused_profile = self.DiffusionStep()
                self.set_u_x(self.u_x + np.multiply(self.dt,diffused_profile) )
                self.set_boundary_conditions()

#------main code----------------------------------
#----------creating the mesh object---------------
Two_D_Mesh = Mesh_Object_1D(Lx,x_step_size)
Two_D_Mesh.set_initial_condition(sigma)
print(Two_D_Mesh.u_x)
print(Two_D_Mesh.nx)
Two_D_Mesh.set_boundary_conditions()
Two_D_Mesh.set_time(length_of_time)
Two_D_Mesh.set_time_step(time_step)
Two_D_Mesh.set_diffusion_coefficient(Diffusion_Coefficient)
Two_D_Mesh.set_velocity(One_D_Velocity)
Two_D_Mesh.advection_diffusion()

#u_analytic = [0 for i in range(0,math.floor(Lx/x_step_size))]
#for i in range(0,math.floor(Lx/x_step_size)):
#    u_analytic[i] = (sigma/np.sqrt(sigma**2 + 2*alpha*length_of_time))*np.exp(-(Two_D_Mesh.x[i] - 1 - One_D_Velocity*length_of_time)/(2*(sigma**2 + 2*alpha*length_of_time)))


#-------plotting--------------------------
print(2*math.sqrt(time_step))
#plt.imshow(Two_D_Mesh.u_xy, cmap='viridis', interpolation='nearest')
#plt.title('2-D Diffusion with ' + str(math.floor(length_of_time/time_step)) + ' time steps')
plt.plot(Two_D_Mesh.x,Two_D_Mesh.u_x)
#plt.plot(Two_D_Mesh.x,u_analytic)
plt.show()

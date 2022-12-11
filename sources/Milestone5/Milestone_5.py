from numpy import linspace, array, cos, sin, transpose, zeros, size, concatenate
from esquemas_temporales import Euler, RK4, Euler_Inv, CN, LF
from Cauchy_problem import Cauchy
from Physics import F_N_body, initial_conditions_Nbody


import matplotlib.pyplot as plt

## Temporal variables 

T = 20        # Integration time 
dt = 0.01     # Integration step 

N = int(T/dt)            # Number of steps
t  = linspace(0,T,N)     # Integration time array

# Integration methods

time_schemes = [RK4]

# Number of bodies and initial conditions

Nb = 4

U0 = zeros(Nb*2*3)      # Number of bodies * position and velocity * dimensions

r, v = initial_conditions_Nbody(Nb)     # Initial conditions of different number of bodies


for i in range(len(r)):     # U0 built with the position and velocities in each direction, aka. [x, x', y, y', z, z']

    U0[2*i] = r[i]
    U0[2*i + 1] = v[i]


for scheme in time_schemes:   # Allows different integration methods

        U = Cauchy(F_N_body,t,U0,scheme)

        # 2D representation in XY plane

        fig, twoD = plt.subplots(figsize=(5,5))
        
        twoD.plot(U[0,:], U[2,:], "b")
        twoD.plot(U[6,:], U[8,:], "r")
        if Nb == 3:
            twoD.plot(U[12,:],U[14,:],"k")
        elif Nb == 4:
            twoD.plot(U[12,:],U[14,:],"k")
            twoD.plot(U[18,:],U[20,:],"purple")
        twoD.set_xlabel("X")
        twoD.set_ylabel("Y",rotation = 0)
        twoD.set_title("Proyección en el plano XY de " + str(Nb) + " cuerpos")

        twoD.grid()

        # 3D representation 

        fig = plt.figure()
        ax1 = fig.add_subplot(111,projection='3d')
        ax1.plot_wireframe(U[0,:].reshape((-1, 1)), U[2,:].reshape((-1, 1)), U[4,:].reshape((-1, 1)), color= "red", label = 'Primer cuerpo')
        ax1.plot_wireframe(U[6,:].reshape((-1, 1)), U[8,:].reshape((-1, 1)), U[10,:].reshape((-1, 1)), color= "blue", label = 'Segundo cuerpo')
        if Nb == 3:
            ax1.plot_wireframe(U[12,:].reshape((-1, 1)), U[14,:].reshape((-1, 1)), U[16,:].reshape((-1, 1)), color= "black", label = 'Tercer cuerpo')
        elif Nb == 4:
            ax1.plot_wireframe(U[12,:].reshape((-1, 1)), U[14,:].reshape((-1, 1)), U[16,:].reshape((-1, 1)), color= "black", label = 'Tercer cuerpo')
            ax1.plot_wireframe(U[18,:].reshape((-1, 1)), U[20,:].reshape((-1, 1)), U[22,:].reshape((-1, 1)), color= "purple", label = 'Cuarto cuerpo')
        
        plt.title(str(Nb) + ' cuerpos con metodo ' + scheme.__name__ + ' y dt = ' + str(dt))
        plt.xlabel("X")
        plt.ylabel("Y",rotation = 0)
        plt.grid()

        plt.legend(loc = 'best')
        plt.show()
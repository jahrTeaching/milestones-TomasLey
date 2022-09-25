from numpy import zeros, linspace, array, sqrt
import matplotlib.pyplot as plt
import metodos as m 

## Initial values

r0 = [1,0]
v0 = [0,1]
U0 = r0 + v0

## Integration steps

N = 2000
dt = 0.01

## Problem definition

U = zeros((len(U0), N+1))
U[:,0] = U0

## Euler method solution and plot

U_Euler = U
U_Euler = m.Euler(U_Euler,N,dt)

fig, ax = plt.subplots(figsize=(5,5))
ax.plot(U_Euler[0,:], U_Euler[1,:])
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Metodo Euler')
ax.legend()
plt.show()


## RK4 method solution and plot

U_RK4 = U
U_RK4 = m.RK4(U_RK4,N,dt)

fig, ay = plt.subplots(figsize=(5,5))
ay.plot(U_RK4[0,:], U_RK4[1,:])
ay.set_xlabel('x')
ay.set_ylabel('y')
ay.set_title('Metodo RK4')
ay.legend()
plt.show()

## Checking solutions for different time-steps

dtvar = [0.1, 0.01, 0.001]

U_Euler_var = zeros((len(dtvar), len(U0), N+1))
U_RK4_var = zeros((len(dtvar), len(U0), N+1))

U_Euler_var[:,:,0] = U0
U_RK4_var[:,:,0] = U0

for i in range(len(dtvar)):

    U_Euler_var[i,:,:] = m.Euler(U_Euler_var[i,:,:],N,dtvar[i])
    U_RK4_var[i,:,:] = m.RK4(U,N,dtvar[i])


fig, axv = plt.subplots(figsize=(5,5))
fig, ayv = plt.subplots(figsize=(5,5))

for i in range(len(dtvar)):

    axv.plot(U_Euler_var[i,0,:], U_Euler_var[i,1,:], label='dt = ' + str(dtvar[i]))
    ayv.plot(U_RK4_var[i,0,:], U_RK4_var[i,1,:], label='dt = ' + str(dtvar[i]))

axv.set_title('Metodo Euler con dt variable')
ayv.set_title('Metodo RK4 con dt variable')
axv.legend()
ayv.legend()

plt.show()

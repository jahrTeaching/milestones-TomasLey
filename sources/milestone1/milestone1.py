from numpy import zeros, linspace, array, sqrt
import matplotlib.pyplot as plt
import metodos as m 

## Initial values

r0 = [1,0]
v0 = [0,1]
U0 = r0 + v0

## Integration steps

N = 1000
dt = 0.01

## Problem definition

U = zeros((len(U0), N+1))
U[:,0] = U0

## Euler method solution and plot

U_Euler = m.Euler(U,N,dt)

fig, ax = plt.subplots(figsize=(4, 4))
ax.plot(U_Euler[0,:], U_Euler[1,:], label='Euler Method')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Euler')
ax.legend()
plt.show()

## RK4 method solution and plot

U_RK4 = m.RK4(U,N,dt)

fig, ay = plt.subplots(figsize=(4, 4))
ay.plot(U_RK4[0,:], U_RK4[1,:], label='RK4 Method')
ay.set_xlabel('x')
ay.set_ylabel('y')
ay.set_title('RK4')
ay.legend()
plt.show()

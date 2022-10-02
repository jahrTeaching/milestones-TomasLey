from numpy import zeros, linspace, array, sqrt
from esquemas_temporales import Euler, RK4, Euler_Inv, CN
from Cauchy_problem import Cauchy
from Orbit import Kepler
import matplotlib.pyplot as plt

## Initial values for kepler orbit

U0 = [1, 0, 0, 1]

## Integration steps

N = 500
t_fin = 10

t = linspace(0, t_fin, N+1)

## Problem definition

temp_sch = [Euler, Euler_Inv, RK4, CN]
title_sch = ["Euler", "Inverse Euler", "Runge-Kutta 4", "Crank-Nicholson"]

i = 0

for scheme in temp_sch:

    U = Cauchy(Kepler, t, U0, scheme)

    x = U[:,0]
    y = U[:,1]

    Fx = U[:,2]
    Fy = U[:,3]

    plt.figure(figsize=(5,5))
    plt.plot(x,y)
    plt.title(title_sch[i])
    plt.show()

    i = i+1





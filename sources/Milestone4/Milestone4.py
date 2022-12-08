from numpy import linspace, array, cos, sin, transpose
from esquemas_temporales import Euler, RK4, Euler_Inv, CN, LF
from Cauchy_problem import Cauchy
from Physics import Oscilator, Kepler
from estabilidad import region_estabilidad

import matplotlib.pyplot as plt

## Initial values for kepler orbit

r0 = [0]
v0 = [1]

U0 = array(r0 + v0)

## Integration steps

t_fin = 30

N = 1000

t = linspace(0,t_fin,N)

d_t = array([1, 0.1, 0.01])

time_schemes = [Euler, RK4, CN, Euler_Inv]
title_sch = ["Euler", "Runge-Kutta 4", "Crank-Nicholson", "Inverse Euler"]

for scheme in time_schemes:

    U = Cauchy( Oscilator, t, U0, scheme)    

    fig, sol = plt.subplots(figsize=(5,5))

    sol.plot(t, U[0,:],color = "m", label = "Approximate solution")
    sol.plot(t, sin(t),"--",color = "r", label = "Analitic solution")
    sol.plot(t, U[1,:],color = "g", label = "Appr. velocity")
    sol.plot(t, cos(t),"--", color = "b", label = "Analit. velocity")
    sol.set_title("Harmonic Oscillator dt = 0.01, for "+ title_sch[time_schemes.index(scheme)]+ " scheme")
    sol.set_xlabel("t")
    sol.set_ylabel("y(t)/dy(t)")
    sol.legend(loc = "lower left")
    sol.grid()
    plt.show()

    #stab_region = region_estabilidad_trampa(x,y,scheme)
    rho, x, y  = region_estabilidad(scheme, N, -4, 4, -4, 4)

    fig, rs = plt.subplots(figsize=(5,5))  
    
    for dt in d_t:
        rs.plot([0,0],[dt,-dt], 'o', label = "Roots for dt = " +str(dt))

    rs.contour(x, y, transpose(rho), levels = [0, 1], colors = "r") 
    rs.contourf(x, y, transpose(rho), levels = [0, 1], colors = '#626262') 
    rs.set_title("Stability Region for " + title_sch[time_schemes.index(scheme)]+ " scheme")
    rs.set_xlabel("Re(r)")
    rs.set_ylabel("Im(r)")
    rs.legend(loc = "lower right")
    rs.grid()
    plt.show()







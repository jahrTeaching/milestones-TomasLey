
from numpy import zeros, linspace, array, sqrt
from esquemas_temporales import Euler, RK4, Euler_Inv, CN
from Cauchy_problem import Cauchy
from Orbit import Kepler
from error_schemes import Error_y_convergencia, Convergence, Richardson
import matplotlib.pyplot as plt

## Initial values for kepler orbit

r0 = [1, 0]
v0 = [0, 1]

U0 = array(r0 + v0)

## Integration steps


t_fin = 20

N1 = 2000                            # Condiciones iniciales para Euler Inverso, ya que 
t1 = linspace(0, t_fin, N1+1)        # da problemas de convergencia para dt = 0.1

N2 = 200                             # Condiciones iniciales para el resto de esquemas
t2 = linspace(0, t_fin, N2+1)

## Problem definition

Error_y_convergencia(Cauchy, Kepler, t1, U0, Euler_Inv, q_order = 1,mpoints = 8)


mpoints = [14, 10, 8]           # Número de puntos evaluados en el gráfico log(Error) - log(N)
temp_sch = [Euler, RK4, CN]     # Esquemas temporales
q_order = [1.,4.,2.]            # Orden de convergencia de los esquemas

for scheme in temp_sch:

    Error_y_convergencia(Cauchy, Kepler, t2, U0, scheme, q_order[temp_sch.index(scheme)],mpoints[temp_sch.index(scheme)])






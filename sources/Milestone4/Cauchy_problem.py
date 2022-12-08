from numpy import zeros, linspace
from esquemas_temporales import Euler, RK4, Euler_Inv, CN, LF

"""
    Inputs:
        Phy: Problema físico a resolver
        t: vector de todos los valores temporales que se van a estudiar donde t[i] = t0 + i*dt 
        U0: Condiciones iniciales
        Esq_temp: Esquema temporal utilizado para integrar el problema de Cauchy

    Outputs:
        U: Vector de estado en tn+1
"""

def Cauchy(Phy , t, U0, Esq_temp):

    N = len(t)-1   # Nº de intervalos   
    Nv = len(U0)    # Nº de variables

    U = zeros((Nv, N+1))

    U[:,0] = U0

    dt = t[1] - t[0]

    if Esq_temp == LF:
        
        U[:,1] = U[:,0] + dt*Phy(U[:,0], t[0])
        
        for n in range(1,N):

#            dt = t[n+1] - t[n]  # Para casos donde dt =/ cte

            U1 = U[:, n-1]
            U2 = U[:, n]

            U[:, n+1] = LF(U1, U2, t[n], dt, Phy)

    else:
        for n in range(N):

#            dt = t[n+1] - t[n]  # Para casos donde dt =/ cte

            U[:,n+1] = Esq_temp(U[:,n], dt, t[n], Phy)

    return U



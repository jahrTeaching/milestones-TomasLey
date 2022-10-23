from numpy import zeros, linspace

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

    for n in range(N):

        dt = t[n+1] - t[n]  # Para casos donde dt =/ cte
        U[:,n+1] = Esq_temp(U[:,n], t[n], dt, Phy)

    return U



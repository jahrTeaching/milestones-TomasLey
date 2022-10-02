from numpy import zeros, linspace

def Cauchy( Phy , t, U0, Esq_temp):

    # t es vector de todos los valores temporales que se van a estudiar donde t[i] = i*dt

    N = len(t)-1   # Nº de intervalos   
    Nv = len(U0)    # Nº de variables

    U = zeros((N+1, Nv))

    U[0,:] = U0

    for n in range(N):

        dt = t[n+1] - t[n]  # Para casos donde dt =/ cte
        U[n+1,:] = Esq_temp(U[n,:], t[n], dt, Phy)

    return U



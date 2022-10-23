from numpy import zeros, linspace, array, sqrt
from scipy.optimize import newton

"""
    Inputs:
        U: Vector de estado en tn
        t: tn
        dt: diferencial de tiempo
        F: Funcion dU/dT = F (En este caso órbita de Kepler)

    Outputs:
        U: Vector de estado en tn+1
"""

def Euler(U, t, dt, F):

    return U + dt * F(U, t)

def RK4(U, t, dt, F):


        k1 = dt * F(U, t)
        k2 = dt * F(U + k1/2, t + dt/2)
        k3 = dt * F(U + k2/2, t + dt/2)
        k4 = dt * F(U + k3, t + dt)

        k = (k1+2*k2+2*k3+k4)/6

        return U + k

def Euler_Inv(U, t, dt, F):

    def Residuo_EI(X):

        res = X - U - dt*F(X, t)

        return res

    return newton(func = Residuo_EI, x0 = U, tol = 1e-04, maxiter = 1000)

def CN(U, t, dt, F):

    def Residuo_CN(X):

        res = X - U - dt * (F(X, t+dt) + F(U, t)) / 2

        return res

    return newton(func = Residuo_CN, x0 = U)


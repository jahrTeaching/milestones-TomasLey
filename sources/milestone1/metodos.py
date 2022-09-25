from numpy import zeros, linspace, array, sqrt


def fun(U):
 
    f1 = U[2]
    f2 = U[3]
    f3 = -U[0] / ((U[0]**2 + U[1]**2)**(3/2))
    f4 = -U[1] / ((U[0]**2 + U[1]**2)**(3/2))
    Fn = [f1,f2,f3,f4]
    Fn = array(Fn)

    return Fn

def Euler(UE,N,dt):

    for i in range(N):

        Fn = fun(UE[:,i])
        UE[:,i+1] = UE[:,i] + dt*Fn

    return UE

def RK4(Urk,N,dt):

    for i in range(N):
        k1 = dt * fun(Urk[:,i])
        k2 = dt * fun(Urk[:,i]+k1/2)
        k3 = dt * fun(Urk[:,i]+k2/2)
        k4 = dt * fun(Urk[:,i]+k3)

        k = (k1+2*k2+2*k3+k4)/6
        Urk[:,i+1] = Urk[:,i] + k

    return Urk


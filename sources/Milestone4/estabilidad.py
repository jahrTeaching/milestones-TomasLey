from numpy import meshgrid, array, zeros, size, absolute, sqrt, linspace, float64
from esquemas_temporales import Euler, RK4, Euler_Inv, CN, LF

## Absolute Stability Region

def region_estabilidad(Temporal_Scheme, N, x0, xf, y0, yf):

    x, y = linspace(x0, xf, N), linspace(y0, yf, N)
    #rho =  zeros( (N, N),  dtype=float64)
    rho =  zeros( (N, N))

    for i in range(N): 
      for j in range(N):

          w = complex(x[i], y[j])
          r = Temporal_Scheme( 1., 1., 0., lambda u, t: w*u )
          rho[i, j] = abs(r) 

    return rho, x, y 



def region_estabilidad_trampa(x,y,Temporal_Scheme):

    N = size(x)
    Z = zeros([N,N],dtype=complex)

    for i in range(N):

        for j in range(N):

            Z[N-1-j,i] = complex(x[i],y[j])

    return absolute(array(polinomio(Temporal_Scheme, Z)))


def polinomio(Temporal_Scheme, w):

    if Temporal_Scheme == Euler:
        return 1 + w

    elif Temporal_Scheme == Euler_Inv:
        return 1/(1-w)

    elif Temporal_Scheme == CN:
        return (1+w/2)/(1-w/2)

    elif Temporal_Scheme == RK4:
        return 1 + w + (w**2)/2 + (w**3)/6 + (w**4)/(4*3*2)

    elif Temporal_Scheme == LF:
        return sqrt(1)


from numpy import array, linspace, zeros, size, log10, sqrt
from numpy.linalg import norm
from scipy.stats import linregress
import matplotlib.pyplot as plt

"""
    Inputs:
        Problema: Cauchy
        Phy: Problema físico a estudiar
        t: vector de todos los valores temporales que se van a estudiar donde t[i] = t0 + i*dt 
        U0: Condiciones iniciales
        Esq_temp: Esquema temporal utilizado para integrar el problema de Cauchy

    Outputs:
        Err
"""

 ## Funciones de representación de resultados

def ploterr(t,Err):

    ErrNorm = zeros(len(t))
    ErrNorm[:] = sqrt(Err[0,:]**2 + Err[1,:]**2)

    plt.subplots(figsize=(5,5))
    plt.plot(t,Err[0,:], label='Error en x')
    plt.plot(t,Err[1,:], label='Error en y')
    plt.plot(t,ErrNorm[:], label='Error total')
    plt.title('Error para dt = {}'.format(round(t[1]-t[0], 4)))
    plt.legend()
    

def plotord(logN,logErr,q):

    plt.subplots(figsize=(5,5))
    plt.plot(logN,logErr)
    plt.plot(logN,q.intercept + q.slope*logN)
    plt.title('Orden del esquema = {}'.format(round(q.slope, 2)))

## Funciones de error

def Richardson(problema, Phy, t, U0, Esq_temp, q_order, iter):

    N = size(t)

    for i in range(iter):

        Err = zeros((len(U0),N))

        t1 = linspace(t[0], t[-1], N)
        t2 = linspace(t[0], t[-1], 2*N)

        Un = problema(Phy, t1, U0, Esq_temp)
        U2n = problema(Phy, t2, U0, Esq_temp)

        for j in range(N):
        
            Err[:,j] = (Un[:,j]-U2n[:,2*j])/(1-1/(2**q_order))

        ploterr(t1,Err)

        N = 10*N

    plt.show()

    return Err


def Convergence(problema, Phy, t, U0, Esq_temp, mpoints): # No se usa la función Richardson para calcular el error porque 
                                                           # el objetivo es obtener el q_order
    N = size(t)

    logErr = zeros(mpoints)
    logN = zeros(mpoints)

    for i in range(mpoints):

        t1 = linspace(t[0], t[-1], N)        
        t2 = linspace(t[0], t[-1], 2*N)

        Un = zeros((len(U0),len(t1)))
        Un2 = zeros((len(U0),len(t2)))

        Un = problema(Phy, t1, U0, Esq_temp)
        U2n = problema(Phy, t2, U0, Esq_temp)

        logErr[i] = log10(norm(Un[:,N-1]-U2n[:,2*N-1]))
        logN[i] = log10(N)
        
        N = 2*N

    q = linregress(logN, logErr)

    plotord(logN,logErr,q)

    plt.show()

    return logErr, logN, q


def Error_y_convergencia(problema, Phy, t, U0, Esq_temp, q_order, mpoints): # Combina ambas funciones para no calcular las soluciones del problema 2 veces

    N = size(t)

    logErr = zeros(mpoints)
    logN = zeros(mpoints)


    for i in range(mpoints):

        t1 = linspace(t[0], t[-1], N)        
        t2 = linspace(t[0], t[-1], 2*N)

        Un = zeros((len(U0),len(t1)))
        Un2 = zeros((len(U0),len(t2)))
        Err = zeros((len(U0),N))

        Un = problema(Phy, t1, U0, Esq_temp)
        U2n = problema(Phy, t2, U0, Esq_temp)

        for j in range(N):
        
            Err[:,j] = (Un[:,j]-U2n[:,2*j])/(1-1/(2**q_order))

        logErr[i] = log10(norm(Un[:,N-1]-U2n[:,2*N-1]))
        logN[i] = log10(N)

        #ploterr(t1,Err)
        
        N = 2*N

    q = linregress(logN, logErr)

    plotord(logN,logErr,q)

    plt.show()





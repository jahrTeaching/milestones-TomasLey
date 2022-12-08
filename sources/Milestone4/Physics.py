from numpy import zeros, linspace, array, sqrt


def Kepler(U, t):      # U vector de estado en un tiempo tn 
    
    d = (U[0]**2 + U[1]**2)**(3/2)

    f1 = U[2]
    f2 = U[3]
    f3 = -U[0] / d
    f4 = -U[1] / d

    return array([f1,f2,f3,f4])

def Oscilator(U,t):

    return array([U[1], -U[0]])

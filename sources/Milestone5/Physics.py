from numpy import zeros, linspace, array, sqrt, reshape, concatenate
from numpy.linalg import norm

"""
    Inputs:
        U: vector de estado en tiempo t

    Outputs:
        F: Funcion de integración
"""

def Kepler(U, t):      
    
    d = (U[0]**2 + U[1]**2)**(3/2)

    f1 = U[2]
    f2 = U[3]
    f3 = -U[0] / d
    f4 = -U[1] / d

    return array([f1,f2,f3,f4])

def Oscilator(U,t):                     # x'' = -x

    return array([U[1], -U[0]])

def F_N_body(U,t):

    (Nb,Nc) = (4,3)    

    Us = reshape(U,(Nb,Nc,2))           
    F  = zeros(len(U))
    Fs = reshape(F,(Nb,Nc,2))           

    r  = reshape(Us[:,:,0],(Nb,Nc))     # N-Body locations
    v  = reshape(Us[:,:,1],(Nb,Nc))     # N-Body velocities

    drdt = reshape(Fs[:,:,0], (Nb,Nc))
    dvdt = reshape(Fs[:,:,1], (Nb,Nc))

    dvdt[:,:] = 0                 

    for i in range(Nb):

        drdt[i,:] = v[i,:]

        for j in range(Nb):
            
            if j != i:                   # Different bodies atract eachother

                d = r[j,:] - r[i,:]
                dvdt[i,:] += d[:]/(norm(d)**3)  

    return F

def initial_conditions_Nbody(Nb):       # Function to determine initial conditions

    if Nb == 2:                     # 2 bodies
        r01 = array([1, 0, 0])
        v01 = array([0, -0.5, 0])

        r02 = array([-0.5, 0, 0])
        v02 = array([0, 0.5, 0])

        r = concatenate((r01,r02),axis=0)
        v = concatenate((v01,v02),axis=0)

    elif Nb == 3:                      # 4 bodies, 

        r01 = array([5, 0, 0])/11
        v01 = array([0, 1, 0])

        r02 = array([-2.5, 4.33012701892, 0])/11
        v02 = array([-sqrt(3)/2, -0.5, 0])

        r03 = array([-2.5, -4.33012701892, 0])/11
        v03 = array([sqrt(3)/2, -0.5, 0])

        r = concatenate((r01,r02,r03),axis=0)
        v = concatenate((v01,v02,v03),axis=0)

    elif Nb == 4:                     # 4 bodies

        r01 = array([2, 2, 1])
        v01 = array([-0.5, 0, 0])

        r02 = array([-2,2,-1])
        v02 = array([0,-0.5,0])

        r03 = array([-2, -2, 1])
        v03 = array([0.5, 0, 0])

        r04 = array([2, -2, -1])
        v04 = array([0, 0.5, 0])

        r = concatenate((r01,r02,r03,r04),axis=0)
        v = concatenate((v01,v02,v03,v04),axis=0)

    return r, v

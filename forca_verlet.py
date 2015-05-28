import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import sys
import time


def distance(x,y,X,Y,l2): 

       
    #calculates the distances
    distx = x[:,np.newaxis] - x[np.newaxis,:] 
    disty = y[:,np.newaxis] - y[np.newaxis,:]
    
    #Periodic condition of contour
    distx = distx -X*np.rint(distx/X)
    disty = disty -Y*np.rint(disty/Y) 
    
    R2 = distx**2 + disty**2

    return R2


def forcas_verlet(x,y,xold,yold,X,Y,l2,R2,rv,rc): 

    #initialize matrices for force and distances
   
    N = len(x)
    f = np.zeros(shape=(N,N))
    V = np.zeros(shape=(N,N))

    #calculates the displacement

    desloc = np.sqrt((xold -x)**2 + (yold-y)**2)

    if (np.any(desloc)>(rv-rc)):

        R2 = distance(x,y,X,Y)
        #calculates verlet lists of idices
        s =  (np.sqrt(R2)<rv)
        vindi = np.nonzero(s)
        vlist = [np.nonzero(s[i])[0] for i in range(len(s))]

    
    else:
        
        for i in range(R2.shape[0]):

            distx = x[i] - x[vlist[i]]
            disty = y[i] - y[vlist[i]]
            
            #Periodic condition of contour
            distx = distx -X*np.rint(distx/X)
            disty = disty -Y*np.rint(disty/Y)

            R2[vindi] = distx**2 + disty**2
 
    
    
    V[vindi] = 4*((R2[vindi]**-6.)-(R2[vindi]**-3.)) 
    
    #calculates the force
    
    f[vindi] = 48*(np.sqrt(R2[vindi])**(-13.) - 0.5*np.sqrt(R2[vindi])**(-7.))
    
    #equals the lower triangle part of the array to the upper
    #f += f.T
    #V += V.T
    
    #separe the force on x and y components
    fx = f*distx
    fy = f*disty
    
    
    #sum them up to have the resulting force in each particle
    V = np.sum(V,axis=1)
    fx = np.sum(fx,axis=1)
    fy = np.sum(fy,axis=1)
    
    return fx,fy,V,R2

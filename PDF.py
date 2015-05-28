import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import sys
import time
import dymol



def gr(X,Y,R2):

    dr = 0.2
    rmax = np.sqrt(X**2 + Y**2)
    nbins = rmax/dr
    R = np.sqrt(R2)
    k = (R[0,:] >0)
    Hist,bord = np.histogram(R[0,:][k],nbins,range=(0,rmax)) 
    for i in range(1,R.shape[0]):
        k = (R[i,:] >0)
        hist,bord = np.histogram(R[i,:][k],nbins,range=(0,rmax))
        Hist +=hist



        c = (bord[1:]-bord[:-1])*0.5 + bord[0:-1]


        #normalizacao da area
        #areas dos aneis
        a = np.pi*(bord[1:]**2 -(bord[:-1])**2)
        
        hist = 2*Hist/a
        #normalizacao de densidade
        hist /= (R.shape[0]**2/(X*Y))
        return c,hist


#plt.plot(c,hist)
#plt.show()

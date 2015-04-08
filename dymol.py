import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import sys
import time


"""Molecular Dynamics code, for a set of initial conditions simulates
the dynamics of particles uniformly distibruted submited to a
potencial (Eg. Lennard Jones).  
Wroted by Anna Barbara 2015, 31 of March
"""

#-----------------------initial conditions------------------------------
def initial():

    X,Y = 32,32               #box's size
    dt = 0.001                #time interval
    N = 256                   #number of particles
    tmax = 4                 #final time of iteraction
    l2 = 25                    #minimum distance for the potencial
    
    
    """Uniformly distribution of positions"""
    
    nx = np.sqrt(N)            #particles in x
    ny = np.sqrt(N)
    
    #tile = repeat the hole array, reapet = repeat each element

    x = np.tile(np.linspace(0.5,X-0.5,nx),nx) 
    y = np.repeat(np.linspace(0.5,Y-0.5,ny),ny)
   
    
    vx = np.random.normal(0,0.1,N) #normal distribution of velocities
    vy = np.random.normal(0,0.1,N)
        
    return x,y,vx,vy,X,Y,l2,tmax,dt,N

x,y,vx,vy,X,Y,l2,tmax,dt,N = initial() #calling function

#------------------Compute the force between all the particles --------
def forcas(x,y,X,Y,l2): 

    #initialize matrices for force and distances

    N=len(x)
    f = np.zeros(shape=(N,N))
    V = np.zeros(shape=(N,N))

    #calculates the distances
    distx = x[:,np.newaxis] - x[np.newaxis,:] 
    disty = y[:,np.newaxis] - y[np.newaxis,:]
    
    #Periodic condition of contour
    distx = distx -X*np.rint(distx/X)
    disty = disty -Y*np.rint(disty/Y) 
    
    #compute de distaces
    R2 = distx*distx + disty*disty
    R2 = np.triu(R2,k=1)
    
    k =(R2>0.0)&(R2<l2)#condition for the minimum distance of
                       #interaction calculates potential
    
    V[k] = 4*((R2[k]**-6.)-(R2[k]**-3.)) 
    
    #calculates the force
    
    f[k] = 48*(np.sqrt(R2[k])**(-13.) - 0.5*np.sqrt(R2[k])**(-7.))
    
    #equals the lower triangle part of the array to the upper
    f += f.T
    V += V.T
    
    #separe the force on x and y components
    fx = f*distx
    fy = f*disty
    
    
    #sum them up to have the resulting force in each particle
    V = np.sum(V,axis=1)
    fx = np.sum(fx,axis=1)
    fy = np.sum(fy,axis=1)
    
    return fx,fy,V,R2


def integrate(x,y,vx,vy,fx,fy,dt):
    ax = fx
    ay = fy
    vx += ax*dt
    vy += ay*dt
    x += vx*dt 
    y += vy*dt
    return x,y

def period(x,y,X,Y):
    x[(x<=0.0)]= X+x[(x<=0.0)]
    y[(y<=0.0)]= Y+y[(y<=0.0)]
    x[(x>= X)] = x[(x>= X)] - X
    y[(y>= Y)] = y[(y>= Y)] - Y
    return x,y


def range_of_steps(s):
    x,y,vx,vy,X,Y,l2,tmax,dt,N = initial()
    for i in range(s):
        fx,fy,V,R2 = forcas(x,y,X,Y,l2)
        x,y = integrate(x,y,vx,vy,fx,fy,dt)
    return x,y,fx,fy,V,R2




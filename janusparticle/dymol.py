import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import sys
import time
import forcas
"""Molecular Dynamics code, for a set of initial conditions simulates
the dynamics of particles uniformly distibruted submited to a
potencial (Eg. Lennard Jones).  
Wroted by Anna Barbara 2015, 31 of March
"""

#-----------------------initial conditions------------------------------
def initial():

    X,Y = 32,32            #box's size
    dt = 0.001             #time interval
    N = 2              #number of particles
    tmax =15.              #final time of iteraction
    l2 = 2.                #minimum distance for the potencial
    
    
    """Uniformly distribution of positions"""
    
    nx = np.sqrt(N)            #particles in x
    ny = np.sqrt(N)
    
    #tile = repeat the hole array, reapet = repeat each element

    #x = np.tile(np.linspace(0.5,X-0.5,nx),nx) 
    #y = np.repeat(np.linspace(0.5,Y-0.5,ny),ny)
    x = np.array([15.,16.])#,17.,18.,19.])
    y = np.array([15.,16.])#,17.,18.,19.])
    
    vx = np.random.normal(0,0.005,N) #normal distribution of velocities
    vy = np.random.normal(0,0.005,N)
        
    return x,y,vx,vy,X,Y,l2,tmax,dt,N

x,y,vx,vy,X,Y,l2,tmax,dt,N = initial() #calling function

#------------------Compute the force between all the particles --------
def forcas(x,y,X,Y,l2): 

    #initialize matrices for force and distances

    N=len(x)
    f = np.zeros(shape=(N,N))
    v = np.zeros(shape=(N,N))

    #calculates the distances
    distx = x[:,np.newaxis] - x[np.newaxis,:] 
    disty = y[:,np.newaxis] - y[np.newaxis,:]
    
    #Periodic condition of contour
    distx = distx -X*np.rint(distx/X)
    disty = disty -Y*np.rint(disty/Y) 
    
    #compute de distaces
    R2 = distx*distx + disty*disty
    #R2 = np.triu(R2,k=1)
    
    k =(R2>0.0)&(R2<l2)#condition for the minimum distance of
                       #interaction calculates potential
    
    v[k] = 4*((R2[k]**-6.)-(R2[k]**-3.)) 
    
    #calculates the force
    
    f[k] = 48*(np.sqrt(R2[k])**(-13.) - 0.5*np.sqrt(R2[k])**(-7.))
    
    #F = f
    #V = v
    #equals the lower triangle part of the array to the upper
    #f += f.T
    #v += v.T
    
    #separe the force on x and y components
    Fx = f*distx
    Fy = f*disty
    
    
    #sum them up to have the resulting force in each particle
    v = np.sum(v,axis=1)
    fx = np.sum(Fx,axis=1)
    fy = np.sum(Fy,axis=1)
    
    return fx,fy,v,R2


#otimizations methods options

"""calculates the indexes of the particles that are in the circle defined by rv"""
def verletlist(R2,rv,rc):
    s =  (R2<rv*rv)&(R2>0.)
    vlist = [np.nonzero(s[i])[0] for i in range(len(R2))]

    return vlist

"""calculates the forces with the verlet method of otimization """

def forcas_verlet(x,y,xold,yold,X,Y,l2,R2,rv,rc,vlist): 

    desloc = np.sqrt((xold -x)**2 + (yold-y)**2)
    
    if (np.any(desloc)>(rv-rc)):
        p = np.c_[x,y]
        size = len(x)
        fx,fy,V,R2 = forcas.lennardjones(size,p,X,Y,l2)
        #calculates verlet lists of idices
        vlist = verletlist(R2,rv,rc)
        
        return fx,fy,V,R2,vlist

    else:
        f = np.zeros(shape=(N,N))
        V = np.zeros(shape=(N,N))
        Fx = np.zeros(shape=(N,N))
        Fy = np.zeros(shape=(N,N))
        for i in range(R2.shape[0]):

            distx = x[i] - x[vlist[i]]
            
            disty = y[i] - y[vlist[i]]
            #print i, vlist[i]
            #Periodic condition of contour
            distx = distx -X*np.rint(distx/X)
            
            disty = disty -Y*np.rint(disty/Y)
            
            R2[i,:][vlist[i]] = distx**2 + disty**2
            r2 = distx**2 + disty**2
            
           
            
            V[i,:][vlist[i]] = 4*((r2**-6.)-(r2**-3.))
            f[i,:][vlist[i]] = 48*(np.sqrt(r2)**(-13.) - 0.5*np.sqrt(r2)**(-7.))
            #separe the force on x and y components
            Fx[i,:][vlist[i]] = f[i,:][vlist[i]]*distx
            Fy[i,:][vlist[i]] = f[i,:][vlist[i]]*disty
            #sum them up to have the resulting force in each particle
        
        fx = np.sum(Fx,axis=1)
        fy = np.sum(Fy,axis=1)
        V = np.sum(V,axis=1)
        return fx,fy,V,R2,vlist
        



def integrate(x,y,vx,vy,fx,fy,dt):
    X = x.copy()
    Y = y.copy()
    ax = fx
    ay = fy
    vx += ax*dt
    vy += ay*dt
    X += vx*dt 
    Y += vy*dt
    return X,Y,vx,vy

def integrate_rot(thetax,thetay,w,tau,sig,dt):

    THETAx = thetax.copy()
    THETAy = thetay.copy()
    W = w.copy()
    I = sig**2/2
    alpha = tau/I
    W += alpha*dt
    THETAx += -W*thetay*dt
    THETAy += W*thetax*dt
    return THETAx,THETAy,W



def period(x,y,X,Y):
    x[(x<=0.0)]= X+x[(x<=0.0)]
    y[(y<=0.0)]= Y+y[(y<=0.0)]
    x[(x>= X)] = x[(x>= X)] - X
    y[(y>= Y)] = y[(y>= Y)] - Y
    return x,y


def range_of_steps(x,y,vx,vy,X,Y,l2,dt,s):
    for i in range(s):
        fx,fy,V,R2 = forcas(x,y,X,Y,l2)
        x,y,vx,vy = integrate(x,y,vx,vy,fx,fy,dt)
    return x,y,vx,vy,fx,fy,V,R2


#0.1, 1, 10
#Testar o mean square displacement
def termo_andersen(vx,vy,dt,Temp,nu):
    if (np.random.uniform(0,1) < nu*dt):
        
        sig = np.sqrt(Temp)
        vx = np.random.normal(0,sig,len(vx))
        vy = np.random.normal(0,sig,len(vy))
                              
    return vx,vy

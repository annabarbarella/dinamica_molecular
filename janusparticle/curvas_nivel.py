import numpy as np
from matplotlib import pyplot as plt
import sys
import time
import dymol
import forcas
from matplotlib import animation
import PDF
#initial conditions-------------------------------------------------
x,y,vx,vy,X,Y,l2,tmax,dt,N = dymol.initial()

theta = 2*np.pi*np.random.rand(len(x))
w = np.ones(len(theta))*0.01
sig = 0.5
a = 1.5/sig
T = 300
kb = 1.38e-5 #cte de boltsmann micrometro X miligramas /s2 X K
C = 2*T*kb
#C = 0.001
mass = 0.02 #miligramas
qx = np.cos(theta)
qy = np.sin(theta)
#qx = np.array([-1.,1.0,1.0]) #np.cos(theta)
#qy = np.array([0.0,0.0,1.0])#/np.sqrt(2.) #np.cos(theta)#np.sin(theta)1

U = qx
V = qy
t = 0
size = len(x)
q = np.c_[qx,qy]
s = 0
C2 = 1.e-0
beta = 1e-2
g = 1e-0
psi = g*np.sqrt(kb*T*mass*beta/dt)*np.random.normal(0,1,1)

F = []
D = []
Ec = []
Er = []
Ep = []
time = []
I = mass*sig**2/2

for T in range(100,900,200):
    for i in range(1000):
    
        q = np.c_[qx,qy]
        p = np.c_[x,y]
        fx,fy,V,R2,tau = forcas.janusparticle(size,p,q,X,Y,l2,sig,a,C,C2)
        x,y,vx,vy = dymol.integrate_bump(x,y,vx,vy,fx,fy,beta,psi,mass,dt)
        qx,qy, w = dymol.integrate_rot(qx,qy,w,tau,sig,mass,dt)
        x,y = dymol.period(x,y,X,Y)
        
        qx,qy = qx/np.sqrt(qx**2+qy**2),qy/np.sqrt(qx**2+qy**2)
        
        psi =  g*np.sqrt(kb*T*mass*beta/dt)*np.random.normal(0,1,1)
    
    
        #Ec.append(np.sum(mass*0.5*(vx*vx+vy*vy)))
        #Ep.append(np.sum(V))
        #Er.append(np.sum(I*w*w*0.5))
        #time.append(i*dt)
        #rscl = 1
        #dscl = 1
        #rbin = (np.max(x)-np.min(x))/rscl
        #dbin = (np.max(y)-np.min(y))/dscl

    c,hist = PDF.gr(X,Y,R2)
    plt.plot(c,hist,label="T="+str(T))
plt.xlabel(r"distance between particles $\mu m$")
plt.ylabel("Pair distribution function g(r)")
plt.legend()
plt.show()
exit()
f, ax = plt.subplots(2, sharex=True)

plt.suptitle("Energy of the system")
ax[0].plot(time,Ec,'g.-',ms=1,label="Kinetic energy")
ax[0].plot(time,Ep,'m.-',ms=1,label="Potential energy")
ax[0].legend(loc=2)
ax[1].plot(time,Er,'r.-',ms=1,label="Rotational energy")
ax[0].set_xlabel(r"time passed ($seconds$)")
ax[0].set_ylabel(r'energy ($ mg \mu m^2 /s^2 $) ')
ax[1].set_xlabel(r"time passed ($seconds$)")
ax[1].set_ylabel(r'energy ($ mg \mu m^2 /s^2 $) ')
ax[1].legend(loc=2)
plt.show()



hist,xlim,ylim = np.histogram2d(x,y,bins=(rbin,dbin))
#print ralim
xlim = xlim[:-1] #+ rscl*0.5
ylim = ylim[:-1] #+ dscl*0.5

imgplot = plt.imshow(hist, extent=[0, X, 0, Y],aspect="auto")         
plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")
#imgplot.set_cmap('hot')

#CS = plt.contourf(xlim, ylim, hist.T)
#plt.savefig("results/CMD")
#imgplot.set_clim(0.0,110.0)
#imgplot.set_interpolation('bicubic')
plt.show()


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import sys
import time
import dymol


x,y,vx,vy,X,Y,l2,tmax,dt,N = dymol.initial()
dr = 0.2
rmax = np.sqrt(X**2 + Y**2)
nbins = rmax/dr
t = 0.
while(t<tmax):
    fx,fy,V,R2 = dymol.forcas(x,y,X,Y,l2)
    x,y = dymol.integrate(x,y,vx,vy,fx,fy,dt)
    t += dt


xi= np.array(x)
yi = np.array(y)
#DELTAR = np.arange(t,tmax,dt)
DELTAR = []
#DELTA.APPEN
tempo = []
t = 0.
while(t<tmax):
    fx,fy,V,R2 = dymol.forcas(x,y,X,Y,l2)
    x,y = dymol.integrate(x,y,vx,vy,fx,fy,dt)
    
    deltar = ((x-xi)**2 + (y-yi)**2)
    DELTAR.append(np.mean(deltar)) 
    tempo.append(t)
    t = t+dt

plt.plot(tempo,DELTAR)
plt.show()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import sys
import time
import dymol
import PDF
import forcas

x,y,vx,vy,X,Y,l2,tmax,dt,N = dymol.initial()
size = len(x)
p = np.c_[x,y]

fx,fy,V,R2 = forcas.lennardjones(size,p,X,Y,l2)

t = 0.
V1 = []
E1 = []
K1 = []
start = time.time()
tempo = []
Temp = 1.
nu = 1.

while(t<tmax):
    
    fx,fy,V,R2 = forcas.lennardjones(size,p,X,Y,l2)
    
       
    x,y,vx,vy = dymol.integrate(x,y,vx,vy,fx,fy,dt)
    
    #vx, vy = dymol.termo_andersen(vx,vy,dt,Temp,nu)
    p = np.c_[x,y]
    k = np.sum(0.5*(vx*vx +vy*vy))
    t += dt
    tempo.append(t)
    #E1.append(v+k) 
    K1.append(k)
    V1.append(np.sum(V))
    #print t
    #print '\r',v,
    #sys.stdout.flush()


print 'time do loop fortran: ', time.time()-start
tempo = np.array(tempo)

plt.plot(tempo,K1,'g-',ms=1,label="Kinetic energy fortran")
plt.plot(tempo,V1,'r-',ms=1,label="Potential energy fortran")
plt.show()

#Pair distributions function
c,hist = PDF.gr(X,Y,R2)
hist1 = hist + 0.0
#mean on time for Pair distribution function

for i in range(0,100):
    c,hist = PDF.gr(X,Y,R2)
    hist1 = hist1 + hist 
    for j in range(0,100):
        fx,fy,V,R2 = forcas.lennardjones(size,p,X,Y,l2)
        x,y,vx,vy = dymol.integrate(x,y,vx,vy,fx,fy,dt)
        p = np.c_[x,y]

HIST = hist1/100.
plt.plot(c,HIST)
plt.title("Pair distribution function density = 1.")
plt.xlabel("dist")
plt.ylabel("g(dist)")
plt.show()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import sys
import time
import dymol
import PDF
import forcas

#verlet radius
rc,rv = 2.,3.

#initialize positions
x,y,vx,vy,X,Y,l2,tmax,dt,N = dymol.initial()
 
#initialize forces and distances and verlet lists
fx,fy,V,R2 = dymol.forcas(x,y,X,Y,l2)
xnew,ynew,vx,vy = dymol.integrate(x,y,vx,vy,fx,fy,dt)
vlist = dymol.verletlist(R2,rv,rc)

t = 0.

E1 = []
K1 = []
V1 = []
start = time.time()
tempo = []
count = 0.
Hist = 0.
temp = 2.0#np.sum(0.5*(vx*vx +vy*vy))/N
while(t<tmax):
    
    fx,fy,V,R2,vlist = dymol.forcas_verlet(xnew,ynew,x,y,X,Y,l2,R2,rv,rc,vlist)
    x = np.array(xnew)
    y = np.array(ynew)
       
    xnew,ynew,vx,vy = dymol.integrate(x,y,vx,vy,fx,fy,dt)
        
    #vx,vy = dymol.termo_andersen(vx,vy,dt,temp,1.)
    
    k = np.sum(0.5*(vx*vx +vy*vy))
    t += dt
    tempo.append(t)
    
    K1.append(k)
    V1.append(np.sum(V))

    #vx,vy = dymol.termo_andersen(vx,vy,dt,temp,0.1)
        

    #print t
    


print 'time do loop: ', time.time()-start
    

tempo = np.array(tempo)


plt.plot(tempo,K1,'g-',ms=1,label="Kinetic energy")
plt.plot(tempo,V1,'r-',ms=1,label="Potential energy")
plt.legend()
plt.show()


#Probability density function
#c,hist = PDF.gr(X,Y,R2)
#hist = Hist

#plt.plot(c,hist,label="Probability density function")
#plt.show()

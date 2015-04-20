import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import sys
import time
import dymol

#condicoes iniciais
x,y,vx,vy,X,Y,l2,tmax,dt,N = dymol.initial()
dr = 0.2
rmax = np.sqrt(X**2 + Y**2)
t = 0.

#para o sistema entrar no equilibrio
while(t<tmax):
    fx,fy,V,R2 = dymol.forcas(x,y,X,Y,l2)
    x,y = dymol.integrate(x,y,vx,vy,fx,fy,dt)
    t += dt
deltam = 0.

gga = []
for j in range(100):
    
    x,y,fx,fy,V,R2 = dymol.range_of_steps(x,y,vx,vy,X,Y,l2,dt,10)
    gga.append(np.array([x,y]))


    


deltam = deltam/100.
tempo = np.arange(4,9,0.05)
print len(deltam),len(tempo)
plt.plot(tempo,deltam,'g-',ms=1)
plt.xlabel("tempo")
plt.ylabel("r^2")
plt.show()

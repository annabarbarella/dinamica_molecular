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

for j in range(100):
    fx,fy,V,R2 = dymol.forcas(x,y,X,Y,l2)
    x,y = dymol.integrate(x,y,vx,vy,fx,fy,dt)
    xi= np.array(x)
    yi = np.array(y)
    x,y,fx,fy,V,R2 = dymol.range_of_steps(50)
    DELTAR = []
    tempo = []
    t = 0.
    xf,yf = xi,yi
    for i in range(100):
        fx,fy,V,R2 = dymol.forcas(xf,yf,X,Y,l2)
        xf,yf = dymol.integrate(x,y,vx,vy,fx,fy,dt)
        deltar = ((xf-xi)**2 + (yf-yi)**2)
        DELTAR.append(np.mean(deltar)) 
        #tempo.append(t)
        #t = t+dt
    DELTAR = np.array(DELTAR)
    deltam += DELTAR
    DELTAR = []

deltam = deltam/100.
tempo = np.arange(4,4.1,0.001)
print len(deltam),len(tempo)
plt.plot(tempo,deltam,'g-',ms=1)
plt.xlabel("tempo")
plt.ylabel("r^2")
plt.show()

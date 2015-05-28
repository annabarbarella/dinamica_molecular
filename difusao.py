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
    x,y,vx,vy = dymol.integrate(x,y,vx,vy,fx,fy,dt)
    t += dt
deltam = 0.
deltamz = 0.

for j in range(100):
    fx,fy,V,R2 = dymol.forcas(x,y,X,Y,l2)
    x,y,vx,vy = dymol.integrate(x,y,vx,vy,fx,fy,dt)
    xi= np.array(x)
    yi = np.array(y)
    vxi = np.array(vx)
    vyi = np.array(vy)
    #x,y,vx,vy,fx,fy,V,R2 = dymol.range_of_steps(x,y,vx,vy,X,Y,l2,dt,10)
    DELTAR = []
    DELTAZ = []
    tempo = []
    t = 0.
    xf,yf,vxf,vyf = np.array(xi),np.array(yi),np.array(vxi),np.array(vyi)
    for i in range(100):
        
        fx,fy,V,R2 = dymol.forcas(xf,yf,X,Y,l2)
        xf,yf,vxf,vyf = dymol.integrate(xf,yf,vxf,vyf,fx,fy,dt)
        deltar = ((xf-xi)**2 + (yf-yi)**2)
        deltaz = (vxf*vxi) + (vyf*vyi)
        DELTAR.append(np.mean(deltar))
        DELTAZ.append(np.mean(deltaz))
        xf,yf,vxf,vyf,fx,fy,V,R2 = dymol.range_of_steps(xf,yf,vxf,vyf,X,Y,l2,dt,50)
           #t = t+dt
    DELTAR = np.array(DELTAR)
    DELTAZ = np.array(DELTAZ)
    deltam += DELTAR
    deltamz += DELTAZ
    DELTAR = []
    DELTAZ= []

deltam = deltam/100.
deltamz = deltamz/100.
tempo = np.arange(15,20,0.05)

print len(deltam),len(tempo)

plt.title("Difusion")
plt.text(60, .025,  "$\rho = 1/4$")
plt.plot(tempo,deltam,'g-',ms=1)
plt.xlabel("tempo")
plt.ylabel("r^2")
plt.show()

plt.title("Difusion")
plt.text(60, .025,  "$\rho = 1/4$")
plt.plot(tempo,deltamz,'r-',ms=1)
plt.xlabel("tempo")
plt.ylabel("Z")
plt.show()

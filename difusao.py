import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import sys
import time
import dymol
import forcas
#condicoes iniciais
x,y,vx,vy,X,Y,l2,tmax,dt,N = dymol.initial()
size = len(x)
p = np.c_[x,y]
dr = 0.2
rmax = np.sqrt(X**2 + Y**2)
t = 0.

#para o sistema entrar no equilibrio
while(t<tmax):
    fx,fy,V,R2 = forcas.lennardjones(size,p,X,Y,l2)

    x,y,vx,vy = dymol.integrate(x,y,vx,vy,fx,fy,dt)
    p = np.c_[x,y]
    t += dt
deltam = 0.
deltamz = 0.

for j in range(50):
    #for k in range(20):
    fx,fy,V,R2 = forcas.lennardjones(size,p,X,Y,l2)
    x,y,vx,vy = dymol.integrate(x,y,vx,vy,fx,fy,dt)
    p = np.c_[x,y]
        
    xi = np.array(x)
    yi = np.array(y)
    pi = np.array(p)
    vxi = np.array(vx)
    vyi = np.array(vy)
    DELTAR = []
    DELTAZ = []
    tempo = []
    t = 0.

    xf = x.copy()
    yf = y.copy()
    pf= np.c_[xf,yf]
    vxf = vx.copy()
    vyf = vy.copy()

    for k in range(50):
        fx,fy,V,R2 = forcas.lennardjones(size,p,X,Y,l2)
        x,y,vx,vy = dymol.integrate(x,y,vx,vy,fx,fy,dt)
        p = np.c_[x,y]
    
    print j    
    for i in range(50):
        for k in range(50):
            fx,fy,V,R2 =  forcas.lennardjones(size,pf,X,Y,l2)
            xf,yf,vxf,vyf = dymol.integrate(xf,yf,vxf,vyf,fx,fy,dt)
            pf = np.c_[xf,yf]
        
        deltar = ((xf-xi)**2 + (yf-yi)**2)
        deltaz = (vxf*vxi) + (vyf*vyi)
        DELTAR.append(np.mean(deltar))
        DELTAZ.append(np.mean(deltaz))
        
    DELTAR = np.array(DELTAR)
    DELTAZ = np.array(DELTAZ)
    deltam += DELTAR
    deltamz += DELTAZ
    DELTAR = []
    DELTAZ= []







deltam = deltam/20.
deltamz = deltamz/20.

tempo = np.arange(15,17.5,0.05)

print len(deltam),len(tempo)

plt.title("Difusion")
#plt.text(60, .025,  "$\rho = 1/4$")
plt.plot(tempo,deltam,'g-',ms=1)
plt.xlabel("tempo")
plt.ylabel("r^2")
plt.show()

plt.title("Difusion")
#plt.text(60, .025,  "$\rho = 1/4$")
plt.plot(tempo,deltamz,'r-',ms=1)
plt.xlabel("tempo")
plt.ylabel("Z")
plt.show()

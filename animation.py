import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import sys
import time
import dymol

x,y,vx,vy,X,Y,l2,tmax,dt,N = dymol.initial()
#======================== defining plots ====================================


# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(-2, X+2), ylim=(-2, Y+2))
line, = ax.plot([], [],'mo', ms=5)
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)


# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text
    

#============================ def animation =======================

t = 0.0 # initial time
#dt = 0.01 # intervals of integration
frame = int(tmax/dt)

#while(t<tmax):
def animate(i):
    start = time.time()
    global t,x,y,vx,vy,X,Y,l2,tmax,dt,N
    fx,fy,V,R2 = dymol.forcas(x,y,X,Y,l2)
    x,y,vx,vy = dymol.integrate(x,y,vx,vy,fx,fy,dt)
    #x,y = period(x,y,X,Y)
    t = t+dt
    print np.sum((vx**2 + vy**2)*0.5) #+sum(V)
    
    
    line.set_data(x, y)
    
    time_text.set_text(time_template%(i*dt))
    #print 'time do loop: ', time.time()-start
    return line, time_text

#-------------------------calling animation -----------------------------------
anim = animation.FuncAnimation(fig, animate, init_func=init,frames=frame, interval=25, blit=True)
plt.show()
#------------------------------------------------------------------------------


#-------------------------calling the simulation without animation-------------
t = 0
E1 = []
K1 = []
V1 = []
tempo=[]
x,y,vx,vy,X,Y,l2,tmax,dt,N = dymol.initial()
start = time.time()
while(t<tmax):
    fx,fy,V,R2 = dymol.forcas(x,y,X,Y,l2)
    x,y,vx,vy = dymol.integrate(x,y,vx,vy,fx,fy,dt)
    t = t+dt
    k = np.sum(0.5*(vx*vx +vy*vy))
    v = np.sum(V)
    tempo.append(t)
    E1.append(v+k) 
    K1.append(k)
    V1.append(v)
    print v
    #print '\r',v,
    #sys.stdout.flush()


print 'time do loop: ', time.time()-start
    
E1 = np.array(E1)
tempo = np.array(tempo)

plt.plot(tempo,E1,'m-',ms=1,label="Total energy")
plt.plot(tempo,K1,'g-',ms=1,label="Kinetic energy")
plt.plot(tempo,V1,'r-',ms=1,label="Potential energy")
plt.legend()
plt.show()


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
#verlet radius
rc,rv = 0.,1.2
#initialize positions
x,y,vx,vy,X,Y,l2,tmax,dt,N = dymol.initial()
#while(t<tmax):
#initialize forces and distances and verlet lists
fx,fy,V,R2 = dymol.forcas(x,y,X,Y,l2)
xnew,ynew,vx,vy = dymol.integrate(x,y,vx,vy,fx,fy,dt)
vlist = dymol.verletlist(R2,rv,rc)

t = 0.


def animate(i):
    start = time.time()
    global x,y,vx,vy,X,Y,l2,tmax,dt,N,fx,fy,V,R2,xnew,ynew,vlist,t
    fx,fy,V,R2,vlist = dymol.forcas_verlet(xnew,ynew,x,y,X,Y,l2,R2,rv,rc,vlist)
    x = np.array(xnew)
    y = np.array(ynew)
       
    xnew,ynew,vx,vy = dymol.integrate(x,y,vx,vy,fx,fy,dt)
        
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

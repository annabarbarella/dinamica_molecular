import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import sys
import time
import dymol
import forcas

x,y,vx,vy,X,Y,l2,tmax,dt,N = dymol.initial()
#======================== defining plots ====================================


# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(-2, X+2), ylim=(-2, Y+2))
line, = ax.plot([], [],'mo', ms=15.)
lines = ax.quiver([],[],[],[],angles="xy")
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)


# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    lines.set_UVC([],[])
    lines.set_offsets([])
    time_text.set_text('')
    return line,lines,time_text
    

#============================ def animation =======================

t = 0.0 # initial time
#dt = 0.01 # intervals of integration
frame = int(tmax/dt)


#initialize positions
x,y,vx,vy,X,Y,l2,tmax,dt,N = dymol.initial()

t = 0.
vx,vy = 0*vx,0*vy
#vx[0] = 0.5

theta = 2*np.pi*np.random.rand(len(x))
w = np.zeros(len(theta))
sig = 0.1
a = 3./sig
C = 1
qx = np.cos(theta)
qy = np.sin(theta)
#qx,qy = np.ones(len(theta)),np.ones(len(theta))
#qx[:len(theta)/2.] = -qx[:len(theta)/2.]
#qy[:len(theta)/2.] = -qy[:len(theta)/2.]
q = np.c_[qx,qy]
size = len(x)
C2 = 1.
def animate(i):
    start = time.time()
    global x,y,vx,vy,X,Y,l2,tmax,dt,N,fx,fy,V,R2,xnew,ynew,vlist,t,q,a,C,sig,size,theta,tau,w,qx,qy,C2
    #fx,fy,V,R2 = dymol.forcas(x,y,X,Y,l2)    
    p = np.c_[x,y]
    #fx,fy,V,R2 = forcas.lennardjones(size,p,X,Y,l2)  
    for i in range(0,2):
        fx,fy,V,R2,tau = forcas.janusparticle(size,p,q,X,Y,l2,sig,a,C,C2)
        x,y,vx,vy = dymol.integrate(x,y,vx,vy,fx,fy,dt)
        qx,qy, w = dymol.integrate_rot(qx,qy,w,tau,sig,dt)
        #qx = np.cos(theta)
        #qy = np.sin(theta)
        q = np.c_[qx,qy]
        #F.append(fx[5]*fx[5]+fy[5]*fy[5])
        #D.append(R2[5])
        t = t+dt
        x,y = dymol.period(x,y,X,Y)
    #t = t+dt
    #print np.sum((vx**2 + vy**2)*0.5),w[5],qx[5],fx[5] #+sum(V)
    print fx
    
    #line.set_data(x, y)
    lines.set_UVC(qx,qy)
    lines.set_offsets(np.array([x, y]).T)
    time_text.set_text(time_template%(i*dt))
    #print 'time do loop: ', time.time()-start
    return line,time_text,




#-------------------------calling animation -----------------------------------
anim = animation.FuncAnimation(fig, animate, init_func=init,frames=frame, interval=25, blit=True)
plt.show()
#------------------------------------------------------------------------------



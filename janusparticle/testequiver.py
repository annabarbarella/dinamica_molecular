import numpy as np
from matplotlib import pyplot as plt
import sys
import time
import dymol
import forcas
from matplotlib import animation

x,y,vx,vy,X,Y,l2,tmax,dt,N = dymol.initial()
theta = 2*np.pi*np.random.rand(len(x))
#w = np.array([0.01,-0.01])#np.ones(len(theta))*0.01
w = np.ones(len(theta))*0.01
sig = 1.
a = 0.03/sig
C = 0.1e-0

qx = np.cos(theta)
qy = np.sin(theta)
#qx = np.array([-1.,1.0])#/np.sqrt(2.) #np.cos(theta)
#qy = np.array([1,0.0])#/np.sqrt(2.) #np.cos(theta)#np.sin(theta)
U = qx
V = qy

frame = int(tmax/dt)

fig, ax = plt.subplots(1,1)
Q = ax.quiver([], [], [], [], color='r', units='inches')
line, = ax.plot([], [],'mo', ms=15.)
ax.set_xlim(-2, X+2)
ax.set_ylim(-2, Y+2)
#Q.set_offsets(np.array([X, Y]))
t = 0
size = len(x)
q = np.c_[qx,qy]
s = 0
print "x",x
F = []
D = []
C2 = 10e-0
def update_quiver(num, Q, x, y,s,vx,vy,X,Y,l2,dt,N,t,q,a,C,C2,sig,size,theta,w,qx,qy):
    s = s+ num
    global xnew, ynew,vxnew,vynew,wnew,qxnew,qynew
    """updates the horizontal and vertical vector components by a
    fixed increment on each frame
    """
    
    if (num<1):
    
        p = np.c_[x,y]

    
        fxnew,fynew,Vnew,R2,tau = forcas.janusparticle(size,p,q,X,Y,l2,sig,a,C,C2)
        xnew,ynew,vxnew,vynew = dymol.integrate(x,y,vx,vy,fxnew,fynew,dt)
    
        qxnew,qynew, wnew = dymol.integrate_rot(qx,qy,w,tau,sig,dt)
        
        xnew,ynew = dymol.period(xnew,ynew,X,Y)
        print "qxnew",qxnew
    q = np.c_[qxnew,qynew]
    p = np.c_[xnew,xnew]
    #t = t+dt
    
    fxnew,fynew,Vnew,R2,tau = forcas.janusparticle(size,p,q,X,Y,l2,sig,a,C,C2)
    #xnew,ynew,vxnew,vynew = dymol.integrate(xnew,ynew,vxnew,vynew,fxnew,fynew,dt)
    
    qxnew,qynew, wnew = dymol.integrate_rot(qxnew,qynew,wnew,tau,sig,dt)
        
    xnew,ynew = dymol.period(xnew,ynew,X,Y)

    for i in range(50):

        q = np.c_[qxnew,qynew]
        p = np.c_[xnew,xnew]
        fxnew,fynew,Vnew,R2,tau = forcas.janusparticle(size,p,q,X,Y,l2,sig,a,C,C2)
        xnew,ynew,vxnew,vynew = dymol.integrate(xnew,ynew,vxnew,vynew,fxnew,fynew,dt)
    
        qxnew,qynew, wnew = dymol.integrate_rot(qxnew,qynew,wnew,tau,sig,dt)
        
        xnew,ynew = dymol.period(xnew,ynew,X,Y)
        #qxnew[1],qynew[1] = 1,0.0
        #xnew[1],ynew[1],ynew[0] = 16.,15.,15
    
    F.append(fxnew[0]**2+fynew[0]**2)
    D.append(R2[0,1])
        
    print "tau",tau,fxnew
    
    #qxnew,qynew = np.cos(num*0.1),np.sin( num*0.1)
    #U = np.cos( num*0.1)
    #V = np.sin( num*0.1)
    

    Q.set_UVC(qxnew,qynew)
    Q.set_offsets(np.array([xnew, ynew]).T)
    line.set_data(xnew,ynew)
    return Q,line,

# you need to set blit=False, or the first set of arrows never gets
# cleared on subsequent frames
anim = animation.FuncAnimation(fig, update_quiver, fargs=(Q, x, y,s,vx,vy,X,Y,l2,dt,N,t,q,a,C,C2,sig,size,theta,w,qx,qy),frames=frame,interval=10, blit=False)

plt.show()
F = np.array(F)
D = np.array(D)


#plt.plot(D[(D>0.05)],F[(D>0.05)],'k.',ms=5.)
plt.plot(D,F,'k.',ms=5.)
plt.show()

import numpy as np
from matplotlib import pyplot as plt
import sys
import time
import dymol
import forcas
from matplotlib import animation

#initial conditions-------------------------------------------------
x,y,vx,vy,X,Y,l2,tmax,dt,N = dymol.initial()

theta = 2*np.pi*np.random.rand(len(x))
w = np.ones(len(theta))*0.01
sig = 0.5
a = 1.5/sig
T = 300
kb = 1.38e-5 #cte de boltsmann micrometro X miligramas /s2 X K
C = 2*T*kb
#C = 0.001
mass = 0.02 #miligramas
qx = np.cos(theta)
qy = np.sin(theta)
#qx = np.array([-1.,1.0,1.0]) #np.cos(theta)
#qy = np.array([0.0,0.0,1.0])#/np.sqrt(2.) #np.cos(theta)#np.sin(theta)1

U = qx
V = qy
t = 0
size = len(x)
q = np.c_[qx,qy]
s = 0
C2 = 1.e-0
beta = 1e-2
g = 1e-10
psi = g*np.sqrt(kb*T*mass*beta/dt)*np.random.normal(0,1,1)

#setting frames and plots---------------------------------------------
frame = int(tmax/dt)
fig, ax = plt.subplots(1,1)
Q = ax.quiver([], [], [], [], color='m', units='inches',width=0.012)
area = 1500*0.16#1500 e uma area de raio 3
line= ax.scatter([], [], s=area, color="b",alpha=0.5)
ax.set_xlim(-2, X+2)
ax.set_ylim(-2, Y+2)
#---------------------------------------------------------------------

F = []
D = []

def update_quiver(num, Q, x, y,s,vx,vy,X,Y,l2,dt,N,t,q,a,C,C2,sig,size,theta,w,qx,qy):
    s = s+ num
    global xnew, ynew,vxnew,vynew,wnew,qxnew,qynew,beta,mass,psi
    """updates the horizontal and vertical vector components by a
    fixed increment on each frame
    """
    
    if (num<1):
    
        p = np.c_[x,y]

    
        fxnew,fynew,Vnew,R2,tau = forcas.janusparticle(size,p,q,X,Y,l2,sig,a,C,C2)
        xnew,ynew,vxnew,vynew = dymol.integrate_bump(x,y,vx,vy,fxnew,fynew,beta,psi,mass,dt)
    
        qxnew,qynew, wnew = dymol.integrate_rot(qx,qy,w,tau,sig,mass,dt)
        
        xnew,ynew = dymol.period(xnew,ynew,X,Y)
        

    q = np.c_[qxnew,qynew]
    p = np.c_[xnew,ynew]
    #t = t+dt
    
    fxnew,fynew,Vnew,R2,tau = forcas.janusparticle(size,p,q,X,Y,l2,sig,a,C,C2)
    xnew,ynew,vxnew,vynew = dymol.integrate_bump(xnew,ynew,vxnew,vynew,fxnew,fynew,beta,psi,mass,dt)
    
    qxnew,qynew, wnew = dymol.integrate_rot(qxnew,qynew,wnew,tau,sig,mass,dt)
        
    xnew,ynew = dymol.period(xnew,ynew,X,Y)
    psi =  g*np.sqrt(kb*T*mass*beta/dt)*np.random.normal(0,1,1)

    for i in range(49):
        
        q = np.c_[qxnew,qynew]
  
        p = np.c_[xnew,ynew]
        fxnew,fynew,Vnew,R2,tau = forcas.janusparticle(size,p,q,X,Y,l2,sig,a,C,C2)
        xnew,ynew,vxnew,vynew = dymol.integrate_bump(xnew,ynew,vxnew,vynew,fxnew,fynew,beta,psi,mass,dt)
    
        qxnew,qynew, wnew = dymol.integrate_rot(qxnew,qynew,wnew,tau,sig,mass,dt)
        
        xnew,ynew = dymol.period(xnew,ynew,X,Y)
        #qxnew[0],qynew[0] = 1.,0.0
        #xnew[1],ynew[1],ynew[0] = 16.,15.,15
        qxnew,qynew = qxnew/np.sqrt(qxnew**2+qynew**2),qynew/np.sqrt(qxnew**2+qynew**2)
        
        psi =  g*np.sqrt(kb*T*mass*beta/dt)*np.random.normal(0,1,1)
    
    F.append(fxnew[0]**2+fynew[0]**2)
    D.append(R2[0,1])
        
    print "tau",tau,fxnew,"vai"
    
    #qxnew,qynew = np.cos(num*0.1),np.sin( num*0.1)
    #U = np.cos( num*0.1)
    #V = np.sin( num*0.1)
    

    Q.set_UVC(qxnew,qynew)
    Q.set_offsets(np.array([xnew, ynew]).T)
    #line.set_data(xnew,ynew)
    line.set_offsets(np.array([xnew, ynew]).T)
    return Q,line

# you need to set blit=False, or the first set of arrows never gets
# cleared on subsequent frames
anim = animation.FuncAnimation(fig, update_quiver, fargs=(Q, x, y,s,vx,vy,X,Y,l2,dt,N,t,q,a,C,C2,sig,size,theta,w,qx,qy),frames=frame,interval=10, blit=False)
anim.save('image3.mp4',fps=15,bitrate=-1,extra_args=['-vcodec', 'libx264'])
plt.show()
#F = np.array(F)
#D = np.array(D)




#plt.plot(D[(D>0.05)],F[(D>0.05)],'k.',ms=5.)
#plt.plot(D,F,'k.',ms=5.)
#plt.show()

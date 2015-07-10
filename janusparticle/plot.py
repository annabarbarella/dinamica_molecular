from pylab import*


distx = np.arange(0.1,3.0,0.001)
disty =0*np.arange(0.1,3.0,0.001)

#parameters

#kb = 1.3e-23 #boltzmann 
#kb = 1.3e-0 #boltzmann 
#T = 273 #temprature
x=1
theta = 2*np.pi*np.random.rand()
w = 0.01
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
s = 0
C2 = 1.e-0
beta = 1e-2
g = 1e-10
#psi = g*np.sqrt(kb*T*mass*beta/dt)*np.random.normal(0,1,1)



#orientation

distqy = 0.
distqx = 2.


R2 = distx**2 + disty**2
R = np.sqrt(R2)

    
fx = +( distqx * distx + distqy * disty) * C * a * distx * exp(-a*(sqrt(R2)-sig))/(R2)**(3./2) \
     +2 * (distqx * distx + distqy * disty) * C * distx * exp(-a*(sqrt(R2)-sig))/(R2)**2 \
     - C * distqx * exp(-a*(sqrt(R2)-sig))/(R2)

fy = +( distqx * distx + distqy * disty) * C * a * disty * exp(-a*(sqrt(R2)-sig))/(R2)**(3./2) \
     +2 * (distqx * distx + distqy * disty) * C * disty * exp(-a*(sqrt(R2)-sig))/(R2)**2 \
     - C * distqy * exp(-a*(sqrt(R2)-sig))/(R2)
 

U =  (distqx*distx + distqy*disty)*C*exp(-a*(sqrt(R2)-sig))/(R2)
#fx=(U[1:]-U[:-1])/1E-3
print C*exp(-a*(sqrt(R2)-sig))/(R2) 

print sqrt(R2)-sig

fxrepul = 48 * e *( (sig/R)**13 )-0.5*( (sig/R)**7 )*distx
fyrepul = 48 * e *( (sig/R)**13 )-0.5*( (sig/R)**7 )*disty


f = np.sqrt(fx**2 + fy**2)
frepul = np.sqrt(fxrepul**2 + fyrepul**2)

dpi = 96
F, ax = plt.subplots(1, 2,figsize=(1600/dpi,500/dpi),dpi=dpi)


fx=-(U[1:]-U[:-1])/1E-3


ax[0].plot(R, U)#'g.--',ms=5.)

#plot(distx, fx2,'g.',ms=1.)
#plot(distx, fx2+fx,'r.',ms=1.)

suptitle("interaction between a pair of particles")
ax[0].set_xlabel("distance")
ax[0].set_ylabel("janus potential")
#ax[0].axis([0,0.4,0,7E7])

ax[1].plot(R[:-1], fx,'g',ms=5.)
ax[1].set_xlabel("distance")
ax[1].set_ylabel("janus force")
#ax[1].axis([0,0.4,0,7E7])


plt.show()


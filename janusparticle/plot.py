from pylab import*


distx = np.arange(0.00001,10.0,0.001)
disty = np.arange(0.00001,10.0,0.001)

#parameters

#kb = 1.3e-23 #boltzmann 
#kb = 1.3e-0 #boltzmann 
#T = 273 #temprature

C = 500.#6*kb*T #(5-7kbT)
sig = 1.e-3 #(particle radius)
a = 3/sig 
e = 1.#/(kb*T)

#orientation

distqy = 0.
distqx = 2.


R2 = distx**2 + disty**2
R = np.sqrt(R2)

    
fx = -( distqx * distx + distqy * disty) * C * a * distx * exp(-a*(sqrt(R2)-sig))/(R2)**(3./2) \
     -2 * (distqx * distx + distqy * disty) * C * distx * exp(-a*(sqrt(R2)-sig))/(R2)**2 \
     + C * distqx * exp(-a*(sqrt(R2)-sig))/(R2)

fy = -( distqx * distx + distqy * disty) * C * a * disty * exp(-a*(sqrt(R2)-sig))/(R2)**(3./2) \
     -2 * (distqx * distx + distqy * disty) * C * disty * exp(-a*(sqrt(R2)-sig))/(R2)**2 \
     + C * distqy * exp(-a*(sqrt(R2)-sig))/(R2)
 

U =  (distqx*distx + distqy*disty)*C*exp(-a*(sqrt(R2)-sig))/(R2) 
print C*exp(-a*(sqrt(R2)-sig))/(R2) 

fxrepul = 48 * e *( (sig/R)**13 )-0.5*( (sig/R)**7 )*distx
fyrepul = 48 * e *( (sig/R)**13 )-0.5*( (sig/R)**7 )*disty


f = np.sqrt(fx**2 + fy**2)
frepul = np.sqrt(fxrepul**2 + fyrepul**2)

dpi = 96
F, ax = plt.subplots(1, 2,figsize=(1600/dpi,500/dpi),dpi=dpi)



ax[0].plot(R, U)#'g.--',ms=5.)

#plot(distx, fx2,'g.',ms=1.)
#plot(distx, fx2+fx,'r.',ms=1.)

suptitle("interaction between a pair of particles")
ax[0].set_xlabel("distance")
ax[0].set_ylabel("janus potential")
ax[0].axis([0,0.4,0,7E7])

ax[1].plot(R, f,'g',ms=5.)
ax[1].set_xlabel("distance")
ax[1].set_ylabel("janus force")
ax[1].axis([0,0.4,0,7E7])


plt.show()


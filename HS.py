import numpy as np


def positions(x,y,vx,vy,dt):

    xnew = x + vx*dt
    ynew = y + yx*dt

    return xnew,ynew


def tempomim(x,y,vx,vy,raio):

    r = np.sqrt(x*x + y*y)
    
    v = np.sqrt(vx*vx + vy*vy)

    deltar = r[:,np.newaxis] - r[np.newaxis,:] 
    deltav = v[:,np.newaxis] - v[np.newaxis,:] 

    deltar[deltar==0] = 1e15
    deltav[deltav==0] = 1e15

    d = (deltav*deltar)**2 - (deltav*deltav)*(deltar*deltar -raio*raio)

    tempo = -(deltar*deltav + np.sqrt(d))/(deltav*deltav)

    return np.argmin(tempo)

    

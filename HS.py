import numpy as np


def positions(x,y,vx,vy,dt):

    xnew = x + vx*dt
    ynew = y + yx*dt

    return xnew,ynew


def tempomim(x,y,vx,vy,raio,X,Y):

    r = np.sqrt(x*x + y*y)

    distx = x[:,np.newaxis] - x[np.newaxis,:] 
    disty = y[:,np.newaxis] - y[np.newaxis,:]

    
    distx = distx -X*np.rint(distx/X)
    disty = disty -Y*np.rint(disty/Y) 
    
    
    distvx = vx[:,np.newaxis] - vx[np.newaxis,:] 
    distvy = vy[:,np.newaxis] - vy[np.newaxis,:]

    distx[distx == 0.]=1E15
    distx[disty == 0.]=1E15

    distx[distvx == 0.]=1E15
    distx[distvx == 0.]=1E15

    
    #deltar = np.sqrt(distx*distx + disty*disty)
    #deltav = np.sqrt(distvx*distvx + distvy*distvy)

    #deltar[deltar==0] = 1e15
    #deltav[deltav==0] = 1e15

    #d = (deltav*deltar)**2 - (deltav*deltav)*(deltar*deltar -raio*raio)
    d=(distx*distvx+disty*distvy)**2- (distvx*distvx+distvy*distvy)*((distvx*distvx+distvy*distvy)-raio*raio)
    d[d<0.]=1E15
    
    tempo=-(distx*distvx+disty*distvy+np.sqrt(d))/(distvx*distvx+distvy*distvy)
    
    #tempo = -(deltar*deltav + np.sqrt(d))/(deltav*deltav)
    tempo[tempo<=0.0] = 1e15    
    
    args = np.where(tempo == tempo.min())
    print args, np.min(tempo)
    exit()
    return args[0],args[1],np.min(tempo)

    

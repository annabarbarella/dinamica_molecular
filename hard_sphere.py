from pylab import*
import numpy as np
import dymol
import HS
import forcas

#initialize positions
x,y,vx,vy,X,Y,l2,tmax,dt,N = dymol.initial()
#calculates the par of particles which is gonna chock first
raio = 0.4



plt.plot(x,y,'g.',ms=3)
plt.quiver(x,y,vx,vy,angles="xy")
plt.show()
#atualizando posicoes e velocidade do par a se chocar

for i in range(0,100):
    arg,arg1,tmin = HS.tempomim(x,y,vx,vy,raio,X,Y)
    print tmin
    #atualizando posicoes e velocidade do par a se chocar
    
    x = x + vx*tmin
    y = y + vy*tmin

    
    vx[arg] = -vx[arg]
    vy[arg] = -vy[arg]
    vx[arg1] = -vx[arg1]
    vy[arg1] = -vy[arg1]

    
    plt.plot(x[arg],y[arg],'g.',ms=15)
    plt.plot(x[arg1],y[arg1],'g.',ms=15)
    plt.show()

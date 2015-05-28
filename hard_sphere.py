import numpy as np
import dymol
import HS


#initialize positions
x,y,vx,vy,X,Y,l2,tmax,dt,N = dymol.initial()
#calculates the par of particles which is gonna chock first
raio = 0.5

arg = HS.tempomim(x,y,vx,vy,raio)

print arg

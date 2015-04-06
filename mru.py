import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

#-----------------------condicoes iniciais------------------------------

def initial():

    X,Y = 16,16 #tamanho da caixa
    dt = 0.01 #delta tempo
    N = 25 #numero de particulas
    tmax = 20 #tempo maximo
    l2 = 9 #distancia maxima de calculo
    #distribuindo uniformemente
    nx = np.sqrt(N) #numero de particulas em x
    dl = X/nx #delta x e y
    x = np.zeros(N) # inicializa matrizes de posicao
    y = np.zeros(N)
    for i in range(N): #armazenas as posicoes iniciais uniformemente
                       #em um quadrado
        z = (i)*dl     
        x[i] = np.mod(z,X)+0.5
        y[i] = int(z/X)*dl + 0.5
    
    vx = np.random.uniform(0,1,N) #inicializa matrizes com velocidades iniciais zeros
    vxmed = np.mean(vx)
    
    vy = np.random.uniform(0,1,N)
    vymed = np.mean(vy)
    vx = vx -vxmed
    vy = vy -vymed
    return x,y,vx,vy,X,Y,l2,tmax,dt,N

x,y,vx,vy,X,Y,l2,tmax,dt,N = initial() #calling function


def integrate(x,y,vx,vy,dt): #integra as posicoes
    
    x += vx*dt #+= operacao mais ele mesmo
    y += vy*dt
    return x,y

x,y = integrate(x,y,vx,vy,dt)    

def period(x,y,X,Y):
    x[(x<=0.0)]= X+x[(x<=0.0)]
    y[(y<=0.0)]= Y+y[(y<=0.0)]
    x[(x>= X)] = x[(x>= X)] - X
    y[(y>= Y)] = y[(y>= Y)] - Y
    return x,y

x,y = period(x,y,X,Y)





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
    

#============================ calling animation =======================

t = 0.0 # initial time
#dt = 0.01 # intervals of integration
frame = int(tmax/dt)
x,y,vx,vy,X,Y,l2,tmax,dt,N = initial()
#while(t<tmax):
def animate(i):
    global t,x,y,vx,vy,X,Y,l2,tmax,dt,N
    
    x,y = period(x,y,X,Y)
    x,y = integrate(x,y,vx,vy,dt)
    
    t = t+dt
    T = (1/(2*N))*np.sum(vx*vx + vy*vy)
    
    line.set_data(x, y)
    
    time_text.set_text(time_template%(i*dt))
    return line, time_text

#def onClick(event):
#    global pause
#    pause ^= True
#fig.canvas.mpl_connect('button_press_event', onClick)

anim = animation.FuncAnimation(fig, animate, init_func=init,frames=frame, interval=25, blit=True)
plt.show()







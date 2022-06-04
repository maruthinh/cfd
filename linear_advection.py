import numpy as np
import matplotlib.pyplot as plt

#number of grid points
N=64

#Domain length
# L = 2*np.pi
L = 10

#Discretization of domain
u=np.linspace(-L, L, N)

#grid spacing
dx=u[1]-u[0]
#array to save initial condition
init = np.zeros(u.shape)
#array to save analytical condition
analytical = np.zeros(u.shape)

#different initial conditions
# init = np.exp(-(u-0.2)**2)
init = np.sin(u)

#solution initialization for upwind scheme
u_int_up = np.zeros(init.shape)
u_int_up[:] = init

#solution initialization for MacCormack scheme
u_int_mc = np.zeros(init.shape)
u_int_mc[:] = init

#Array to store upwind solution
u_up = np.zeros(u.shape)
#Array to store intermediate solution in solution in MacCormack method
u_mc_in = np.zeros(u.shape)
#Array to store solution in solution in MacCormack method
u_mc = np.zeros(u.shape)

t=0
a=0.5
period = 2*L/abs(a)
dt=0.8*(dx/abs(a))
t=0

#main iteration loop
while(t<period):
  #####MacCormack method######
  u_int_mc[0] = u_int_mc[N-3]
  u_int_mc[N-1] = u_int_mc[2]
  
  #Predictor step
  u_mc_in[:-1] = u_int_mc[:-1] - a*(dt/dx)*(u_int_mc[1:]-u_int_mc[:-1])
  
  #Corrector step
  u_mc[1:] = 0.5*(u_int_mc[1:]+u_mc_in[1:]) - 0.5*a*(dt/dx)*(u_mc_in[1:]-u_mc_in[:-1])
  
  #update solution
  u_int_mc[:]=u_mc[:]
  
  ##########Upwind method#####
  u_int_up[0]   = u_int_up[N-3]
  u_int_up[N-1] = u_int_up[2]
  
  if (a>=0):
      u_up[1:-1] = u_int_up[1:-1] - a*(dt/dx)*(u_int_up[1:-1]-u_int_up[:-2])
  elif(a<0):
      u_up[1:-1] = u_int_up[1:-1] - a*(dt/dx)*(u_int_up[2:]-u_int_up[1:-1])
      
  #update solution
  u_int_up[:]=u_up[:]
  
  t+=dt
  
  #analytical solution
  if(t<=0.25*period):
      t_temp = t
      # analytical = np.exp(-((u-a*t_temp)-0.2)**2)
      analytical = np.sin(u-a*t_temp)
  else:
      t_temp=t-period
      # analytical = np.exp(-((u-a*t_temp)-0.2)**2)
      analytical = np.sin(u-a*t_temp)
  
  #plotting
  plt.plot(u, init, label='Initial cond')
  plt.plot(u, analytical, label='Analytical')
  plt.plot(u[1:-1], u_up[1:-1], label='Upwind')
  plt.plot(u[1:-1], u_mc[1:-1], label='MacCormack')
  plt.xlabel('x', fontsize=16)
  plt.ylabel('$u(x)$', fontsize=16)
  plt.xticks(fontsize=16)
  plt.yticks(fontsize=16)
  plt.title('Linear advection: time: '+str("{:.2f}".format(t-dt)), fontsize=16)
  plt.legend(fontsize=12)
  plt.show() 
 
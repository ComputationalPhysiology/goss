#import os
from gotran import load_ode
from goss import *
import numpy as np

oscilator = jit(load_ode("oscilator"))


#solver1 = ExplicitEuler(oscilator)
solver1 = RL1(oscilator)

tstop = 10.
n_steps = 1000
time = np.linspace(0,tstop,n_steps+1)

dt = time[1]-time[0]
u0 = np.array([1.0,0.0])

u_exact = np.array([np.cos(time),np.sin(time)])
u1 = np.zeros_like(u_exact)
u2 = np.zeros_like(u_exact)
u1[:,0] = u0
u2[:,0] = u0

print(u1)
u = u0
for step in range(1,n_steps+1):
    #u = u1[step]
    t = time[step]
    solver1.forward(u, t, dt)
    u1[:,step] = u

np.savetxt('time.txt',time)
np.savetxt('solution.txt',u1)

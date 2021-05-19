import matplotlib.pyplot as plt
import numpy as np

time = np.loadtxt('time.txt')
u = np.loadtxt('solution.txt')

#print(u)

plt.plot(time, u[0,:])
plt.plot(time, u[1,:])
plt.show()

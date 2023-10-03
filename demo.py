import random
import numpy as np
import pylab as plt

Parent = []
Daughter = []
T_half_life = 28.8

lambdaa = np.log(2)/T_half_life # decay constant
num_nuclei_i = 1000 # initial num of radioactive nuclei

Np, Nd = num_nuclei_i, 0 
time = np.arange(0,1000,1)

for t in time:
    Parent.append(Np)
    Daughter.append(Nd)
    decay = 0 

    for n in range(Np):
        p = random.random()
        if p < lambdaa:
            decay +=1
    Np -= decay
    Nd += decay

plt.plot(Parent, linestyle='-', color='red', label='Parent Nucleus')
plt.plot(Daughter, linestyle='-', color='blue', label='Daughter Nucleus')
plt.legend()
plt.xlabel('Time (years)', fontsize=14)
plt.ylabel('N(t)', fontsize=14)
plt.show()

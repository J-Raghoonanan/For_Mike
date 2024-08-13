import numpy as np
import random
from matplotlib import pyplot as plt

def propose(m):
    random_num = random.random()
    if random_num <= 0.4:
        n = m+1
    elif 0.4 < random_num <= 0.8:
        n = m-1
    else:
        n=m
    return n

def update(n,m,T):
    if n<=m:
        if m==0:
            return m
        else:
            return n
    else:
        random_num = random.random()
        if random_num < np.exp(-1/T):
            return n
        else:
            return m
        
def generate(start,nmax,T):
    list_states = [start,]
    list_energies = [start+0.5,]
    m = start
    for i in range(n_max):
        proposal = propose(m)
        n = update(proposal, m, T)
        list_states.append(m)
        list_energies.append(n+0.5)
        m=n
    return (list_states,list_energies)

n_max = 2500
T=3
for n_start in [10,0]:
    states, energies = generate(n_start,n_max, T)
    plt.plot(energies, label = f"T={T},n_start={n_start}")

plt.ylabel(r'Energy')
plt.xlabel(r'n steps')
plt.title('Energy')
plt.legend()
plt.show()

fig = plt.gcf()
T=0.1
for n_start in [10,0]:
    states, energies = generate(n_start,n_max, T)
    plt.plot(energies, label = f"T={T},n_start={n_start}")

plt.ylabel(r'Energy')
plt.xlabel(r'n steps')
plt.title('Energy')
plt.legend()
plt.show()

def average_energy(n,old_ave,new_energy):
    return (n-1)*old_ave/n + new_energy/n

# Maxwell-boltzmann speed distribution

import numpy as np
from numpy import sqrt, pi,exp
import numpy as np
import matplotlib.pyplot as plt

kB = 1.38064852e-23
T = 1.5e7
m = 1.67e-27

v_rms = lambda T: sqrt(3*kB*T/m)
v = np.linspace(0.01*v_rms(T), 5*v_rms(T), 1000000)

def MB_speed(T, v, m):
    return (m/(2*pi*kB*T))**(3/2)*4*pi*(v**2)*exp(-m*v**2/(2*kB*T))

fusjon_fart = np.linspace(kritisk_hastighet,100*kritisk_hastighet,10000)

print("Andel som har større enn kritisk hastighet: ", sum(MB_speed(T,fusjon_fart,m)))

plt.plot(v, MB_speed(T, v, m))
plt.xlabel("Hastighet [m/s]")
plt.ylabel("Sannsynlighet")
plt.title("Maxwell-Boltzmann hastighetsfordeling")
plt.show()
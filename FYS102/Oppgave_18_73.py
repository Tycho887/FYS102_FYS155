"""Oppgave 18.73"""

import numpy as np

atm = 101300; L = 0.001
n = 1.0;      R = 8.314

P1 = 2.0*atm; P2 = 1.0*atm
V1 = 11.4*L;  V2 = 23*L

# Finner temperaturen og alpha
T = V1*P1/(n*R)
alpha = (P2-P1)/(V2-V1)

# Setter opp funksjonen for linje
def P_line(V):
    return alpha*(V-V2)+P2

# Definerer isotermisk kurve
def P_log(V):
    return n*R*T/V

# Finner arealet ved å løse integralen numerisk
dt = 0.001
V = np.arange(V1,V2,dt)
I = dt*np.sum(P_line(V)-P_log(V))

print(f"Arbeidet utført er: {I:.2g} Joule")


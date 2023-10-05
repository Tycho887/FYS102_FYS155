# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 08:25:48 2023

@author: Michael
"""

import numpy as np


mass = 0.267+0.250
T = 1.79
w_ = 2*np.pi/T
B = 0.101

def k(m,T):
    return 4*m*(np.pi)**2/(T)**2

def omega(O_,B,m):
    return np.sqrt(O_**2+(B/(2*m))**2)

w = omega(w_,B,mass)


# print(f'Gjennomsnitt er: {data2.mean():.3f}')
# print(f'Målt k-verdi: {k(mass,data2.mean()):.3f}')
print(f'Målt omega: {omega(w_,B,mass):.3f}')
print(f"delta_w: {w-w_:.4f}")
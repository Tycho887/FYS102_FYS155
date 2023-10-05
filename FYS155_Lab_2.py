# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 09:31:59 2023

@author: Michael
"""

import numpy as np

m = 0.208
r = 0.0142
g = 9.81
alpha1 = np.mean([2.48,2.39,2.36,2.40])
alpha2 = np.mean([2.09,2.08,2.09,2.08,2.09])

def I(a):
    return m*r*(g-r*a)/a

print(f"{I(alpha2):.4f} kgm^2")
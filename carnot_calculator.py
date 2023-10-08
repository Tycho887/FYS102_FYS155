import termoPy as TP
import matplotlib.pyplot as plt
import numpy as np

atm = 101325

gamma = 1.4
n = 1

P1 = 1*atm
P2 = 2*atm

T_hot = 300
T_cold = 200

V1 = n*TP.R*T_hot/P1



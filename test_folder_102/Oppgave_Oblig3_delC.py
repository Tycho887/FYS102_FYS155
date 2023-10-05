# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 18:09:42 2023

@author: Michael
"""

Specific_heat_vann = 4186
Specific_heat_bly  = 129
Latent_heat = 2260*1000

T1_bly = 327.3
T1_vann = 75.00
M_bly = 1.250
M_vann = 0.5000


E_vann = M_vann*Specific_heat_vann*(100-T1_vann)
E_bly  = M_bly*Specific_heat_bly*(100-T1_bly)

Mass_loss = (E_vann+E_bly)/Latent_heat

print(f"Massen til systemet er {M_vann+M_bly-Mass_loss:.3f} kg")
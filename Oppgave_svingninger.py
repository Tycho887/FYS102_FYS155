import numpy as np
import matplotlib.pyplot as plt


B = 1
m = 100
X0 = 0.2
V0 = 0
w0 = 0.1

# sjekker betingelsen for dempet svingning

assert B < 2*m*w0, "Betingelsen for dempet svingning er ikke oppfylt"

# Finner de spesifikke variablene for den harmoniske svingningen

w_ = np.sqrt(w0**2-(B/(2*m))**2)
A  = np.sqrt(X0**2+(V0/w_+B*X0/(2*m*w_)**2))
theta = np.arctan(w_*X0/(V0+B*X0/(2*m)))

def x(t):
    return A*np.exp(-B*t/(2*m))*np.sin(w_*t+theta)

def Amplitude(t):
    return A*np.exp(-B*t/(2*m))

time = np.linspace(0,1000,1000)


# make a cyan plot of x(t)
plt.plot(time,Amplitude(time),color="red")
#plt.plot(time,x(time),color="cyan")
plt.title("Amplitude plot, logaritmisk")
plt.yscale("log")
plt.xlabel("time [seconds]")
plt.ylabel("x(t) [radians]")
# add a grid
plt.grid()
plt.savefig("Underdempet_svingning03.png",dpi=1024)
plt.show()

import termoPy as TP
import matplotlib.pyplot as plt
import numpy as np

atm = 101300; L = 0.001

initial_volume = TP.IdealGas(n=3,P1=1*atm,T1=273).V1
final_volume = initial_volume*0.50

isoterm_prosess = TP.Isothermal(n=3,T=273,P1=1*atm)
isoterm_prosess.generate_graph_from_dV(V2=final_volume)

adiabatisk_prosess = TP.Adiabatic(n=3,gamma=7/5,T1=273,P1=1*atm)
adiabatisk_prosess.generate_graph_from_dV(V2=final_volume)

def plot_processes(A,B):

    plot_title = A.title + " og " + B.title

    plt.title(plot_title)
    plt.xlabel("Volum [m^3]")
    plt.ylabel("Trykk [Pa]")
    plt.plot(A.volume,A.pressure,label=A.title)
    plt.plot(B.volume,B.pressure,label=B.title)
    plt.legend()
    plt.grid()
    plt.show()

plot_processes(isoterm_prosess,adiabatisk_prosess)

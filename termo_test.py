import termoPy as TP
import matplotlib.pyplot as plt

atm = 1.01325e5
C = lambda c: c + 273.15

# system_states = [{"P":1*atm,    "V":1.0,    "T":300,    "n":1},
#                  {"P":None,     "V":3.0,    "T":None,   "n":1},
#                  {"P":None,     "V":None,   "T":600,    "n":1},
#                  {"P":2*atm,    "V":None,   "T":None,   "n":1}]

# cycle = TP.Carnot(system_states)

# #print(cycle.system_types)
# print(cycle.processes)


P1 = 1.0*atm
T1 = C(0)
n = 1
V1 = n*TP.R*T1/P1

T_high = C(600)
T_low = C(300)



def print_values(V,P,T):
    print(f"V = {V:.3f} m^3   T = {T:.3f} C    P = {P/1e3:.3f} kPa")

isothermal_expansion = TP.Isothermal(n=n,P1=P1,T1=T_high,diatomic=True)
isothermal_expansion.generate_data_from_dV(isothermal_expansion.V1*5)

V2 = isothermal_expansion.volume[-1]
P2 = isothermal_expansion.pressure[-1]
T2 = isothermal_expansion.temperature[-1]

adiabatic_expansion = TP.Adiabatic(n=n,P1=P2,T1=T2,diatomic=True)
adiabatic_expansion.generate_data_from_dT(T_low)

V3 = adiabatic_expansion.volume[-1]
P3 = adiabatic_expansion.pressure[-1]
T3 = adiabatic_expansion.temperature[-1]

isothermal_compression = TP.Isothermal(n=n,P1=P3,T1=T3,diatomic=True)
isothermal_compression.generate_data_from_dP(P1*0.5)

V4 = isothermal_compression.volume[-1]
P4 = isothermal_compression.pressure[-1]
T4 = isothermal_compression.temperature[-1]

adiabatic_compression = TP.Adiabatic(n=n,P1=P4,T1=T4,diatomic=True)
adiabatic_compression.generate_data_from_dT(T_high)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(isothermal_expansion.pressure, isothermal_expansion.volume, isothermal_expansion.temperature, label="Isothermal Expansion")
ax.plot(adiabatic_expansion.pressure, adiabatic_expansion.volume, adiabatic_expansion.temperature, label="Adiabatic Expansion")
ax.plot(isothermal_compression.pressure, isothermal_compression.volume, isothermal_compression.temperature, label="Isothermal Compression")
ax.plot(adiabatic_compression.pressure, adiabatic_compression.volume, adiabatic_compression.temperature, label="Adiabatic Compression")
ax.legend()
ax.set_xlabel("Pressure (Pa)")
ax.set_ylabel("Volume (m^3)")
ax.set_zlabel("Temperature (K)")

print(adiabatic_expansion.consistency)
plt.show()

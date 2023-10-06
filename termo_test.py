import termoPy as TP

atm = 1.01325e5

isobaric = TP.Isobaric(n=1,P1=1.5*atm,T1=300,diatomic=True)
# isobaric.generate_data_from_dT(600,show=True)
# print(isobaric.entropy)

system_states = [{"P":1*atm,    "V":1.0,    "T":300,    "n":1},
                 {"P":None,     "V":3.0,    "T":None,   "n":1},
                 {"P":None,     "V":None,   "T":600,    "n":1},
                 {"P":2*atm,    "V":None,   "T":None,   "n":1}]

cycle = TP.Carnot(system_states)

#print(cycle.system_types)
print(cycle.processes)
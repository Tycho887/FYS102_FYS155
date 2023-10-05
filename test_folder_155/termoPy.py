import numpy as np
import matplotlib.pyplot as plt

atm = 101300; L = 0.001; R = 8.314; k = 10000

class IdealGas:
    def __init__(self,n=None,P1=None,V1=None,T1=None):
        """
        n: number of moles
        gamma: ratio of specific heats
        P1: initial pressure (Pa)
        V1: initial volume (M^3)
        T1: initial temperature (K)
        """
        self.n = n
        self.P1 = P1
        self.V1 = V1
        self.T1 = T1
        


        self.volume = np.array([])
        self.pressure = np.array([])
        self.temperature = np.array([])
        self.internal_energy = np.array([])
        self.work_done_by = 0
        self.heat_absorbed = 0
        self.title = ""

        self._find_missing()
    
    def P(self,V,T):
        return self.n*R*T/V
    def V(self,P,T):
        return self.n*R*T/P
    def T(self,P,V):
        return P*V/(self.n*R)
    def get_n(self,P,V,T):
        return P*V/(R*T)
    
    def _plot_PV(self):
        plt.plot(self.volume,self.pressure)
        plt.xlabel("Volum [m^3]")
        plt.ylabel("Trykk [Pa]")
        plt.grid()
        plt.title(self.title)
        plt.show()

    def _find_missing(self):
        assert self.P1 != None or self.V1 != None or self.T1 != None or self.n != None, "En av P1,V1 eller T1 må være definert"
        if self.P1 == None:
            self.P1 = self.P(self.V1,self.T1)
        elif self.V1 == None:
            self.V1 = self.V(self.P1,self.T1)
        elif self.T1 == None:
            self.T1 = self.T(self.P1,self.V1)
        elif self.n == None:
            self.n = self.get_n(self.P1,self.V1,self.T1)

class Isothermal(IdealGas):
    def __init__(self,T,n,V1=None,P1=None):
        """
        n: number of moles
        T: temperature (K)
        """
        super().__init__(n,P1=P1,V1=V1,T1=T)
        self.T1 = T
        self.V1 = V1
        self.P1 = P1
        self.title = "Isotermisk prosess"
        self._find_missing()

    def work_by(self):
        self.work_done_by = self.n*R*self.T1*np.log(self.volume[-1]/self.volume[0])
        self.heat_absorbed = self.work_done_by
        return self.work_done_by

    def generate_data_from_dV(self,V2,show=False,steps=k):
        self.volume = np.linspace(self.V1,V2,steps)
        self.pressure = self.P(self.volume,self.T1)
        if show:
            self._plot_PV()
        return self.volume,self.pressure

    def generate_data_from_dP(self,P2,show=False,steps=k):
        self.pressure = np.linspace(self.P1,P2,steps)
        self.volume = self.V(self.pressure,self.T1)
        if show:
            self._plot_PV()
        return self.volume,self.pressure


class Adiabatic(IdealGas):
    def __init__(self,n,gamma,P1=None,V1=None,T1=None):
        super().__init__(n,P1=P1,V1=V1,T1=T1)
        self.P1 = P1; self.V1 = V1; self.T1 = T1
        self.gamma = gamma
        self.heat_absorbed = 0
        self.title = "Adiabatisk prosess"
        self._find_missing()

    
    def P2_from_V2(self,V2):
        assert self.P1 != None and self.V1 != None, "P1,V1 må være definert"
        return self.P1*(self.V1/V2)**self.gamma
    
    def T2_from_V2(self,V2):
        assert self.V1 != None and self.T1 != None, "V1,T1 må være definert"
        return self.T1*(self.V1/V2)**(self.gamma-1)
    
    def V2_from_P2(self,P2):
        assert self.P1 != None and self.V1 != None, "P1,V1 må være definert"
        return self.V1*(self.P1/P2)**(1/self.gamma)
    
    def T2_from_P2(self,P2):
        assert self.P1 != None and self.T1 != None, "P1,T1 må være definert"
        return self.T1*(self.P1/P2)**((self.gamma-1)/self.gamma)
    
    def work_by(self):
        self.work_done_by = self.n*R*self.T1*np.log(self.volume[-1]/self.volume[0])
        return self.work_done_by
    
    def generate_data_from_dV(self,V2,show=False, steps = k):
        self.volume = np.linspace(self.V1,V2,steps)
        self.pressure = self.P2_from_V2(self.volume)
        self.temperature = self.T2_from_V2(self.volume)
        if show:
            self._plot_PV()
        return self.volume,self.pressure

    def generate_data_from_dP(self,P2,show=False, steps = k):
        self.pressure = np.linspace(self.P1,P2,steps)
        self.volume = self.V2_from_P2(self.pressure)
        self.temperature = self.T2_from_P2(self.pressure)
        if show:
            self._plot_PV()
        return self.volume,self.pressure


# Ny test
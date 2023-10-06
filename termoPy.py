import numpy as np
import matplotlib.pyplot as plt

atm = 101300; L = 0.001; R = 8.314; k = 10000

class IdealGas:
    def __init__(self,n=None,P1=None,V1=None,T1=None,monatomic=False,diatomic=False,Cv=None):
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
        if monatomic:
            self.Cv = 3/2
        elif diatomic:
            self.Cv = 5/2
        else:
            self.Cv = None
        
        self.volume = None
        self.pressure = None
        self.temperature = None
        self.internal_energy = None
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
    
    def entropy_change(self,Cv=None):
        assert self.volume != None and self.pressure != None and self.temperature != None, "Volum, trykk og temperatur må være definert"
        if Cv != None:
            self.Cv = Cv
        else:
            assert self.Cv != None, "Cv må være definert"
        return self.n*R*np.log(self.V2/self.V1)+self.Cv*np.log(self.T2/self.T1)
    
    def _plot_PV(self):
        plt.plot(self.volume,self.pressure)
        plt.xlabel("Volum [m^3]")
        plt.ylabel("Trykk [Pa]")
        plt.grid()
        plt.title(self.title)
        plt.savefig(f"data/{self.title}_{self.n}.png")
        plt.show()

    def _find_missing(self):
        assert (self.P1 != None) + (self.V1 != None) + (self.T1 != None) + (self.n != None) > 2, "Tre av P1,V1,T1 eller n må være definert"
        if self.P1 == None:
            self.P1 = self.P(self.V1,self.T1)
        elif self.V1 == None:
            self.V1 = self.V(self.P1,self.T1)
        elif self.T1 == None:
            self.T1 = self.T(self.P1,self.V1)
        elif self.n == None:
            self.n = self.get_n(self.P1,self.V1,self.T1)

class Isothermal(IdealGas):
    def __init__(self,n=None,T1=None,V1=None,P1=None):
        """
        n: number of moles
        T: temperature (K)
        """
        super().__init__(n,P1=P1,V1=V1,T1=T1)
        self.title = "Isotermisk prosess"
        self._find_missing()

    def calculate_work_done_by(self):
        self.work_done_by = self.n*R*self.T1*np.log(self.volume[-1]/self.volume[0])
        self.heat_absorbed = self.work_done_by
        return self.work_done_by

    def generate_data_from_dV(self,V2,show=False,steps=k):
        self.volume = np.linspace(self.V1,V2,steps)
        self.pressure = self.P(self.volume,self.T1)
        self.temperature = self.T1*np.ones(len(self.volume))
        if show:
            self._plot_PV()
        return self.volume,self.pressure

    def generate_data_from_dP(self,P2,show=False,steps=k):
        self.pressure = np.linspace(self.P1,P2,steps)
        self.volume = self.V(self.pressure,self.T1)
        self.temperature = self.T1*np.ones(len(self.pressure))
        if show:
            self._plot_PV()
        return self.volume,self.pressure

class Isobaric(IdealGas):
    def __init__(self,n=None,P1=None,T1=None,V1=None,monatomic=False,diatomic=False):
        super().__init__(n,P1=P1,V1=V1,T1=T1,monatomic=monatomic,diatomic=diatomic)
        self.title = "Isobar prosess"
        self._find_missing()
    
    def calculate_work_done_by(self):
        self.work_done_by = self.P1*(self.volume[-1]-self.volume[0])
        return self.work_done_by
    
    def calculate_heat_absorbed(self):
        self.heat_absorbed = self.n*(self.Cv+1)*R*(self.temperature[-1]-self.temperature[0])
        return self.heat_absorbed
    
    def generate_data_from_dV(self,V2,show=False, steps = k):
        self.volume = np.linspace(self.V1,V2,steps)
        self.temperature = self.T(self.P1,self.volume)
        self.pressure = self.P1*np.ones(len(self.volume))
        if show:
            self._plot_PV()
        return self.volume,self.pressure
    
    def generate_data_from_dT(self,T2,show=False, steps = k):
        self.temperature = np.linspace(self.T1,T2,steps)
        self.volume = self.V(self.P1,self.temperature)
        self.pressure = self.P1*np.ones(len(self.temperature))
        if show:
            self._plot_PV()
        return self.volume,self.pressure

class Isochoric(IdealGas):
    def __init__(self,n=None,V1=None,T1=None,P1=None):
        super().__init__(n,P1=P1,V1=V1,T1=T1)
        self.title = "Isokor prosess"
        self._find_missing()
    
    def calculate_work_done_by(self):
        self.work_done_by = 0
        return self.work_done_by
    
    def calculate_heat_absorbed(self):
        self.heat_absorbed = self.n*self.Cv*R*(self.temperature[-1]-self.temperature[0])
        return self.heat_absorbed
    
    def generate_data_from_dT(self,T2,show=False, steps = k):
        self.temperature = np.linspace(self.T1,T2,steps)
        self.pressure = self.P(self.V1,self.temperature)
        self.volume = self.V1*np.ones(len(self.temperature))
        if show:
            self._plot_PV()
        return self.volume,self.pressure
    
    def generate_data_from_dP(self,P2,show=False, steps = k):
        self.pressure = np.linspace(self.P1,P2,steps)
        self.temperature = self.T(self.pressure,self.V1)
        self.volume = self.V1*np.ones(len(self.pressure))
        if show:
            self._plot_PV()
        return self.volume,self.pressure
    
class Adiabatic(IdealGas):
    def __init__(self,gamma,n=None,P1=None,V1=None,T1=None,monatomic=False,diatomic=False):
        super().__init__(n,P1=P1,V1=V1,T1=T1,monatomic=monatomic,diatomic=diatomic)
        self.gamma = gamma
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
    
    def P2_from_T2(self,T2):
        assert self.T1 != None and self.P1 != None, "T1,P1 må være definert"
        return self.P1*(self.T1/T2)**(self.gamma/(self.gamma-1))
    
    def V2_from_T2(self,T2):
        assert self.T1 != None and self.V1 != None, "T1,V1 må være definert"
        return self.V1*(self.T1/T2)**(1/(self.gamma-1))
    
    def calculate_heat_absorbed(self):
        self.heat_absorbed = 0
        return self.heat_absorbed
    
    def calculate_work_done_by(self):
        self.work_done_by = self.Cv * self.n * R * (self.temperature[-1]-self.temperature[0])
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
    
    def generate_data_from_dT(self,T2,show=False, steps = k):
        self.temperature = np.linspace(self.T1,T2,steps)
        self.volume = self.V2_from_T2(self.temperature)
        self.pressure = self.P2_from_T2(self.temperature)
        if show:
            self._plot_PV()
        return self.volume,self.pressure
    
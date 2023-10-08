import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
# we need to know what time it is
import datetime

atm = 101300; L = 0.001; R = 8.314; k = 1.38064852e-23 # Boltzmanns konstant
K = 10000

class IdealGas:
    def __init__(self,n=None,P1=None,V1=None,T1=None,monatomic=False,diatomic=False,Cv=None,specific_heat=None,molar_mass=None,atomic_mass=None,diameter=1e-10,dummy=False):
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

        self.volume = None
        self.pressure = None
        self.temperature = None

        if monatomic:
            self.Cv = 3/2
            self.gamma = 5/3
        elif diatomic:
            self.Cv = 5/2
            self.gamma = 7/5
        elif Cv != None:
            self.Cv = Cv
            self.gamma = Cv+1
        else:
            self.Cv = None
        
        self.specific_heat = specific_heat
        self.molar_mass = molar_mass
        self.atomic_mass = atomic_mass

        self.internal_energy = None
        self.entropy = None


        self.work_done_by = 0
        self.heat_absorbed = 0
        
        self.diameter = diameter # 1 angstrom by default
        self.nv = None # number of particles per volume

        self.rms_speed = None
        self.mean_free_path = None
        self.mean_free_time = None
        
        self.title = ""

        if not dummy:
            self._find_missing()

        self.consistency = None

    def _generate_extra_data(self,show):
        self.internal_energy = self.Cv*self.n*R*self.temperature
        self.entropy = self.n*R*np.log(self.volume)+self.Cv*np.log(self.temperature)
        self.consistency = self.pressure*self.volume-(self.n*R*self.temperature)
        
        if self.molar_mass != None:
            self.rms_speed = np.sqrt(3*self.temperature*R/self.molar_mass)
            self.nv = self.n*6.022e23/self.volume
            self.atomic_mass = self.molar_mass/6.022e23
            self.mean_free_path = 1/(np.sqrt(2)*np.pi*self.diameter**2*self.nv)
            self.mean_free_time = self.mean_free_path/self.rms_speed
        
        if show: self.plot_PV()

    def P(self,V,T):
        return self.n*R*T/V
    
    def V(self,P,T):
        return self.n*R*T/P
    
    def T(self,P,V):
        return P*V/(self.n*R)
    
    def get_n(self,P,V,T):
        return P*V/(R*T)
    
    def generate_data_from_dV(self,V2,show=False,steps=K):
        self.volume      = self.V1*np.ones(steps)
        self.pressure    = self.P1*np.ones(steps)
        self.temperature = self.T1*np.ones(steps)
        self._generate_extra_data(show)
        return self.volume,self.pressure
    
    def generate_data_from_dP(self,P2,show=False,steps=K):
        return self.generate_data_from_dV(self.V1,show=show,steps=steps)
    
    def generate_data_from_dT(self,T2,show=False,steps=K):  
        return self.generate_data_from_dV(self.V1,show=show,steps=steps)
    
    def maxwell_boltzmann_speed_distribution(self,T=None):
        if T == None:
            T = np.max(self.temperature)
        v_max = np.sqrt(2*k*T/self.atomic_mass)
        standard_deviation = np.sqrt(k*T/self.atomic_mass)
        v = np.linspace(v_max-3*standard_deviation,v_max+3*standard_deviation,1000)
        return 4/np.sqrt(np.pi)*(self.atomic_mass/(k*T))**(3/2)*v**2*np.exp(-self.atomic_mass*v**2/(2*k*T))
    
    def calculate_work_done_by(self):
        self.work_done_by = self.P1*(self.volume-self.V1)
        return self.work_done_by
    
    def calculate_heat_absorbed(self):
        self.heat_absorbed = self.n*(self.Cv+1)*R*(self.temperature-self.T1)
        return self.heat_absorbed

    def _show_picture(self,save=False,name=None,xlabel=None,ylabel=None,type=""):
        plt.grid()
        plt.legend()
        plt.title(f"{self.title}")
        if xlabel != None and ylabel != None:
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
        # we want to save the fig with the date as a name
        
        if save and name==None: 
            now = datetime.datetime.now()
            plt.savefig(f"{self.title}_{type}_{now.strftime('%Y_%m_%d')}.png",dpi=1024)
        elif save: plt.savefig(f"{name}.png",dpi=1024)
        plt.show()

    def plot_PV(self,save=False,name=None):
        plt.plot(self.volume,self.pressure,label=self.title)
        plt.scatter(self.volume[0],self.pressure[0],label="Startpunkt")
        plt.scatter(self.volume[-1],self.pressure[-1],label="Sluttpunkt")
        self._show_picture(save,name,"Volum [m^3]", "Trykk [Pa]","PV")

    def plot_PT(self,save=False,name=None):
        plt.plot(self.temperature,self.pressure,label=self.title)
        plt.scatter(self.temperature[0],self.pressure[0],label="Startpunkt")
        plt.scatter(self.temperature[-1],self.pressure[-1],label="Sluttpunkt")
        self._show_picture(save,name,"Temperatur [K]", "Trykk [Pa]", "PT")

    def plot_VT(self,save=False,name=None):
        plt.plot(self.temperature,self.volume,label=self.title)
        plt.scatter(self.temperature[0],self.volume[0],label="Startpunkt")
        plt.scatter(self.temperature[-1],self.volume[-1],label="Sluttpunkt")
        self._show_picture(save,name,"Temperatur [K]", "Volum [m^3]", "VT")

    def plot_ST(self,save=False,name=None):
        plt.plot(self.temperature,self.entropy,label=self.title)
        plt.scatter(self.temperature[0],self.entropy[0],label="Startpunkt")
        plt.scatter(self.temperature[-1],self.entropy[-1],label="Sluttpunkt")
        self._show_picture(save,name,"Temperatur [K]", "Entropi [J/K]", "ST")

    def plot_PVT(self,save=False,name=None):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(self.volume,self.pressure,self.temperature,label=self.title)
        ax.scatter(self.volume[0],self.pressure[0],self.temperature[0],label="Startpunkt")
        ax.scatter(self.volume[-1],self.pressure[-1],self.temperature[-1],label="Sluttpunkt")
        ax.set_xlabel("Volum [m^3]")
        ax.set_ylabel("Trykk [Pa]")
        ax.set_zlabel("Temperatur [K]")
        self.show_picture(save,name,type="PVT")

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
    
    def __str__(self): # the __str__ method is used when printing the object
        return f"n: {self.n}\nP1: {self.P1}\nV1: {self.V1}\nT1: {self.T1}\nCv: {self.Cv}\ngamma: {self.gamma}\n"
class Isothermal(IdealGas):
    def __init__(self,n=None,T1=None,V1=None,P1=None,monatomic=False,diatomic=False):
        """
        n: number of moles
        T: temperature (K)
        """
        super().__init__(n,P1=P1,V1=V1,T1=T1,monatomic=monatomic,diatomic=diatomic)
        self.title = "Isotermisk prosess"
        self._find_missing()

    def calculate_heat_absorbed(self):
        self.heat_absorbed = self.calculate_work_done_by()
        return self.heat_absorbed

    def generate_data_from_dV(self,V2,show=False,steps=K):
        self.volume = np.linspace(self.V1,V2,steps)
        self.pressure = self.P(self.volume,self.T1)
        self.temperature = self.T1*np.ones(len(self.volume))
        self._generate_extra_data(show)
        return self.volume,self.pressure

    def generate_data_from_dP(self,P2,show=False,steps=K):
        self.pressure = np.linspace(self.P1,P2,steps)
        self.volume = self.V(self.pressure,self.T1)
        self.temperature = self.T1*np.ones(len(self.pressure))
        self._generate_extra_data(show)
        return self.volume,self.pressure
class Isobaric(IdealGas):
    def __init__(self,n=None,P1=None,T1=None,V1=None,monatomic=False,diatomic=False):
        super().__init__(n,P1=P1,V1=V1,T1=T1,monatomic=monatomic,diatomic=diatomic)
        self.title = "Isobar prosess"
        self._find_missing()
    
    def generate_data_from_dV(self,V2,show=False, steps = K):
        self.volume = np.linspace(self.V1,V2,steps)
        self.temperature = self.T(self.P1,self.volume)
        self.pressure = self.P1*np.ones(len(self.volume))
        self._generate_extra_data(show)
        return self.volume,self.pressure
    
    def generate_data_from_dT(self,T2,show=False, steps = K):
        self.temperature = np.linspace(self.T1,T2,steps)
        self.volume = self.V(self.P1,self.temperature)
        self.pressure = self.P1*np.ones(len(self.temperature))
        self._generate_extra_data(show)
        return self.volume,self.pressure
class Isochoric(IdealGas):
    def __init__(self,n=None,V1=None,T1=None,P1=None,monatomic=False,diatomic=False):
        super().__init__(n,P1=P1,V1=V1,T1=T1,monatomic=monatomic,diatomic=diatomic)
        self.title = "Isokor prosess"
        self._find_missing()
    
    def calculate_work_done_by(self):
        self.work_done_by = 0
        return self.work_done_by
    
    def generate_data_from_dT(self,T2,show=False, steps = K):
        self.temperature = np.linspace(self.T1,T2,steps)
        self.pressure = self.P(self.V1,self.temperature)
        self.volume = self.V1*np.ones(len(self.temperature))
        self._generate_extra_data(show)
        return self.volume,self.pressure
    
    def generate_data_from_dP(self,P2,show=False, steps = K):
        self.pressure = np.linspace(self.P1,P2,steps)
        self.temperature = self.T(self.pressure,self.V1)
        self.volume = self.V1*np.ones(len(self.pressure))
        self._generate_extra_data(show)
        return self.volume,self.pressure
class Adiabatic(IdealGas):
    def __init__(self,gamma=None,n=None,P1=None,V1=None,T1=None,monatomic=False,diatomic=False):
        super().__init__(n,P1=P1,V1=V1,T1=T1,monatomic=monatomic,diatomic=diatomic)
        if gamma != None:
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
        return self.T1*(P2/self.P1)**((self.gamma-1)/self.gamma)
    
    def P2_from_T2(self,T2):
        assert self.T1 != None and self.P1 != None, "T1,P1 må være definert"
        return self.P1*(T2/self.T1)**(self.gamma/(self.gamma-1))
    
    def V2_from_T2(self,T2):
        assert self.T1 != None and self.V1 != None, "T1,V1 må være definert"
        return self.V1*(self.T1/T2)**(1/(self.gamma-1))
    
    def calculate_heat_absorbed(self):
        self.heat_absorbed = np.zeros(len(self.volume))
        return self.heat_absorbed
    
    def generate_data_from_dV(self,V2,show=False, steps = K):
        self.volume = np.linspace(self.V1,V2,steps)
        self.pressure = self.P2_from_V2(self.volume)
        self.temperature = self.T2_from_V2(self.volume)
        self._generate_extra_data(show)
        return self.volume,self.pressure

    def generate_data_from_dP(self,P2,show=False, steps = K):
        self.pressure = np.linspace(self.P1,P2,steps)
        self.volume = self.V2_from_P2(self.pressure)
        self.temperature = self.T2_from_P2(self.pressure)
        self._generate_extra_data(show)
        return self.volume,self.pressure
    
    def generate_data_from_dT(self,T2,show=False, steps = K):
        self.temperature = np.linspace(self.T1,T2,steps)
        self.volume = self.V2_from_T2(self.temperature)
        self.pressure = self.P2_from_T2(self.temperature)
        self._generate_extra_data(show)
        return self.volume,self.pressure
    
class Cycle:
    def __init__(self,process_types,system_states,monatomic=False,diatomic=False,Cv=None,specific_heat=None,molar_mass=None,diameter=1e-10,gamma=None):
        
        self.system_states = system_states
        self.process_types = process_types
        self.processes = []

        self.work_done_by = 0
        self.heat_absorbed = 0

        self.monatomic = monatomic
        self.diatomic = diatomic

        self.Cv = Cv
        self.gamma = gamma
        self.specific_heat = specific_heat
        self.molar_mass = molar_mass
        self.diameter = diameter

        assert len(self.system_states) == len(self.process_types), "system_states og process_types må ha samme lengde"

        self.title = ""

    def find_missing_variables(self,process,system_state):
        
        testingmode = False

        if system_state["P"] == None:
            try: process.generate_data_from_dV(system_state["V"],show=testingmode)
            except: process.generate_data_from_dT(system_state["T"],show=testingmode)
        elif system_state["V"] == None:
            try: process.generate_data_from_dP(system_state["P"],show=testingmode)
            except: process.generate_data_from_dT(system_state["T"],show=testingmode)
        elif system_state["T"] == None:
            try: process.generate_data_from_dV(system_state["V"],show=testingmode)
            except: process.generate_data_from_dP(system_state["P"],show=testingmode)
        elif system_state["n"] == None:
            try: process.get_n(system_state["P"],system_state["V"],system_state["T"],show=testingmode)
            except: assert False, "Det er ikke definert nok variabler i system_state"
        else:
            print(system_state)
        
        return process.volume[-1],process.pressure[-1],process.temperature[-1],process.n

    def define_processes(self,process_type,system_state):
        
        if process_type == "Isothermal":
            process = Isothermal(n=system_state["n"],T1=system_state["T"],V1=system_state["V"],P1=system_state["P"],monatomic=self.monatomic,diatomic=self.diatomic)
        elif process_type == "Isobaric":
            process = Isobaric(n=system_state["n"],T1=system_state["T"],V1=system_state["V"],P1=system_state["P"],monatomic=self.monatomic,diatomic=self.diatomic)
        elif process_type == "Isochoric":
            process = Isochoric(n=system_state["n"],T1=system_state["T"],V1=system_state["V"],P1=system_state["P"],monatomic=self.monatomic,diatomic=self.diatomic)
        elif process_type == "Adiabatic":
            process = Adiabatic(n=system_state["n"],T1=system_state["T"],V1=system_state["V"],P1=system_state["P"],monatomic=self.monatomic,diatomic=self.diatomic,gamma=self.Cv+1)
        else:
            assert False, "Prosessen er ikke definert"
        
        return process
            
    def generate_processes(self):
        for process_name,system_state in zip(self.process_types,self.system_states):
            try:
                # find pressure, volume and temperature at the end of the last process, and use these as the initial values for the next process
                system_state["P"] = self.processes[-1].pressure[-1]
                system_state["V"] = self.processes[-1].volume[-1]
                system_state["T"] = self.processes[-1].temperature[-1]
                system_state["n"] = self.processes[-1].n
                process = self.define_processes(process_name,system_state)
                self.find_missing_variables(process,system_state)
                self.processes.append(process)
            except:
                # if this is the first process, use the initial values from the system states
                process = self.define_processes(process_name,system_state)
                self.find_missing_variables(process,system_state)
                self.processes.append(process)


    def _plot_PV(self):
        for process in self.processes:
            plt.plot(process.volume,process.pressure,label=process.title)
        plt.legend()
        plt.xlabel("Volum [m^3]")
        plt.ylabel("Trykk [Pa]")
        plt.grid()
        plt.title(self.title)
        plt.show()
    
    def _calculate_work_done_by(self):
        for process in self.processes:
            self.work_done_by += process.calculate_work_done_by()

    def _calculate_heat_absorbed(self):

        for process in self.processes:
            # we only want to add heat absorbed, meaning positive values
            heat_absorbed = process.calculate_heat_absorbed()
            if heat_absorbed > 0:
                self.heat_absorbed += heat_absorbed
    
    def _calculate_efficiency(self):
        self.efficiency = self.work_done_by/self.heat_absorbed
class Carnot(Cycle):
    def __init__(self,system_states,monatomic=False,diatomic=False,Cv=None,specific_heat=None,molar_mass=None,diameter=1e-10):
        super().__init__(system_states=system_states,
                         monatomic=monatomic,
                         diatomic=diatomic,
                         Cv=Cv,
                         specific_heat=specific_heat,
                         molar_mass=molar_mass,
                         diameter=diameter,
                         process_types=["Isothermal","Adiabatic","Isothermal","Adiabatic"])
        self.title = "Carnot syklus"
        self.generate_processes()
class Otto(Cycle): 
    pass


def __chceck_consisency(process,type):
#        if process.title == "Adiabatisk prosess":
        assert (np.mean(process.consistency)) < 1e-10, "Consistency is not zero"
        print(f"{type} method is consistent. mean consistency: {np.mean(process.consistency):.2e}")

def __test_processes(processes):
    for process in processes:
        print("Testing",process.title)
        try:
            process.generate_data_from_dV(process.V1*10)
            __chceck_consisency(process,"+dV")
            process.generate_data_from_dV(process.volume[-1]*0.1)
            __chceck_consisency(process,"-dV")
        except: assert False, f"The from_dV method does not work as intended on {process.title}"
        try:
            process.generate_data_from_dP(process.P1*10)
            __chceck_consisency(process,"+dP")
            process.generate_data_from_dP(process.pressure[-1]*0.1)
            __chceck_consisency(process,"-dP")
        except: assert False, f"The from_dP method does not work as intended on {process.title}"
        try: 
            process.generate_data_from_dT(process.T1*10)
            __chceck_consisency(process,"+dT")
            process.generate_data_from_dT(process.temperature[-1]*0.1)
            __chceck_consisency(process,"-dT")
        except: assert False, f"The from_dT method does not work as intended on {process.title}"
            
if __name__ == "__main__":
    IS = {"P":1*atm,"V":0.001,"T":300,"n":1}
    processes = [Isothermal(n=IS["n"],V1=IS["V"],P1=IS["P"],monatomic=True),
                 Isochoric(n=IS["n"],T1=IS["T"],P1=IS["P"],monatomic=True),
                 Isobaric(n=IS["n"],V1=IS["V"],T1=IS["T"],monatomic=True),
                 Adiabatic(n=IS["n"],V1=IS["V"],T1=IS["T"],monatomic=True)]
    __test_processes(processes)
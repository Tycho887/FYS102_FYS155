import numpy as np
import matplotlib.pyplot as plt

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
        elif diatomic:
            self.Cv = 5/2
        elif Cv != None:
            self.Cv = Cv
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
    
    def P(self,V,T):
        return self.n*R*T/V
    def V(self,P,T):
        return self.n*R*T/P
    def T(self,P,V):
        return P*V/(self.n*R)
    def get_n(self,P,V,T):
        return P*V/(R*T)
    
    def _generate_internal_energy(self):
        self.internal_energy = self.Cv*self.n*R*self.temperature

    def _generate_entropy(self):
        self.entropy = self.n*R*np.log(self.volume)+self.Cv*np.log(self.temperature)

    def _generate_rms_speed(self):
        self.rms_speed = np.sqrt(3*self.temperature*R/self.molar_mass)

    def generate_particle_density(self):
        self.nv = self.n*6.022e23/self.volume

    def _generate_atomic_mass(self):
        assert self.molar_mass != None, "molar_mass må være definert"
        self.atomic_mass = self.molar_mass/6.022e23
    
    def _generate_mean_free_path(self):
        self.mean_free_path = 1/(np.sqrt(2)*np.pi*self.diameter**2*self.nv)

    def _generate_mean_free_time(self):
        assert self.rms_speed != None and self.mean_free_path != None, "rms_speed eller mean_free_path er ikke definert"
        self.mean_free_time = self.mean_free_path/self.rms_speed

    def _generate_extra_data(self,show):
        self._generate_internal_energy()
        self._generate_entropy()
        if self.molar_mass != None:
            self._generate_rms_speed()
            self.generate_particle_density()
            self._generate_atomic_mass()
            self._generate_mean_free_path()
            self._generate_mean_free_time()

        if show: self._plot_PV()

    def maxwell_boltzmann_speed_distribution(self,T=None):
        if T == None:
            T = np.max(self.temperature)
        v_max = np.sqrt(2*k*T/self.atomic_mass)
        standard_deviation = np.sqrt(k*T/self.atomic_mass)
        v = np.linspace(v_max-3*standard_deviation,v_max+3*standard_deviation,1000)
        return 4/np.sqrt(np.pi)*(self.atomic_mass/(k*T))**(3/2)*v**2*np.exp(-self.atomic_mass*v**2/(2*k*T))
    
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
    
    def calculate_work_done_by(self):
        self.work_done_by = self.P1*(self.volume[-1]-self.volume[0])
        return self.work_done_by
    
    def calculate_heat_absorbed(self):
        self.heat_absorbed = self.n*(self.Cv+1)*R*(self.temperature[-1]-self.temperature[0])
        return self.heat_absorbed
    
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
        self.work_done_by = self.Cv * self.n * R * (self.temperature[0]-self.temperature[-1])
        return self.work_done_by
    
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
    def __init__(self,process_types,system_states,monatomic=False,diatomic=False,Cv=None,specific_heat=None,molar_mass=None,diameter=1e-10):
        
        self.system_states = system_states
        self.process_types = process_types
        self.processes = []

        self.work_done_by = 0
        self.heat_absorbed = 0

        self.monatomic = monatomic
        self.diatomic = diatomic

        self.Cv = Cv
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
            process = Isothermal(n=system_state["n"],T1=system_state["T"],V1=system_state["V"],P1=system_state["P"])
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
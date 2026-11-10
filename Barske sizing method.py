import math
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI
#Barske sizing method

class impeller:
    def __init__(self, fluid, mdot, a_in, P_in, t, n, dP, k_d):
        self.fluid = fluid
        self.mdot = mdot
        self.a_in = a_in
        self.P_in = P_in
        self.t = t
        self.n = n
        self.dP = dP
        self.k_d = k_d
        self.g = 9.80665

    def size(self):
        rho = PropsSI('D','T|liquid',self.t,'P',self.P_in,self.fluid)
        H_ini = self.dP/(rho*self.g) #Total differential head required for dP of design
        c_0 = self.mdot/(rho*self.a_in)
        d_0 = math.sqrt((4*self.mdot)/(math.pi*rho*c_0)) #initial pump inlet diameter
        d_1 = 1.1*d_0 #The impeller blade inlet diameter is enlarged by a small amount as some pre-rotation is desirable for increased head.
        d_2 = math.sqrt(((1/(1+self.k_d))*((7200*self.dP)/(rho*(self.n**2)*(math.pi**2))+d_1**2))) #Rearranged mass continuity equation to find impeller tip diameter.
        u_2 = (math.pi*d_2*self.n)/60
        psi = self.g*H_ini/(u_2**2) #Head coefficient
        print('''
Pump inlet diameter: ''', d_1*1000, '''mm'''
              '''
Impeller diamater: ''', d_2*1000, '''mm''')
        s_ax = 1.1*(d_2/100)
        b_1 = (1.1*(1/4))*d_1
        b_2 = (1.1*3*s_ax)
        B = 2*s_ax+b_2
        H = B/2
        if 1.5 < B/H < 2.5:
            print("B/H check passed")
        else:
            print("The calculation is not valid. Adjust factors")
        d_3 = 0.85*u_2
        d_4 = math.sqrt(4*(3.5*((math.pi*(d_3**2))/4))/math.pi)
        return (c_0, rho, d_0, d_1, d_2, s_ax, b_1, b_2, d_3, d_4, u_2)

    #def losses(self):
        #P_RR = 1956*self.rho*self.v**0.2*(self.n/1000)**2.8*(d_2**4.6+4.6*d_1**3.6*b_1)

    def locks_refined(self, c_0, rho, d_0, d_1, d_2, s_ax, b_1, b_2, d_3, d_4, u_2):
        #This model has been shown to underpredict static head rise supposedly due to physical rounding of peripheral leading edges on the impeller which reduce the effective head.
        H_array = []
        Q_array = []
        h_0 =  0.7   #Sci-Hub the recommended values from "Thomas Lock. A forced vortex pump for high speed, high pressure, low flow applications. SAE Transactions, 74:729–737, 1966."
        C_H =  0.85   #Sci-Hub the recommended values from "Thomas Lock. A forced vortex pump for high speed, high pressure, low flow applications. SAE Transactions, 74:729–737, 1966."
        Q = 0
        H = 0
        K_1 = 0.45
        P_vap = CP.PropsSI("P", "T", self.t, "Q", 0, self.fluid)
        zeta_D = 0.194 #Diffuser loss coefficient (might change this idk really what it is)
        Q_opt = math.sqrt(h_0*(u_2)**2)/((24/(math.pi**2))*((zeta_D/(d_3**4))*(1-((d_3/d_4)**2))+1/(d_4**4)-(1/(d_0**4)))) #Optimum volumetric flow-rate capacity to extract fluid through diffuser without turbulent backflow.
        H_vap = (1/(rho*self.g))*(self.P_in-P_vap) #Vapour pressure head
        while H <= H_vap:
             H_T2 = h_0*(u_2**2/self.g)-C_H*K_1*(1-(Q/Q_opt)*(d_1/d_2)*((u_2**2)/self.g)) #Refined total head developed by the impeller.
             c_3 = Q/((math.pi*(d_3**2))/4)
             c_4 = Q/((math.pi*(d_4**2))/4)
             H_NL = zeta_D*((1-(d_3/d_4)**2)**2)*(((c_3)**2)/(2*self.g)) #Diffuser static head losses
             H = H_T2 - H_NL - (c_4**2)/(2*self.g) + (c_0**2)/(2*self.g) #Static head developed across the pump.
             if Q == Q_opt:
                H_opt = H
             H_array.append(H)
             Q_array.append(Q)
             Q+=0.1
        Q_cut = math.sqrt((1/8)*self.g*(math.pi**2)*(H_vap-H))/((zeta_D/(d_3**4))*(1-(d_3/d_4)**2)+1/(d_4**4)-1/(d_3**4)) #Flowrate at which cavitation cut-out occurs. (Cavitation in diffuser throat)
        plt.plot(H_array, Q_array)
        plt.ylabel('H (m)')
        plt.xlabel('Q (L/s)')
        plt.axvline(x = Q_opt, color = 'r', linestyle = '--', label = 'Optimum Q (L/s)')
        plt.axvline(x = Q_cut, color = 'r', linestyle = '--', label = 'Q Cut-out (L/s)')
        plt.legend()
        plt.show()
        return(H, Q_opt, Q_cut)

    #def performance(self):

           
     
       

IPA_pump = impeller(fluid = "Ethanol", mdot = 0.7, a_in = 0.00009369, P_in = 101325, t = 293.15, n = 30000, dP = 3000000, k_d = 0.15)

print(IPA_pump.size())
size = IPA_pump.size()
print(size)
print(IPA_pump.locks_refined(*size))
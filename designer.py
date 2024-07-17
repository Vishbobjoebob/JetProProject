class engineDesigner():
    
    def __init__(self, Ta, pa, M, Pr_c, Pr_f, beta, b, f, f_ib, f_ab):
        # Flight Conditions
        self.Ta = Ta # Ambient Temperature
        self.pa = pa # Ambient Pressure
        self.M = M # Mach Number
        self.Pr_c = Pr_c # Compressor Pressure Ratio
        self.Pr_f = Pr_f # Fan Pressure Ratio
        self.beta = beta # Bypass Ratio
        self.b = b # Bleed Air to Core Air Ratio
        self.f = f # Fuel-Air Mixture Ratio
        self.f_ib = f_ib # Fuel-Air Mixture Ratio for Interturbine
        self.f_ab = f_ab # Fuel-Air Mixture Ratio for Interturbine
  
    def diffuser(self):
        if self.M <= 1:
            self.rd = 1
        elif self.M > 1 and self.M < 5:
            self.rd = 1 - (0.075*((self.M-1)**1.35))
        
        ad_eff_d = 0.94 # Adiabatic Efficiency of Diffuser
        gamma_ratio = 3.5 # cp/R or gamma/(gamma-1)

        gamma_d = gamma_ratio/(gamma_ratio-1)

        self.p01 = self.rd*self.pa*(1+ad_eff_d*((gamma_d-1)/2)*(self.M**2))**(gamma_d/(gamma_d-1))
        self.T01 = self.Ta*(1+(gamma_d-1)/2*self.M**2)
    
    def bypass_fan(self):
        MW = 0.0298 # Molecular weight of gas in fan (kg/mol)
        gamma_ratio = 3.5 # cp/R or gamma/(gamma-1)
        poly_eff_f = 0.92 # Polytropic Efficiency of Fan

        cp_f = gamma_ratio * (8.314/MW) 
        ad_eff_f = (self.Pr_f**(1/gamma_ratio)-1)/(self.Pr_f**((1/gamma_ratio)/poly_eff_f)-1)

        self.p02 = self.Pr_f*self.p01
        self.T02 = self.T01*(1+1/ad_eff_f*(self.Pr_f**(1/gamma_ratio)-1))
        self.w_f = (1+self.beta)*cp_f*(self.T02-self.T01)
    
    def compressor(self):
        MW = 0.0298 # Molecular weight of gas in compressor (kg/mol)
        gamma_ratio = 3.62 # cp/R or gamma/(gamma-1)
        poly_eff_c = 0.91 # Polytropic Efficiency of Compressor

        cp_c = gamma_ratio * (8.314/MW) 
        ad_eff_c = (self.Pr_c**(1/gamma_ratio)-1)/(self.Pr_c**((1/gamma_ratio)/poly_eff_c)-1)

        self.p03 = self.Pr_c*self.p02
        self.T03 = self.T02*(1+1/ad_eff_c*(self.Pr_c**(1/gamma_ratio)-1))
        self.w_c = cp_c*(self.T03-self.T02)
    
    def combustor_main(self):
        self.Pr_b = 0.95
        MW = 0.0298 # Molecular weight of gas in combustor (kg/mol)
        eff_b = 0.99 # Efficiency of Main Combustor

        cp_b = (8.314/MW) * (3.70 + 0.66*((self.T03/1000)**2) - 0.20*((self.T03/1000)**3))
        delta_hr = 43.52 * 10**6 # Heating Rate of Fuel (J/kg)

        print(self.f)
        print(cp_b)
        print(self.T03)
        print(self.b)

        self.p04 = self.Pr_b * self.p03
        self.T04 = (self.f*(eff_b*delta_hr/cp_b)+self.T03*(1-self.b))/(1-self.b+self.f)
    
    def turbine(self): # Powers the compressor
        MW = 0.0298 # Molecular weight of gas in turbine (kg/mol)
        poly_eff_t = 0.94 # Polytropic Efficiency of Turbine
        Tmax_0 = 1500
        bmax = 0.1
        n = 0.6
        Cbl = 500

        cp_t = (8.314/MW) * (3.38 + 0.70*(self.T04/1000)**2 - 0.20*(self.T04/1000)**3)
        self.T051 = self.T04 - self.w_c/(cp_t*(1-self.b+self.f))
        Tr_t = self.T051/self.T04
        ad_eff_t = (Tr_t-1)/(Tr_t**(1/poly_eff_t)-1)
        self.p051 = self.p04 * (1-1/ad_eff_t*(1-Tr_t))**(cp_t/(8.314/MW))
        self.Tmax_t = Tmax_0 + Cbl * (self.b/bmax)**n

    def mixer_turbine(self):
       MW = 0.0298 # Molecular weight of gas in turbine (kg/mol)

       cp_tm_n = (3.43 + 0.72*(self.T051/1000)**2 - 0.21*(self.T051/1000)**3) # Normalized cp for turbine mixer (cp/R)

       self.T051m = self.T051 + self.b*(self.T03-self.T051)/(1+self.f)
       self.p051m = self.p051 * ((self.T051m/self.T051)/((self.T03/self.T051)**(self.b/(1+self.f)))**cp_tm_n)
    
    def combustor_interturbine(self):
        self.Pr_ib = 0.94 # Interturbine Combustor Pressure Ratio
        MW = 0.0298 # Molecular weight of gas in Interturbine Combustor (kg/mol)
        eff_ib = 0.99 # Efficiency of Interturbine Combustor

        cp_ib = (8.314/MW) * (3.82 + 0.32*(self.T051m/1000)**2 - 0.060*(self.T051m/1000)**3)
        delta_hr = 43.52 * 10**6 # Heating Rate of Fuel (J/kg)

        self.p0514 = self.Pr_ib * self.p051m
        self.T0514 = (self.f_ib*(eff_ib*delta_hr/cp_ib)+(1+self.f)*self.T051m)/(1+self.f+self.f_ib)
        self.T_max_ib = 1500

    def turbine_fan(self): # Powers the fan 
        MW = 0.0298 # Molecular weight of gas in turbine (kg/mol)
        poly_eff_ft = 0.94 # Polytropic Efficiency of Turbine

        cp_ft = (8.314/MW) * (3.40 + 0.63*(self.T0514/1000)**2 - 0.18*(self.T0514/1000)**3)
        self.T052 = self.T0514 - self.w_f/cp_ft

        Tr_ft = self.T052/self.T0514
        ad_eff_ft = (Tr_ft-1)/(Tr_ft**(1/poly_eff_ft)-1)

        self.p052 = self.p0514 * (1-1/ad_eff_ft*(1-Tr_ft))**(cp_ft/(8.314/MW))
    
    def combustor_afterburner(self):
        self.Pr_ab = 0.97 # Interturbine Combustor Pressure Ratio

        MW = 0.0298 # Molecular weight of gas in Interturbine Combustor (kg/mol)
        eff_ab = 0.96 # Efficiency of Interturbine Combustor

        cp_ab = (8.314/MW) * (3.50 + 0.72*(self.T052/1000)**2 - 0.21*(self.T052/1000)**3)
        delta_hr = 43.52 * 10**6 # Heating Rate of Fuel (J/kg)

        self.p06 = self.Pr_ab * self.p052
        self.T06 = (self.f_ab*(eff_ab*delta_hr/cp_ab)+(1+self.f+self.f_ib)*self.T052)/(1+self.f+self.f_ib+self.f_ab)
        self.T_max_ab = 2300

    def nozzle_core(self):
        ad_eff_n = 0.96
        cp_n_n = 3.45 + 0.55*(self.T06/1000)**2 - 0.15*(self.T06/1000)**3 # Normalized cp for core nozzle (cp/R)

        self.pe = self.pa
        self.Te = self.T06*(1-ad_eff_n*(1-(self.pe/self.p06)**(1/cp_n_n)))

    def nozzle_fan(self):
        ad_eff_fn = 0.97
        cp_fn_n = 3.5 # Normalized cp for fan nozzle (cp/R)

        self.pef = self.pa
        self.Tef = self.T02*(1-ad_eff_fn*(1-(self.pef/self.p02)**(1/cp_fn_n)))
    
    def mixer_nozzle(self):
        C_nm = 0.029
        Pr_nm = 1-C_nm*self.beta

        self.fuel_sum = 1+self.f+self.f_ib+self.f_ab

        self.T07 = (self.beta*self.T02 + self.T06*self.fuel_sum)/(self.fuel_sum+self.beta)
        gamma_nm = 1.44 - 0.139*(self.T07/1000) + 0.0357*(self.T07/1000)**2 - 0.004*(self.T07/1000)**3
        gamma_ratio = (gamma_nm)/(gamma_nm-1)

        self.p07 = Pr_nm*(self.T07**gamma_ratio)*(self.p06**(self.fuel_sum/(self.fuel_sum+self.beta)))*(self.p02**(self.beta/(self.fuel_sum+self.beta)))/(self.T06**(gamma_ratio*self.fuel_sum/(self.fuel_sum+self.beta)))/(self.T02**(self.beta*gamma_ratio/(self.fuel_sum+self.beta)))

    def nozzle_combined(self):
        ad_eff_nc = 0.96
        cp_nc_n = 3.45 + 0.55*(self.T07/1000)**2 - 0.15*(self.T07/1000)**3 # Normalized cp for core nozzle (cp/R)

        self.pec = self.pa
        self.Tec = self.T07*(1-ad_eff_nc*(1-(self.pec/self.p07)**(1/cp_nc_n)))
    
        



    




        
        
        
    


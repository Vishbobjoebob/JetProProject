classdef engineDesigner
    properties
        % % Flight Conditions %%
        Ta % Ambient Temperature
        pa % Ambient Pressure
        M % Mach Number
        Pr_c % Compressor Pressure Ratio
        Pr_f % Fan Pressure Ratio
        beta % Bypass Ratio
        b % Bleed Air to Core Air Ratio
        f % Fuel-Air Mixture Ratio
        f_ib % Fuel-Air Mixture Ratio for Interturbine
        f_ab % Fuel-Air Mixture Ratio for Afterburner

        % Intermediate variables
        rd
        p01
        T01
        p02
        T02
        w_f
        p03
        T03
        w_c
        p04
        T04
        T051
        p051
        Tmax_t
        T051m
        p051m
        p0514
        T0514
        T_max_ib
        T052
        p052
        T06
        p06
        pe
        Te
        pef
        Tef
        fuel_sum
        T07
        p07
        pec
        Tec
        Pr_b
        Pr_ib
        Pr_ab
        T_max_ab
        w_p
        w_t
        p_f1
        p_f2
        f_max_main
        f_max_ib
        f_max_ab
        ST_cn
        TSFC_cn
        nth_cn
        n0_cn
        np_cn
        ST_sn
        TSFC_sn
        nth_sn
        n0_sn
        np_sn
        ue
        Me
        uef
        Mef
        uec
        Mec

    end
    
    methods
        function obj = engineDesigner(ambient, Pr_c, Pr_f, beta, b, f, f_ib, f_ab)
            obj.Ta = ambient(1);
            obj.pa = ambient(2);
            obj.M = ambient(3);
            obj.Pr_c = Pr_c;
            obj.Pr_f = Pr_f;
            obj.beta = beta;
            obj.b = b;
            obj.f = f;
            obj.f_ib = f_ib;
            obj.f_ab = f_ab;

            obj = obj.diffuser();
            obj = obj.bypass_fan();
            obj = obj.compressor();
            obj = obj.fuel_pump();
            obj = obj.combustor_main();
            obj = obj.turbine();
            obj = obj.mixer_turbine();
            obj = obj.combustor_interturbine();
            obj = obj.turbine_fan();
            obj = obj.combustor_afterburner();
            obj = obj.nozzle_core();
            obj = obj.nozzle_fan();
            obj = obj.mixer_nozzle();
            obj = obj.nozzle_combined();
            obj = obj.get_performance_parameters_cn();
            obj = obj.get_performance_parameters_sn();
        end

        
        function obj = diffuser(obj)
            if obj.M <= 1
                obj.rd = 1;
            elseif obj.M > 1 && obj.M < 5
                obj.rd = 1 - (0.075*((obj.M-1)^1.35));
            end
            
            ad_eff_d = 0.94; % Adiabatic Efficiency of Diffuser
            gamma_ratio = 3.5; % cp/R or gamma/(gamma-1)
            gamma_d = gamma_ratio/(gamma_ratio-1);
            
            obj.p01 = obj.rd * obj.pa * (1 + ad_eff_d * ((gamma_d - 1) / 2) * (obj.M^2))^(gamma_d / (gamma_d - 1));
            obj.T01 = obj.Ta * (1 + (gamma_d - 1) / 2 * obj.M^2);
        end
        
        function obj = bypass_fan(obj)
            MW = 0.0289; % Molecular weight of gas in fan (kg/mol)
            gamma_ratio = 3.5; % cp/R or gamma/(gamma-1)
            poly_eff_f = 0.92; % Polytropic Efficiency of Fan
            
            cp_f = gamma_ratio * (8.314 / MW);
            ad_eff_f = (obj.Pr_f^(1 / gamma_ratio) - 1) / (obj.Pr_f^((1 / gamma_ratio) / poly_eff_f) - 1);
            
            obj.p02 = obj.Pr_f * obj.p01;
            obj.T02 = obj.T01 * (1 + 1 / ad_eff_f * (obj.Pr_f^(1 / gamma_ratio) - 1));
            obj.w_f = (1 + obj.beta) * cp_f * (obj.T02 - obj.T01);
        end
        
        function obj = compressor(obj)
            MW = 0.0289; % Molecular weight of gas in compressor (kg/mol)
            gamma_ratio = 3.62; % cp/R or gamma/(gamma-1)
            poly_eff_c = 0.91; % Polytropic Efficiency of Compressor
            
            cp_c = gamma_ratio * (8.314 / MW);
            ad_eff_c = (obj.Pr_c^(1 / gamma_ratio) - 1) / (obj.Pr_c^((1 / gamma_ratio) / poly_eff_c) - 1);
            
            obj.p03 = obj.Pr_c * obj.p02;
            obj.T03 = obj.T02 * (1 + 1 / ad_eff_c * (obj.Pr_c^(1 / gamma_ratio) - 1));
            obj.w_c = cp_c * (obj.T03 - obj.T02);
        end

        function obj = fuel_pump(obj)
            rho_f = 780;  % kg/m^3
            delta_p_f1 = 20.7; % kPa
            delta_p_inj = 540; % kPa

            ad_eff_p = 0.45;

            obj.p_f1 = obj.pa + delta_p_f1;
            obj.p_f2 = delta_p_inj + obj.p03;

            delta_pump = obj.p_f2 - obj.p_f1;

            obj.w_p = (obj.f+obj.f_ib+obj.f_ab)*delta_pump/(ad_eff_p*rho_f);

        end

        function obj = combustor_main(obj)
            obj.Pr_b = 0.95;
            MW = 0.0289; % Molecular weight of gas in combustor (kg/mol)
            eff_b = 0.99; % Efficiency of Main Combustor

            Tmax_0 = 1500;
            bmax = 0.1;
            Cbl = 500;
            n = 0.6;
            
            cp_b = (8.314 / MW) * (3.70 + 0.66 * ((obj.T03 / 1000)^2) - 0.20 * ((obj.T03 / 1000)^3));
            delta_hr = 43.52 * 10^6; % Heating Rate of Fuel (J/kg)
            
            obj.p04 = obj.Pr_b * obj.p03;
            obj.T04 = (obj.f * (eff_b * delta_hr / cp_b) + obj.T03 * (1 - obj.b)) / (1 - obj.b + obj.f);

            obj.Tmax_t = Tmax_0 + Cbl * (obj.b / bmax)^n;
            obj.f_max_main = ((1 - obj.b)* (obj.Tmax_t/obj.T03 - 1))/ ((eff_b*delta_hr)/(cp_b*obj.T03) - (obj.Tmax_t/obj.T03));
        end
        
        function obj = turbine(obj)
            MW = 0.0289; % Molecular weight of gas in turbine (kg/mol)
            poly_eff_t = 0.94; % Polytropic Efficiency of Turbine

            obj.w_t = obj.w_c+obj.w_p;
            
            cp_t = (8.314 / MW) * (3.38 + 0.70 * (obj.T04 / 1000)^2 - 0.20 * (obj.T04 / 1000)^3);
            obj.T051 = obj.T04 - obj.w_t / (cp_t * (1 - obj.b + obj.f));
            Tr_t = obj.T051 / obj.T04;
            ad_eff_t = (Tr_t - 1) / (Tr_t^(1 / poly_eff_t) - 1);
            obj.p051 = obj.p04 * (1 - 1 / ad_eff_t * (1 - Tr_t))^(cp_t /(8.314/MW));
        end
        
        function obj = mixer_turbine(obj)
            cp_tm_n = (3.43 + 0.72 * (obj.T051 / 1000)^2 - 0.21 * (obj.T051 / 1000)^3); % Normalized cp for turbine mixer (cp/R)
            
            obj.T051m = obj.T051 + obj.b * (obj.T03 - obj.T051) / (1 + obj.f);
            obj.p051m = obj.p051 * ((obj.T051m / obj.T051) / ((obj.T03 / obj.T051)^(obj.b / (1 + obj.f))))^cp_tm_n;
        end
        
        function obj = combustor_interturbine(obj)
            obj.Pr_ib = 0.94; % Interturbine Combustor Pressure Ratio
            MW = 0.0289; % Molecular weight of gas in Interturbine Combustor (kg/mol)
            eff_ib = 0.99; % Efficiency of Interturbine Combustor
            
            cp_ib = (8.314 / MW) * (3.82 + 0.32 * (obj.T051m / 1000)^2 - 0.060 * (obj.T051m / 1000)^3);
            delta_hr = 43.52 * 10^6; % Heating Rate of Fuel (J/kg)
            
            obj.p0514 = obj.Pr_ib * obj.p051m;
            obj.T0514 = (obj.f_ib * (eff_ib * delta_hr / cp_ib) + (1 + obj.f) * obj.T051m) / (1 + obj.f + obj.f_ib);
            obj.T_max_ib = 1500;
            obj.f_max_ib = ((1 + obj.f )* (obj.T_max_ib/obj.T051m - 1))/ ((eff_ib*delta_hr)/(cp_ib*obj.T051m) - (obj.T_max_ib/obj.T051m));
        end
        
        function obj = turbine_fan(obj)
            MW = 0.0289; % Molecular weight of gas in turbine (kg/mol)
            poly_eff_ft = 0.94; % Polytropic Efficiency of Turbine
            
            cp_ft = (8.314 / MW) * (3.40 + 0.63 * (obj.T0514 / 1000)^2 - 0.18 * (obj.T0514 / 1000)^3);
            cp_ft_n = 3.40 + 0.63*(obj.T0514 / 1000)^2 - 0.18*(obj.T0514 / 1000)^3;
            
            obj.T052 = obj.T0514 - obj.w_f / (cp_ft*(1+obj.f+obj.f_ib));
            
            Tr_ft = obj.T052 / obj.T0514;
            ad_eff_ft = (Tr_ft - 1) / ((Tr_ft^(1/poly_eff_ft)) - 1);
            
            obj.p052 = obj.p0514 * (1 - 1/ad_eff_ft * (1 - Tr_ft))^(cp_ft_n);
        end
        
        function obj = combustor_afterburner(obj)
            obj.Pr_ab = 0.97; % Interturbine Combustor Pressure Ratio
            MW = 0.0289; % Molecular weight of gas in Interturbine Combustor (kg/mol)
            eff_ab = 0.96; % Efficiency of Interturbine Combustor
            
            cp_ab = (8.314 / MW) * (3.50 + 0.72 * (obj.T052 / 1000)^2 - 0.21 * (obj.T052 / 1000)^3);
            delta_hr = 43.52 * 10^6; % Heating Rate of Fuel (J/kg)
            
            obj.p06 = obj.Pr_ab * obj.p052;
            obj.T06 = (obj.f_ab * (eff_ab * delta_hr / cp_ab) + (1 + obj.f + obj.f_ib) * obj.T052) / (1 + obj.f + obj.f_ib + obj.f_ab);
            obj.T_max_ab = 2300;
            obj.f_max_ab = ((1 + obj.f + obj.f_ib)* (obj.T_max_ab/obj.T052 - 1))/ ((eff_ab*delta_hr)/(cp_ab*obj.T052) - (obj.T_max_ab/obj.T052));
        end
        
        function obj = nozzle_core(obj)
            MW = 0.0289;

            ad_eff_n = 0.96;
            cp_n_n = 3.45 + 0.55 * (obj.T06 / 1000)^2 - 0.15 * (obj.T06 / 1000)^3; % Normalized cp for core nozzle (cp/R)
            cp_n = cp_n_n*(8.314/MW);
            
            obj.pe = obj.pa;
            obj.Te = obj.T06 * (1 - ad_eff_n * (1 - (obj.pe / obj.p06)^(1 / cp_n_n)));
            
            obj.ue = sqrt(2*cp_n*(obj.T06-obj.Te));
            gamma_n = cp_n_n/(cp_n_n-1);
            obj.Me = obj.ue/sqrt(gamma_n*(8.314/MW)*obj.Te);
        end
        
        function obj = nozzle_fan(obj)
            MW = 0.0289;
            ad_eff_fn = 0.97;
            cp_fn_n = 3.5; % Normalized cp for fan nozzle (cp/R)
            cp_fn = cp_fn_n*(8.314/MW);
            
            obj.pef = obj.pa;
            obj.Tef = obj.T02 * (1 - ad_eff_fn * (1 - (obj.pef / obj.p02)^(1 / cp_fn_n)));

            obj.uef = sqrt(2*cp_fn*(obj.T02-obj.Tef));
            gamma_fn = cp_fn_n/(cp_fn_n-1);
            obj.Mef = obj.uef/sqrt(gamma_fn*(8.314/MW)*obj.Tef);
        end
        
        function obj = mixer_nozzle(obj)
            C_nm = 0.029;
            Pr_nm = 1 - C_nm * obj.beta;
            
            obj.fuel_sum = 1 + obj.f + obj.f_ib + obj.f_ab;
            
            obj.T07 = (obj.beta * obj.T02 + obj.T06 * obj.fuel_sum) / (obj.fuel_sum + obj.beta);
            gamma_nm = 1.44 - 0.139 * (obj.T07 / 1000) + 0.0357 * (obj.T07 / 1000)^2 - 0.004 * (obj.T07 / 1000)^3;
            gamma_ratio = gamma_nm / (gamma_nm - 1);
            
            obj.p07 = Pr_nm * (obj.T07^gamma_ratio) * (obj.p06^(obj.fuel_sum / (obj.fuel_sum + obj.beta))) * ...
                      (obj.p02^(obj.beta / (obj.fuel_sum + obj.beta))) / (obj.T06^(gamma_ratio * obj.fuel_sum / (obj.fuel_sum + obj.beta))) / ...
                      (obj.T02^(obj.beta * gamma_ratio / (obj.fuel_sum + obj.beta)));
        end
        
        function obj = nozzle_combined(obj)
            MW = 0.0289;
            ad_eff_nc = 0.96;
            cp_nc_n = 3.45 + 0.55 * (obj.T07 / 1000)^2 - 0.15 * (obj.T07 / 1000)^3; % Normalized cp for core nozzle (cp/R)
            
            obj.pec = obj.pa;
            obj.Tec = obj.T07 * (1 - ad_eff_nc * (1 - (obj.pec / obj.p07)^(1 / cp_nc_n)));

            obj.uec = sqrt(2*cp_nc_n*(8.314/MW)*(obj.T07 - obj.Tec));
            gammaec = cp_nc_n/(cp_nc_n - 1);
            obj.Mec = obj.uec / sqrt(gammaec*(8.314/MW)*obj.Tec);
        end

        function obj = get_performance_parameters_cn(obj)
            delta_hr = 43.52*10^6;
            C_b1 = 263;
            MW = 0.0289;
            u = obj.M*sqrt(1.4*(8.314/MW)*obj.Ta);
            drag = C_b1*(obj.M^2)*(obj.pa/101.325)*(obj.beta^(1.5));

            obj.ST_cn = ((1+ obj.beta+obj.f + obj.f_ib + obj.f_ab)*obj.uec - (1+obj.beta)*u-drag)/1000;
            obj.TSFC_cn = ((obj.f + obj.f_ib + obj.f_ab)/obj.ST_cn)*3600;
            obj.n0_cn = ((1/(obj.TSFC_cn/3600))*(u/delta_hr))*1000;
            obj.np_cn = (2*obj.ST_cn*1000*u)/((1+obj.beta+obj.f+obj.f_ib+obj.f_ab)*obj.uec^2-(1+obj.beta)*u^2);
            obj.nth_cn = obj.n0_cn/obj.np_cn;
        end

        function obj = get_performance_parameters_sn(obj)
            delta_hr = 43.52*10^6;
            C_b1 = 263;
            MW = 0.0289;
            u = obj.M*sqrt(1.4*(8.314/MW)*obj.Ta);
            drag = C_b1*(obj.M^2)*(obj.pa/101.325)*(obj.beta^(1.5));

            obj.ST_sn = ((1 + obj.f + obj.f_ib + obj.f_ab)*obj.ue - (1+obj.beta)*u+obj.beta*obj.uef-drag)/1000;
            obj.TSFC_sn = ((obj.f + obj.f_ib + obj.f_ab)/obj.ST_sn)*3600;
            obj.n0_sn = ((1/(obj.TSFC_sn/3600))*(u/delta_hr))*1000;
            obj.np_sn = (2*obj.ST_sn*1000*u)/((1 + obj.f + obj.f_ib + obj.f_ab)*obj.ue^2 - (1+obj.beta)*u^2+obj.beta*obj.uef^2);
            obj.nth_sn = obj.n0_sn/obj.np_sn;
        end
        
    end
end
function [P_ch_avg_kPa,P_ch_avg_psi,thrust_avg_N,thrust_avg_lb,m_dot_avg_ox,m_dot_avg_fuel,m_dot_avg,OF_avg,impulse_sp,impulse_tot_Ns,L_c_m,P_delta_avg_kPa,P_delta_avg_psi,P_drop_kPa,P_drop_avg_kPa,c_star_mps,c_star_avg_mps,T_th_K,T_th_avg_K,gamma_avg,MW_avg_g,R_sp_avg_J,C_p_prop_avg,burn_time,net_efficiency,FUEL_P,P_ch_kPa ,T_ch_K,T_ex_K,V_ex,P_ex_kPa,M_ex,gamma,MW_g,R_sp_J,CdA_OX,CdA_FUEL,OF,m_dot_ox ,m_dot_fuel,m_dot,OXGEN_mass,FGEN_mass,OX_Status,burn_status,t,P_ch_psi ,P_delta_kPa ,P_delta_psi,P_ex_psi,impulse_Ns,impulse_del_Ns,thrust_N,thrust_lb,M_ex_avg,T_ex_avg_K,T_ch_avg_K] = simulate_burn(dt,P_amb_Pa,A_ex_m,c_star_eff,noz_eff,dP,P_ox_kPa,OX_mass_kg,FUEL_mass_kg,OX_rho_kg,OX_D_orifice_in,OX_N_orifice,OX_OD_annulus_in,OX_ID_annulus_in,OX_Discharge_Coeff,FUEL_rho_kg,FUEL_D_orifice1_in,FUEL_N_orifice1,FUEL_D_orifice2_in,FUEL_N_orifice2,FUEL_OD_annulus_in,FUEL_ID_annulus_in,FUEL_Discharge_Coeff,AeAt,A_th_m,D_ch_m,L_ch_m,g)
%% VARIABLES - BASED ON INPUTS
net_efficiency=c_star_eff*noz_eff; %[-] - Overall Efficiency
FUEL_P(1)=P_ox_kPa(1)-dP/(6894.76*1000);%[kPa] - Fuel Initial Pressure
%% CONSTANTS
t(1)=0;%[s] - Initial Time
P_ch_kPa(1)=101.3;%[kPa] - Chamber Initial Pressure
T_ch_K=[500 1000];%[K] - Chamber Initial Temperatures
OX_Status="Liquid"; % - Initial Oxidizer State
M_ex=[0.5,1.0]; %[-] - Initial Mach Number at Exit
gamma(1:2)=1.2; % - Gamma
MW_g(1:2)=28.9647; % - Initial (1-2) Molar Weight
R_sp_J=8314./MW_g; % - Initial (1-2) Gas Constant
%% OXIDIZER INJECTOR AREAS AND CdA
% Calculations - Orifice and Annulus Areas, CdA_OX
OX_A_orifice = (pi*(OX_D_orifice_in/2)^2)/1550.003; % [m^2]
OX_A_annulus = (pi*(OX_OD_annulus_in/2)^2-pi*(OX_ID_annulus_in/2)^2)/1550.003; % [m^2]

CdA_OX = (OX_A_orifice*OX_N_orifice+OX_A_annulus)*OX_Discharge_Coeff; % [m^2]
%% FUEL INJECTOR AREAS AND CdA
% Calculations - Orifice, Annulus Areas, and CdA_FUEL
FUEL_A_orifice1 = (pi*(FUEL_D_orifice1_in/2)^2)/1550.003; % [m^2]
FUEL_A_orifice2 = (pi*(FUEL_D_orifice2_in/2)^2)/1550.003; % [m^2]
FUEL_A_annulus = (pi*(FUEL_OD_annulus_in/2)^2-pi*(FUEL_ID_annulus_in/2)^2)/1550.003; % [m^2]

CdA_FUEL = (FUEL_A_orifice1*FUEL_N_orifice1+FUEL_A_orifice2*FUEL_N_orifice2+FUEL_A_annulus)*FUEL_Discharge_Coeff; % [m^2]
%% INITIAL VALUES
i=1;
% Calculate Initial Flowrates based on CHAMBER_P(1)=0[kPa]
m_dot_ox(i) = CdA_OX*sqrt(2*OX_rho_kg*(P_ox_kPa(i)-P_ch_kPa(i))*1000); % [kg/s]
m_dot_fuel(i) = CdA_FUEL*sqrt(2*FUEL_rho_kg*(FUEL_P(i)-P_ch_kPa(i))*1000); % [kg/s]
m_dot(i)=m_dot_ox(i)+m_dot_fuel(i);

OF(i)=m_dot_ox(i)/m_dot_fuel(i);

% Calculate Initial OXGEN_mass based on OX_flowrate(1)
OXGEN_mass(i)=m_dot_ox(i)*dt;
FGEN_mass(i)=m_dot_fuel(i)*dt;

% Initalize Loop Control
i=2;
while(OXGEN_mass(i-1)>0 && FGEN_mass(i-1)>0)
    %% TIME
    t(i)=t(i-1)+dt;
    %% gam0, MW, R
    if(i>2)
        % Calculate Constant Values from Lookup
        gamma(i)=calc_gamma(P_ch_kPa(i-1),OF(i-1));
        MW_g(i)=calc_MW_g(P_ch_kPa(i-1),OF(i-1));
        R_sp_J(i)=8314/MW_g(i);
    end
    %% OX_status
    if((OX_mass_kg(i-1)/OX_mass_kg(1))>0.1)
        OX_Status(i)="Liquid";
    else
        OX_Status(i)="Gas";
    end

    if(OX_mass_kg(i-1)>0 && FUEL_mass_kg(i-1)>0)
        P_ch_kPa(i)=m_dot(i-1)*sqrt(T_ch_K(i-1))/(A_th_m*sqrt(gamma(i)/R_sp_J(i))*((gamma(i)+1)/2)^(-(gamma(i)+1)/(2*(gamma(i)-1))))/1000;
    end
    %% CHAMBER_temp (T0)
    if(i>2)
        T_ch_0_K=calc_T_ch_0_K(OF(i-1),P_ch_kPa(i-1));%[K] - Chamber Temperature
        % Calculates by which factor to multiply the result, depending on if
        % the previous pressure is greater than the current, or not.
        if(P_ch_kPa(i)<P_ch_kPa(i-1))
            factor = P_ch_kPa(i)/P_ch_kPa(i-1);
        else
            factor = 1;
        end
        %T_ch_K(i)=c_star_eff*T_ch_0_K*factor;
        T_ch_K(i)=T_ch_0_K;
    end
    %% OX_P
    if(strcmp(OX_Status(i-1),"Liquid"))
        if(OX_mass_kg(i-1)>0 && FUEL_mass_kg(i-1)>0)
            P_ox_kPa(i)=P_ox_kPa(1)*(OX_mass_kg(i-1)/OX_mass_kg(1)*0.3+0.7); % [kPa]
        end
    elseif(strcmp(OX_Status(i-1),"Gas"))
        P_ox_kPa(i)=OX_mass_kg(i-1)/exp(0.05*(t(i)-t(i-1))); % [kPa]
    end
    %% FUEL_P
    if(OX_mass_kg(i-1)>0 && FUEL_mass_kg(i-1)>0)
        if(P_ox_kPa(i)-dP/6.895<0)
            FUEL_P(i)=0;%[kPa]
        else
            FUEL_P(i)=P_ox_kPa(i)-dP/6.895; % [kPa]
        end
    end
    %% OX,FUEL,NET_flowrate
    % Calculate Step Flowrates based on CHAMBER_P(i)
    if(FUEL_P(i)>0 && P_ox_kPa(i)>0)
        m_dot_ox(i) = CdA_OX*sqrt(2*OX_rho_kg*(P_ox_kPa(i)-P_ch_kPa(i))*1000); % [kg/s]
        m_dot_fuel(i) = CdA_FUEL*sqrt(2*FUEL_rho_kg*(FUEL_P(i)-P_ch_kPa(i))*1000); % [kg/s]
    else
        m_dot_ox(i)=0;
        m_dot_fuel(i)=0;
    end
    m_dot(i)=m_dot_ox(i)+m_dot_fuel(i);
    %% MR
    % Calculate Step MR
    if(m_dot_ox(i)==0 || m_dot_fuel(i)==0)
        if(OF(i-1) ~= 0)
            OF(i)=0; %[-]
        end
    else
        OF(i)=m_dot_ox(i)/m_dot_fuel(i); %[-]
    end
    %% OX_mass, OXGEN_mass
    % Calculate Step OX_mass and OXGEN_mass
    if((OX_mass_kg(i-1)-OXGEN_mass(i-1))<0)
        OX_mass_kg(i)=0; %[kg]
    else
        OX_mass_kg(i)=OX_mass_kg(i-1)-OXGEN_mass(i-1); %[kg]
    end
    if(OX_mass_kg(i)==0)
        OXGEN_mass(i)=0; %[kg]
    else
        OXGEN_mass(i)=m_dot_ox(i)*dt;%[kg]
    end
    %% FUEL_mass, FGEN_mass
    % Calculate Step FUEL_mass and FGEN_mass
    if((FUEL_mass_kg(i-1)-FGEN_mass(i-1))<0)
        FUEL_mass_kg(i)=0;%[kg]
    else
        FUEL_mass_kg(i)=FUEL_mass_kg(i-1)-FGEN_mass(i-1);%[kg]
    end
    if(FUEL_mass_kg(i)==0)
        FGEN_mass(i)=0;
    else
        FGEN_mass(i)=m_dot_fuel(i)*dt;
    end
    i=i+1;
end
%% OUTSIDE LOOP CONTROL CALCULATIONS
P_ch_kPa(end)=P_amb_Pa/1000;%[kPa]
P_ch_psi=P_ch_kPa/6.895;%[psi]

for i=1:length(P_ch_kPa)
    if(P_ox_kPa(i)-P_ch_kPa(i)>0)
        P_delta_kPa(i)=P_ox_kPa(i)-P_ch_kPa(i);%[kPa]
    else
        P_delta_kPa(i)=0;
    end
end

P_delta_psi=P_delta_kPa/6.895;%[psi]

% Calculate exit mach matrix using gam0 values
M_ex(3)=net_efficiency*calc_M_ex(AeAt,gamma(3));
for i=4:length(gamma)
    M_ex(i)=calc_M_ex(AeAt,gamma(i));%[-]
end

% Calculate Exit Pressure
P_ex_kPa=[0,P_ch_kPa(2)/10.0];%[kPa]
for i=3:length(P_ch_kPa)
    P_ex_kPa(i) = P_ch_kPa(i)*(1+(gamma(i)-1)/2*M_ex(i)^2)^(-gamma(i)/(gamma(i)-1));%[kPa]
end
P_ex_psi=P_ex_kPa/6.895;%[psi]

% Calculate Exit Temperature
for i=1:length(T_ch_K)
    T_ex_K(i)=T_ch_K(i)/(1+(gamma(i)-1)/2*M_ex(i)^2);%[K]
end

% Calculate Exit Velocity
for i=1:length(M_ex)
    V_ex(i)=M_ex(i)*sqrt(gamma(i)*R_sp_J(i)*T_ex_K(i));%[m/s]
end

% Calculate Thrust
thrust_N=[m_dot(1)*V_ex(1),m_dot(2)*V_ex(2)];
for i=3:length(V_ex)
    if(m_dot(i)*V_ex(i)+(P_ex_kPa(i)*1000-P_amb_Pa)*A_ex_m > 0)
        thrust_N(i)=m_dot(i)*V_ex(i)+(P_ex_kPa(i)*1000-P_amb_Pa)*A_ex_m;%[N]
    else
        thrust_N(i)=0;
    end
end
thrust_lb=thrust_N*0.2248090795;

% Calculate Impulse
impulse_Ns=thrust_N.*dt;

% Calculate Delivered Impulse
impulse_del_Ns(1)=0;
for i=2:length(impulse_Ns)
    impulse_del_Ns(i)=impulse_del_Ns(i-1)+impulse_Ns(i-1);
end

% Calculate Fuel/OX Status
burn_status=["Burning"];
for i=2:length(FUEL_mass_kg)
    if(strcmp(burn_status(i-1),"FD") || strcmp(burn_status(i-1),"OD"))
        burn_status(i)=burn_status(i-1);
    else
        if(FUEL_mass_kg(i)==0)
            burn_status(i)="FD";
        else
            if(OX_mass_kg(i)==0)
                burn_status(i)="OD";
            else
                burn_status(i)="Burning";
            end
        end
    end
end

P_ch_avg_kPa=mean(P_ch_kPa(7:14));%[kPa]
P_ch_avg_psi=P_ch_avg_kPa/6.895;%[psi]

thrust_avg_N=mean(thrust_N(7:14));%[N]
thrust_avg_lb=thrust_avg_N*0.2248090795;%[lbf]

m_dot_avg_ox=mean(m_dot_ox(7:14));%[kg/s]
m_dot_avg_fuel=mean(m_dot_fuel(7:14));%[kg/s]
m_dot_avg=mean(m_dot(7:14));%[kg/s]

OF_avg=mean(OF(7:14));%[-]

impulse_sp=thrust_avg_N/(m_dot_avg*g);
impulse_tot_Ns=sum(impulse_Ns,"all");

L_c_m=(((pi*D_ch_m^2)/4*L_ch_m))/A_th_m; % characteristic length

P_delta_avg_kPa=mean(P_delta_kPa(7:14));%[kPa]
P_delta_avg_psi=mean(P_delta_psi(7:14));%[psi]
P_drop_kPa=P_delta_kPa./P_ch_kPa;%[-]
P_drop_avg_kPa=mean(P_drop_kPa(7:14));%[-]

c_star_mps=(P_ch_kPa*10^3).*A_th_m./m_dot;
c_star_avg_mps=mean(c_star_mps(7:14));

T_th_K=(2./(gamma+1)).*T_ch_K; % [K], Estimates throat temp with gamma
T_th_avg_K=mean(T_th_K(7:14));

gamma_avg=mean(gamma(7:14));%[-]
MW_avg_g=mean(MW_g(7:14));%[g/mol]
R_sp_avg_J=mean(R_sp_J(7:14));%[J/kg*K]
C_p_prop_avg=(gamma_avg/(gamma_avg-1))*(R_sp_avg_J);

M_ex_avg=mean(M_ex(7:14));
T_ex_avg_K=mean(T_ex_K(7:14));
T_ch_avg_K=mean(T_ch_K(7:14));

count=0;
for i=1:length(burn_status)
    if(burn_status(i)=="Burning")
        count=count+1;
    end
end
burn_time=count*dt;%[s]


end

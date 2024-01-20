function [T_max_wall_K,T_melt_wall_K,rho_wall_kg,mat_name] = heat_transfer(delta_m,D_ch_m,D_th_m,T_guess_K,T_th_avg_K,T_amb_K,P_ch_avg_kPa,c_star_avg_mps,C_p_prop_avg,burn_time,mu_prop_Pa,mat_type)
%% CALCULATE HEAT TRANSFER COEFFICIENT
T_avg_K=(T_th_avg_K+T_guess_K)/2; 

h_in=(0.026/(D_th_m^0.2))*(((P_ch_avg_kPa*10^3)/c_star_avg_mps)^0.8)*(C_p_prop_avg*mu_prop_Pa^0.2)*((T_th_avg_K/T_avg_K)^(0.8-0.2*0.95));%[W/m^2/K]
h_ex=0;%[W/m^2/K] Assuming No Loss to Atmosphere
%% DETERMINE MATERIAL PROPERTIES
% K_wall: Thermal Conductivity
% C_p: Specific Heat Capacity
% rho_wall_kg: Density
% alpha_wall: Thermal Diffusivity
switch mat_type
    case 1 % Aluminum 6061
        mat_name = "Aluminum 6061";
        K_wall = 152; %[W/(m*K)] at 300K
        C_p_wall = 897; %[J/(kg*K)]
        rho_wall_kg = 2720; %[kg/m^3]
            alpha_wall =  K_wall/(rho_wall_kg*C_p_wall); %[m^2/s]
        T_melt_wall_K = 855.15;%[K]
    case 2 % Copper 110
        mat_name = "Copper 110";
        K_wall = 390.8; %[W/(m*K)]
        C_p_wall = 385; %[J/(kg*K)]
        rho_wall_kg = 8910; %[kg/m^3]
            alpha_wall =  K_wall/(rho_wall_kg*C_p_wall); %[m^2/s]
        T_melt_wall_K = 1338.15; %[K]
    case 3 % Stainless Steel 304
        mat_name = "Stainless Steel 304";
        K_wall = 21.5; %[W/(m*K)] at 100C
        C_p_wall = 500; %[J/(kg*K)]
        rho_wall_kg = 8000; %[kg/m^3]
            alpha_wall =  K_wall/(rho_wall_kg*C_p_wall); %[m^2/s]
        T_melt_wall_K = 1673.15;%[K]
    case 4 % Steel 1018
        mat_name = "Steel 1018";
        K_wall = 51.9; %[W/(m*K)]
        C_p_wall = 486; %[J/(kg*K)]
        rho_wall_kg = 7870; %[kg/m^3]
            alpha_wall =  K_wall/(rho_wall_kg*C_p_wall); %[m^2/s]
        T_melt_wall_K = 1693.15;%[K]
    case 5 % Inconel 625
        mat_name = "Inconel 625";
        K_wall = 9.80; %[W/(m*K)]
        C_p_wall = 410; %[J/(kg*K)]
        rho_wall_kg = 8440; %[kg/m^3]
            alpha_wall =  K_wall/(rho_wall_kg*C_p_wall); %[m^2/s]
        T_melt_wall_K = 1563.15;%[K]
    otherwise
        input("Enter custom material name: ", mat_name);
        fprintf("\n");
        input("Enter custom value for K: ",K);
        fprintf("\n");
        input("Enter custom value for alpha_wall: ", alpha_wall);
        fprintf("\n");
        input("Enter custom value for melting temperature: ", T_melt_wall_K);
        fprintf("\n");
end
    %% FINITE DIFFERENCE SET UP
    %numerical Inputs
    intervals = 10; %number of intervals to divide the wall into
    
    nodes = intervals+1;
    del_R = delta_m/(intervals);
    
    %The following terms are defined to reduce clutter in the following for
    %loops
    A = D_ch_m/(2*del_R);
    H_a = del_R*h_in/K_wall;
    H_b = del_R*h_ex/K_wall;
    %Note: Coincidentally, H_a and H_b can be interpreted as the numerical Biot
    %numbers
    
    
    %NOTE: "r" is used in MECE-725 (CFD), I'm choosing using s instead to
    %differentiate it from the radial coordinate R
    % s <=0.5 is required for numerical stability
    % s = alpha*del_t/del_R^2
    % this can be interpreted as the numerical Fourier Number
    s_wanted = 0.45;
    
    del_t = s_wanted*del_R^2/alpha_wall;
    steps = round(burn_time/del_t);
    s = alpha_wall*del_t/del_R^2;
    
    T_th = zeros(nodes,steps+1);
    T_th(:,1) = T_amb_K;
    %% F-D CALCULATION
    %Using a simple Explicit Forward-Time, Central Space scheme
    % 1-D Cylindrical Coordinates
    %
    % Governing Equation:
    % 1/alpha * dT/dt = ( d^2 T /dr^2 ) + 1/r * dT/dr
    %
    % Approximating the time derivative with a forward difference
    % dT/dt = ( T(i,n+1) - T(i,n) )/delta_t
    %
    % Approximating the two space derivatives with central differences
    % These are both O(delta_r^2)
    %
    % d^2 T /dr^2 = ( T(i+1,n) - 2 * T(i,n) + T(i-1,n) )/delta_r^2 
    % dT/dr = ( T(i+1,n) - T(i-1,n) ) / (2*delta_r)
    % 
    % The order of error for this is O(delta_t, delta_r^2)
    % 
    
    for n = 1:steps
        for i = 1:nodes
            B = A+i-1;
            switch i
                case 1
                    %case for the interior wall
                    T_th(i,n+1) = 2*s*T_th(i+1,n) + (1-2*s*(1+(1-1/B)*H_a))*T_th(i,n) + 2*s*(1-1/B)*H_a*T_th_avg_K;
                case nodes
                    %case for the exterior wall
                    T_th(i,n+1) = (1-2*s*(1+(1+1/B)*H_b))*T_th(i,n) + 2*s*T_th(i-1,n) + 2*s*(1+1/B)*H_b*T_amb_K;
                    
                otherwise
                    %General case for interior points
                    T_th(i,n+1) = s*(1+1/B)*T_th(i+1,n) + (1-2*s)*T_th(i,n) + s*(1-1/B)*T_th(i-1,n);
            end
            
        end
    end
    dx = 1000*del_R*linspace(0,nodes,nodes);
    T_final = T_th(:,steps+1);
    HT_table = table(dx',T_final,'VariableNames',{'x_mm','T_degK'});
    HT_throat = T_th(1,:);

    T_max_wall_K=HT_table{1,2};
end


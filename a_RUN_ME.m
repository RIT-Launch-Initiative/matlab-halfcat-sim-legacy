%% INTRODUCTION
% Title:             ISOPROPANOL & NITROUS OXIDE ENGINE BUILDER & SIMULATOR
% Developer:         Evan Olsen
% Date:              11/13/2023
%
% Project Iteration: Nitron I
% Team:              RIT Liquid Propulsion
clc; clear; close all; tic;
%% INPUTS --------------------------------------------------------- GENERAL
engine_name = "Nitron I";
team_name = "RIT Liquid Propulsion";
dt=0.025; % --------------------------------------------------------------- [s] Timestep
dx_mm=0.2; % -------------------------------------------------------------- [mm] Space Between Points On Output
figure_style = 1; % ------------------------------------------------------- [-] Figure Color Scheme
                                                                            % 1 = Light Mode
                                                                            % 2 = Dark Mode
                                                                            % 3 = Blueprint
grid_style = 4; % --------------------------------------------------------- [-] Grid Style
                                                                            % 1 = Major Only
                                                                            % 2 = Minor Only
                                                                            % 3 = Major & Minor
                                                                            % 4 = Off
%% INPUTS --------------------------------------------------- HEAT TRANSFER
mat_type=3; % ------------------------------------------------------------- [-] Material Type
                                                                            % 1 = Aluminum 6061
                                                                            % 2 = Copper 110
                                                                            % 3 = Stainless Steel 304
                                                                            % 4 = Steel 1018
                                                                            % 5 = Inconel 625
                                                                            
T_guess_K=903; % ---------------------------------------------------------- [K] Estimated Final Throat Wall Temperature
mu_prop_Pa=0.000051420000000000006; % ------------------------------------- [Pa*s] Viscosity of Products of Combustion
%% INPUTS ------------------------------------------------- ENGINE GEOMETRY
nozzle_type_int=3; % ------------------------------------------------------ [-] Nozzle Type
                                                                            % 1 = 15-Degree Half-Angle
                                                                            % 2 = Bezier Curve
                                                                            % 3 = Method of Characteristics

delta_in=0.35; % ---------------------------------------------------------- [in] Wall Thickness
D_th_in=0.51; % ------------------------------------------------------------ [in] Throat Diameter
D_ex_in=1.09; % --------------------------------------------------------- [in] Nozzle Exit Diameter
D_ch_in=2; % -------------------------------------------------------------- [in] Combustion Chamber Diameter
noz_eff=0.95; % ----------------------------------------------------------- [%] Nozzle Efficiency
c_star_eff=0.80; % -------------------------------------------------------- [%] C-star Efficiency
%% INPUTS ------------------------------------------------------- CONSTANTS
P_amb_Pa=101352.9; % -------------------------------------------------------- [Pa] Ambient Pressure
T_amb_K=298.15; % -------------------------------------------------------- [K] Ambient Temperature
g=9.81; % ----------------------------------------------------------------- [m/s^2] Accel. Due to Gravity
R_uni_J=8.3144598; % ------------------------------------------------------ [J/mol*K] Universal Gas Constant
%% INPUTS --------------------------------------------- FUEL SPECIFICATIONS
% FUEL TANK PARAMETERS
FUEL_type = "Isopropanol";

FUEL_TANK_OD=4.000; % ------------------------------------------------------ [in] Fuel Tank OD
FUEL_TANK_ID=3.750; % ----------------------------------------------------- [in] Fuel Tank ID
FUEL_TANK_len=6.000; % ------------------------------------------------------- [in] Fuel Tank Length
FUEL_displacement=3.750; % --------------------------------------------------- [in] Fuel Tank Displacement

dP=15; % ------------------------------------------------------------------ [psi] Fuel->OX Pressure Difference
FUEL_rho_kg=786; % -------------------------------------------------------- [kg/m^3] Fuel Density

% FUEL INJECTOR PARAMETERS
FUEL_D_orifice1_in = 0.03130; % -------------------------------------------- [in] Fuel Orifice Diameter #1
FUEL_N_orifice1 = 14; % ---------------------------------------------------- [-] Fuel Orifice Diameter #1 Count
FUEL_D_orifice2_in = 0; % ------------------------------------------------- [in] Fuel Orifice Diameter #2
FUEL_N_orifice2 = 0; % ---------------------------------------------------- [-] Fuel Orifice Diameter #2 Count
FUEL_OD_annulus_in = 0; % ------------------------------------------------- [in] Fuel Annulus OD
FUEL_ID_annulus_in = 0; % ------------------------------------------------- [in] Fuel Annulus ID
FUEL_Discharge_Coeff = 0.25; % -------------------------------------------- [-] Fuel Discharge Coefficient

%% INPUTS ----------------------------------------- OXIDIZER SPECIFICATIONS
% OXIDIZER TANK PARAMETERS
OX_type="Nitrous Oxide [N20]";
TANK_style="Stacked"; % --------------------------------------------------- [-] Tank Style
                                                                            % "Concentric"
                                                                            % "Stacked"

OX_TANK_OD=4.000; % ----------------------------------------------------------- [in] Oxidizer Tank OD
OX_TANK_ID=3.750; % ------------------------------------------------------- [in] Oxidizer Tank ID
OX_displacement=8.5;% ------------------------------------------------------ [in] Oxidizer Tank Displacement
SIPHON_dia=0; % ----------------------------------------------------------- [in] Oxidizer Tank Siphone Diameter

OX_temp=T_amb_K; % -------------------------------------------------------- [K] Oxidizer Temperature
OX_rho_kg=744; % ---------------------------------------------------------- [kg/m^3] Oxidizer Density
P_ox_kPa(1)=5660; % ------------------------------------------------------- [kPa] Initial Oxidizer Pressure

% OXIDIZER INJECTOR PARAMETERS
OX_D_orifice_in = 0.03130; % ------------------------------------------------ [in] Oxidizer Orifice Diameter
OX_N_orifice = 28; % ------------------------------------------------------- [-] Oxidizer Orifice Count
OX_OD_annulus_in = 0; % --------------------------------------------------- [in] Oxidizer Annulus OD
OX_ID_annulus_in = 0; % --------------------------------------------------- [in] Oxidizer Annulus ID
OX_Discharge_Coeff = 0.25; % ---------------------------------------------- [-] Oxidizer Discharge Coefficient

%% RUN PROGRAM
% CALCULATE ADDITIONAL GEOMETRY SPECS
[dx_m,delta_m,nozzle_type_str,D_th_m,A_th_m,D_ex_m,A_ex_m,AeAt,D_ch_m,    ...
 L_ch_m,L_ch_in]=calc_geometry_specs(dx_mm,nozzle_type_int,delta_in,      ...
 D_th_in,D_ex_in,D_ch_in);

% CALCULATE ADDITIONAL FUEL SPECS
[FUEL_TANK_vol,FUEL_vol,FUEL_mass]=calc_fuel_specs(FUEL_TANK_OD,          ...
 FUEL_TANK_len,FUEL_TANK_ID,FUEL_displacement,FUEL_rho_kg);

% CALCULATE ADDITIONAL OX SPECS
[OX_vol_in,OX_usable_vol_in,OX_total_mass_kg,OX_usable_mass_kg,OX_mass_kg]...
 =calc_ox_specs(OX_rho_kg,OX_TANK_ID,OX_displacement,FUEL_displacement,   ...
 FUEL_TANK_OD,FUEL_TANK_len,SIPHON_dia,TANK_style);

% SIMULATE BURN
[P_ch_avg_kPa,P_ch_avg_psi,thrust_avg_N,thrust_avg_lb,m_dot_avg_ox,       ...
 m_dot_avg_fuel,m_dot_avg,OF_avg,impulse_sp,impulse_tot_Ns,L_c_m,         ...
 P_delta_avg_kPa,P_delta_avg_psi,P_drop_kPa,P_drop_avg_kPa,c_star_mps,    ...
 c_star_avg_mps,T_th_K,T_th_avg_K,gamma_avg,MW_avg_g,R_sp_avg_J,          ...
 C_p_prop_avg,burn_time,net_efficiency,FUEL_P,P_ch_kPa ,T_ch_K,T_ex_K,    ...
 V_ex,P_ex_kPa,M_ex,gamma,MW_g,R_sp_J,CdA_OX,CdA_FUEL,OF,m_dot_ox ,       ...
 m_dot_fuel,m_dot,OXGEN_mass ,FGEN_mass,OX_Status,burn_status,t,P_ch_psi, ...
 P_delta_kPa ,P_delta_psi,P_ex_psi,impulse_Ns,impulse_del_Ns,thrust_N,    ...
 thrust_lb,M_ex_avg,T_ex_avg_K,T_ch_avg_K] = simulate_burn(dt,P_amb_Pa,A_ex_m,c_star_eff,       ...
 noz_eff,dP,P_ox_kPa,OX_mass_kg,FUEL_mass,OX_rho_kg,OX_D_orifice_in,      ...
 OX_N_orifice,OX_OD_annulus_in,OX_ID_annulus_in,OX_Discharge_Coeff,       ...
 FUEL_rho_kg,FUEL_D_orifice1_in,FUEL_N_orifice1,FUEL_D_orifice2_in,       ...
 FUEL_N_orifice2,FUEL_OD_annulus_in,FUEL_ID_annulus_in,                   ...
 FUEL_Discharge_Coeff,AeAt,A_th_m,D_ch_m,L_ch_m,g);

% GENERATE GEOMETRY
[engine_contour_m,engine_contour_in,x_tangent_in] = calc_geometry(D_ch_m, ...
 D_ex_m,D_th_m,L_ch_m,dx_m,gamma_avg,P_ch_avg_kPa,P_amb_Pa,T_ex_avg_K,    ...
 T_ch_avg_K,R_sp_avg_J,M_ex_avg,nozzle_type_int);

% 1D TRANSIENT HEAT TRANSFER ANALYSIS @ THROAT
[T_max_wall_K,T_melt_wall_K,rho_wall_kg,mat_name] = heat_transfer(delta_m,...
 D_ch_m,D_th_m,T_guess_K,T_th_avg_K,T_amb_K,P_ch_avg_kPa,c_star_avg_mps,  ...
 C_p_prop_avg,burn_time,mu_prop_Pa,mat_type);

% ENGINE VOLUME & WEIGHT
[engine_volume_m,engine_weight_kg] = calc_volume(engine_contour_m,        ...
 delta_m,rho_wall_kg);

% COMMAND WINDOW OUTPUT
output_command_window(P_ch_avg_psi,thrust_avg_lb,thrust_avg_N,m_dot_avg,  ...
 OF_avg,impulse_sp,P_delta_avg_psi,P_delta_avg_kPa,burn_time,             ...
 impulse_tot_Ns,mat_name,delta_in,engine_weight_kg,T_max_wall_K,          ...
 T_melt_wall_K);

% FIGURE WINDOW OUTPUT
[output_geometry] = output_geometry(engine_contour_in,delta_in,L_ch_in,   ...
 x_tangent_in,D_th_in,engine_name,team_name,nozzle_type_str,figure_style, ...
 grid_style);
[output_simulation] = output_simulation(t,thrust_lb,P_ch_psi,m_dot,       ...
 P_delta_psi,figure_style);

% SAVE FIGURES
savefig(output_simulation,'output_simulation.fig');
savefig(output_geometry,'output_geometry.fig');
toc
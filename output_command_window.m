function [] = output_simulation(P_ch_avg_psi,thrust_avg_lb,thrust_avg_N,m_dot_avg,OF_avg,impulse_sp,P_delta_avg_psi,P_delta_avg_kPa,burn_time,impulse_tot_Ns,mat_name,delta_in,engine_weight_kg,T_max_wall_K,T_melt_wall_K)
fprintf('<strong>PERFORMANCE SPECIFICATIONS></strong>\n')
fprintf('CHAMBER PRESSURE:   %.2f  [psi]\n',P_ch_avg_psi)
fprintf('INITIAL THRUST:     %.2f  [lbs]\n                    %.2f  [N]\nMASS FLOWRATE:      %.3f   [kg/s]\n',thrust_avg_lb,thrust_avg_N,m_dot_avg);
fprintf('MIXTURE RATIO(O/F): %.3f   [-]\nSPECIFIC IMPULSE:   %.2f  [s]\nPRESSURE DELTA:     %.2f  [psi]\n                    %.2f [kPa]\n',OF_avg,impulse_sp,P_delta_avg_psi,P_delta_avg_kPa);
fprintf('BURN TIME (FUEL):   %.2f    [s]\nTOTAL IMPULSE:      %.2f [Ns]\n\n',burn_time,impulse_tot_Ns);
fprintf('<strong>HEAT TRANSFER SPECIFICATIONS></strong>\n')
fprintf('SELECTED MATERIAL:  %s\n', upper(mat_name))
fprintf('WALL THICKNESS:     %.2f    [in]\n',delta_in);
fprintf('WEIGHT:             %.2f    [kg]\n                    %.2f    [lb]\n',engine_weight_kg,engine_weight_kg*2.205);
fprintf('MAX WALL TEMP:      %.2f  [K]\n',T_max_wall_K);
fprintf('MELTING TEMP:       %2.2f [K]\n',T_melt_wall_K);
end


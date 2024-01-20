function [dx_m,delta_m,nozzle_type_str,D_th_m,A_th_m,D_ex_m,A_ex_m,AeAt,D_ch_m,L_ch_m,L_ch_in] = calc_geometry_specs(dx_mm,nozzle_type_int,delta_in,D_th_in,D_ex_in,D_ch_in)
dx_m=dx_mm/1000;%[m]
delta_m=delta_in/39.37;%[m] Wall Thickness

if nozzle_type_int == 1
    nozzle_type_str = "15-Degree Half-Angle";
elseif nozzle_type_int == 2
    nozzle_type_str = "Bezier Curve";
elseif nozzle_type_int == 3
    nozzle_type_str = "Method of Characteristics";
else
    warning("Invalid Nozzle Type!")
end

    D_th_m=D_th_in/39.37;%[m]
    A_th_m=(pi/4)*D_th_in^2/1550;%[m^2]

    D_ex_m=D_ex_in/39.37;%[m]
    A_ex_m=(pi/4)*D_ex_in^2/1550;%[m^2]

    AeAt=A_ex_m/A_th_m;

    D_ch_m=D_ch_in/39.37;%[m]

    L_ch_cm=exp(0.029*log(D_th_in*2.54)^2+0.47*log(D_th_in*2.54)+1.94); %[cm] Chamber Length
    L_ch_in=L_ch_cm/2.54;%[in]
    L_ch_m=L_ch_in/39.37;%[m]
end


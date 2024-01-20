function [engine_volume_m,engine_weight_kg] = calc_volume(engine_contour_m,delta_m,rho_wall_kg)
    [top_curve_x,top_curve_y] = curve_offset(engine_contour_m(1,:),engine_contour_m(2,:),delta_m);
    top_curve=[top_curve_x;top_curve_y];
    bottom_curve = engine_contour_m;
    V_1=0; V_2=0;
    for i=1:length(bottom_curve(1,:))-1
        dx_1=bottom_curve(1,i+1)-bottom_curve(1,i);
        dV_1=pi*bottom_curve(2,i)^2*dx_1;
        V_1=V_1+dV_1;
    end
    for j=1:length(top_curve(1,:))-1
        dx_2=top_curve(1,j+1)-top_curve(1,j);
        dV_2=pi*top_curve(2,j)^2*dx_2;
        V_2=V_2+dV_2;
    end
    engine_volume_m=V_2-V_1;%[m^3] Engine Volume
    engine_weight_kg=engine_volume_m*rho_wall_kg;%[kg] Engine Weight
end


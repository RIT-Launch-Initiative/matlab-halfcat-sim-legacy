function [FUEL_TANK_vol,FUEL_vol,FUEL_mass] = calc_fuel_specs(FUEL_TANK_OD,FUEL_TANK_len,FUEL_TANK_ID,FUEL_displacement,FUEL_rho_kg)
    FUEL_TANK_vol=pi*FUEL_TANK_OD^2*FUEL_TANK_len/4;%[in^3]

    FUEL_vol=pi*FUEL_TANK_ID^2*FUEL_displacement/4;%[in^3]
    FUEL_mass(1)=FUEL_rho_kg*FUEL_vol/61024;%[kg]
end


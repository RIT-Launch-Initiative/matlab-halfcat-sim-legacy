function [OX_vol_in,OX_usable_vol_in,OX_total_mass_kg,OX_usable_mass_kg,OX_mass_kg] = calc_ox_specs(OX_rho_kg,OX_TANK_ID,OX_displacement,FUEL_displacement,FUEL_TANK_OD,FUEL_TANK_len,SIPHON_dia,TANK_style)
if(TANK_style=="Concentric")
    OX_vol_in=pi*OX_TANK_ID^2*OX_displacement/4-pi*FUEL_TANK_OD^2*            ...
        FUEL_TANK_len/4-pi*SIPHON_dia^2*(OX_displacement-FUEL_TANK_len)/4;     % [in^3]
elseif(TANK_style=="Stacked")
    OX_vol_in=pi*OX_TANK_ID^2*OX_displacement/4;
end
    if(TANK_style=="Concentric")
        if(SIPHON_dia==0)
            if(FUEL_TANK_len<OX_displacement/2)
                OX_usable_vol_in=OX_vol_in-pi*FUEL_TANK_OD^2*FUEL_displacement/4;%[in^3]
            else
                OX_usable_vol_in=pi*(OX_TANK_ID^2-FUEL_TANK_OD^2)*OX_displacement/4;%[in^3]
            end
        else
            OX_usable_vol_in=OX_vol_in;%[in^3]
        end
    elseif(TANK_style=="Stacked")
        OX_usable_vol_in=OX_vol_in;%[in^3]
    else
        warning("Invalid Tank Style\n");
    end
    OX_total_mass_kg=OX_vol_in*OX_rho_kg/61024;%[kg]
    OX_usable_mass_kg=OX_usable_vol_in*OX_rho_kg/61024;%[kg]
    OX_mass_kg(1)=OX_usable_mass_kg;%[kg]
end


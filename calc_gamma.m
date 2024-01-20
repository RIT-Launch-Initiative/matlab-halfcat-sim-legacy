function [gam0] = calc_gamma(CHAMBER_P,MR)
    % Iterates within the main control loop within the calcBurn() function,
    % does not need to iterate within itself to calculate the whole gam0
    % matrix for all of burntime, since other values are dependant on it
    % within the main function during the control loop.
    %% INITALIZE
    % Import RAW CEA data
    RAW_gam=[];
    load('RAW_gam.mat');
    
    % Create Table Bounds
    MR_table=0.25:0.25:10;
    CHAMBER_P_PSI_table=50:50:1000;

    % Convert [kPa]->[psi]
    CHAMBER_P_PSI=CHAMBER_P/6.895;
    %% MAIN CODE
    ind=1;
    for i=1:length(RAW_gam)
        if(rem(i,2)==0)
            % if i is even, do nothing
        else
            % if i is odd, add number to filtered_T0 matrix
            filtered_gam0(ind)=RAW_gam(i);
            ind=ind+1;
        end
    end

    for j=1:length(MR_table)
        for i=1:length(CHAMBER_P_PSI_table)
            if(i+(j-1)*20 ~= length(filtered_gam0))
                gam0_table(i,j)=filtered_gam0(i+(j-1)*20);
            end
        end
    end

    CHAMBER_P_PSI_nearest=round(CHAMBER_P_PSI/50)*50;

    ind_P=1;
    for i=1:length(CHAMBER_P_PSI_table)
        if(CHAMBER_P_PSI_table(i)<=CHAMBER_P_PSI_nearest)
            ind_P=i;
        end
    end

    MR_nearest=round(MR/0.25)*0.25;
    ind_MR=1;
    for i=1:length(MR_table)
        if(MR_table(i)<=MR_nearest)
            ind_MR=i;
        end
    end

    gam0=gam0_table(ind_P,ind_MR);
end
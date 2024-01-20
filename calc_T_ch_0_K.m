%% Title Section
% Description: Code is meant to replicate the work done by HalfCatSimv1.2.1
%              in calculating T0 values based on the precalculated MR values.
% Developer: Evan Olsen
% Team: Launch Initiative - Liquid Propulsion
% Project: Nitron I
function [T0]=T_ch_0_K_calc(MR,CHAMBER_P)
    %% Constants
    % Imports the needed matrices full of values which are "constant" relative
    % to calculating T0 values. 
    RAW_T=[];
    load('RAW_T.mat');
   
    CHAMBER_P_PSI=CHAMBER_P/6.895;%[psi]

    % Creates O/F and Pc_psi matrix rounded table values
    MR_table = 0.25:0.25:10;
    CHAMBER_P_PSI_table = 50:50:1000;

    ind=1;
    for i=1:length(RAW_T)
        if(rem(i,2)==0)
            % if i is even, do nothing
        else
            % if i is odd, add number to filtered_T0 matrix
            filtered_T0(ind)=RAW_T(i);
            ind=ind+1;
        end
    end

    for j=1:length(MR_table)
        for i=1:length(CHAMBER_P_PSI_table)
            if(i+(j-1)*20 ~= length(filtered_T0))
                T_chamber_table(i,j)=filtered_T0(i+(j-1)*20);
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

    % OF_vs_Pc_propellant_val grabs a certain value within the matrix for
    % OF_vs_Pc_propellant depending on the rounded values for Pc and OF
    % during a certain iteration.
    T0 = T_chamber_table(ind_P,ind_MR);
end

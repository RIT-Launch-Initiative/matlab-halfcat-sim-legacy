function [EXIT_mach] = calc_M_ex(AeAt,gamma)
    % Generate mach_table
    resolution=100;
    %mach_table=linspace(1,6,resolution);
    mach_table=1:.01:5.99;

    gam0=round(gamma/0.05)*0.05;

    AeAt_table=zeros(1,length(mach_table));
    for i=1:length(mach_table)

        AeAt_table(i)=1/mach_table(i)*((1+0.5*(gam0-1)*mach_table(i)^2)/(1+0.5*(gam0-1)))^((gam0+1)/(2*(gam0-1)));
        if(AeAt_table(i)<=AeAt)
            EXIT_mach=mach_table(i);
        end
    end
end
%% Function Declaration
function [engine_contour_m,engine_contour_in,x_tangent_in] = calc_geometry(D_ch_m,D_ex_m,D_th_m,L_ch_m,dx_m,gamma_avg,P_ch_avg_kPa,P_amb_Pa,T_ex_avg_K,T_ch_avg_K,R_sp_avg_J,M_ex_avg,nozzle_type_int)
%% NOZZLE GENERATION
% ALL MEASUREMENTS IN METERS
if nozzle_type_int == 1
    % 15[deg] half-angle nozzle
    R_th_m = D_th_m/2;
    R_ex_m = D_ex_m/2;

    theta_ex=15; %deg
    x_noz_end=(R_ex_m-R_th_m)/tand(theta_ex);
    
    nozzle_contour(1,:)=0:dx_m:x_noz_end;
    noz_pts=length(nozzle_contour(1,:));

    nozzle_contour(2,:)=linspace(R_th_m,R_ex_m,noz_pts);
elseif nozzle_type_int == 2
    % Rao Approximation - shortened length 15-degree half angle
    R_ex_m=D_ex_m/2;
    R_th_m=D_th_m/2;

    A_th=pi*R_th_m^2;
    A_ex=pi*R_ex_m^2;
    AeAt=A_ex/A_th;

    theta_ex=14;
    theta_ent=22.5;

    % Quadratic Bezier Curve
    N_x = 0.382*R_th_m*cosd(theta_ent - 90); 
    N_y = 0.382*R_th_m*sind(theta_ent - 90) + 0.382*R_th_m + R_th_m; 
    E_x = 0.8*((sqrt(AeAt) - 1)*R_th_m/tand(15)); 
    E_y = sqrt(AeAt)*R_th_m; 
    m1 = tand(theta_ent); 
    m2 = tand(theta_ex); 
    C1 = N_y - m1*N_x; 
    C2 = E_y - m2*E_x; 
    Q_x = (C2 - C1)/(m1 - m2); 
    Q_y = (m1*C2 - m2*C1)/(m1 - m2); 
    
    x_noz_end = (1 - 1).^2*N_x + 2*(1 - 1).*1*Q_x + 1.^2*E_x;
    noz_pts=round(x_noz_end/dx_m);

    t=linspace(0,1,noz_pts);

    x_curve = (1 - t).^2*N_x + 2*(1 - t).*t*Q_x + t.^2*E_x; 
    y_curve = (1 - t).^2*N_y + 2*(1 - t).*t*Q_y + t.^2*E_y; 
    nozzle_contour=[x_curve;y_curve];
elseif nozzle_type_int == 3

% [raw_noz_wall,~,~]=AxSymMLN(M_ex_avg,gamma_avg,1000,20,5,1e-11,1);
% 
% nozzle_contour=transpose(raw_noz_wall)*D_th_m*0.5;

nozzle_contour=csvread('MOC_contour_temp.csv');
nozzle_contour(:,3)=[];
nozzle_contour=transpose(nozzle_contour)*D_th_m*0.5;
end
%% CHAMBER & THROAT CURVES GENERATION
% Calculated Nozzle Curve
x_noz = nozzle_contour(1,:);
y_noz = nozzle_contour(2,:);
% Finding Tangent Circle at Throat
entrant_angle=45; % deg


slope_nozzle = (y_noz(2) - y_noz(1)) / (x_noz(2) - x_noz(1));

r_th=y_noz(1);
r1=0.382*r_th; %[m] FROM RPE

num_iter = 1000;
ent_x=linspace(x_noz(1),x_noz(2),num_iter);
ent_y=linspace(y_noz(1),y_noz(2),num_iter);

theta_end1=3*pi/2+atan(slope_nozzle);
th_circle1_theta=linspace(3*pi/2,theta_end1,num_iter);
th_circle1_y=r1*sin(th_circle1_theta)+y_noz(1)+r1;
th_circle1_x=r1*cos(th_circle1_theta);

for i=length(th_circle1_y):-1:1
    if(th_circle1_y(i)>y_noz(2))
        th_circle1_y(i)=[];
        th_circle1_x(i)=[];
    end
end

ind=1;
dy_th=zeros(1,length(ent_y));
for i=2:num_iter
    dy_th(i)=abs(th_circle1_y(end)-ent_y(i));
    if(dy_th(i)<=dy_th(i-1))
        ind=i;
    end
end
dx_th=th_circle1_x(end)-ent_x(ind);
th_circle1_x=th_circle1_x-dx_th;

r2=1.5*r_th; % [m]
theta_end2=3*pi/2;
th_circle2_theta=linspace((270-entrant_angle)*pi/180,theta_end2,num_iter);
th_circle2_y=r2*sin(th_circle2_theta(1:end-1))+y_noz(1)+r2;
th_circle2_x=r2*cos(th_circle2_theta(1:end-1))-abs(th_circle1_x(1));

nozzle=[x_noz(3:end);y_noz(3:end)];
th_circle1=[th_circle1_x;th_circle1_y];
th_circle2=[th_circle2_x;th_circle2_y];
th_circle=[th_circle2,th_circle1];
% Generating throat exit curve
r_chamber=D_ch_m/2;
dy_ex=r_chamber-th_circle(2,1);
dx_ex=dy_ex/tand(entrant_angle);
exit_x=linspace(th_circle(1,1),th_circle(1,1)-dx_ex,num_iter);
exit_y=linspace(th_circle(2,1),r_chamber,num_iter);
% Generating chamber interior filet
ch_circle_theta=linspace(deg2rad(90-entrant_angle),pi/2,num_iter);
ch_circle_y=r2*sin(ch_circle_theta)+r_chamber-r2;
ch_circle_x=r2*cos(ch_circle_theta);
for i=2:num_iter
    dy_ch(i)=abs(ch_circle_y(1)-exit_y(i));
    if(dy_ch(i)<dy_ch(i-1))
        min=dy_ch(i);
        ind=i;
    end
end
dx_ch=ch_circle_x(1)-exit_x(ind);
ch_circle_x=ch_circle_x-dx_ch;
ch_circle=[flip(ch_circle_x);flip(ch_circle_y)];
exit=[flip(exit_x(1:ind-1));flip(exit_y(1:ind-1))];
% Generating chamber length
len_x=ch_circle(1,1)-L_ch_m;

ch_len_x=len_x:dx_m:ch_circle(1,1);
ch_len_y=linspace(r_chamber,r_chamber,length(ch_len_x));

ch_len=[ch_len_x;ch_len_y];
% outputs
engine_contour_m=[ch_len,ch_circle,exit,th_circle,nozzle];
x_tangent_in=ch_circle(1,1)*39.37;

% Get unique X values and their indices
[~, unique_indices] = unique(engine_contour_m(1, :));

% Create a new matrix with the unique X values and their corresponding Y values
engine_contour_m = engine_contour_m(:, unique_indices);
engine_contour_in=engine_contour_m*39.37;%[in]

%%  Export Full Contour Data to '.csv' in [cm]:
    % A limitation of the python program is that the data it exports is
    % understood to use [cm] as the units. This issue continues into
    % exporting it as a .csv, as Fusion360 only understands the XYZ data
    % using [cm] as units. As long as we are consistent with our units,
    % data transfer from MATLAB->Python->Fusion360 should operate smoothly.
    z_data=zeros(1,size(engine_contour_m,2));
    output_contour=transpose([engine_contour_m(1,:);engine_contour_m(2,:);z_data]);
    
    writematrix(output_contour,'engine_contour_m.csv');
    
    % Outputs file path where the 'contour.csv' file was generated
    filepath = fileparts(which('engine_contour_m.csv'));
    fprintf('''engine_contour_m.csv'' generated at %s\n',filepath);
end


function [output_geometry] = output_geometry(engine_contour_in,delta_in,L_ch_in,x_tangent_in,D_th_in,engine_name,team_name,engine_type_str,figure_style,grid_style)
%%  Axes Configuration
    output_geometry=figure('Name','Engine Schematic','Toolbar','none','Menubar','none','WindowStyle','Docked');
    ax=gca;
    if grid_style == 1 % Major Only
        grid on
        ax.GridAlpha = 0.25;
    elseif grid_style == 2 % Minor Only
        grid off
        grid(ax,'minor')
        ax.MinorGridLineStyle = '-';
        ax.MinorGridAlpha = 0.1;
    elseif grid_style == 3 % Major & Minor
        grid on
        grid(ax,'minor')
        ax.MinorGridLineStyle = '-';
        ax.GridAlpha = 0.25;
        ax.MinorGridAlpha = 0.1;
    elseif grid_style == 4 % Off
        grid off;
    else
        error("Invalid Grid Style\n");
    end

    % FIGURE STYLE
    if figure_style == 1 % Light Mode
        bg_color='w';
        att_color='k';
    elseif figure_style == 2 % Dark Mode
        bg_color='k';
        att_color='w';
    elseif figure_style == 3 % Blueprint
        bg_color='#3057E1';
        att_color='w';
    else
        error("Invalid Figure Style\n");
    end

    ax.Color=bg_color;
    ax.GridColor = att_color;
    ax.MinorGridColor = att_color;

% FIGURE ASPECT RATIO
    set(ax,'position',[0 0 1 1])
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    ax.LineWidth = 1.5;


    % Sets our figure aspect ratio to 1:1:1, which is important to maintain
    % a real-to-life visualization of the engine geometry.
    pbaspect([1,1,1]);
    % Related to previous point on aspect ratio, in order for the aspect
    % ratio to remain true, the axes limits must also be equal.
    edge_offset = (engine_contour_in(1,end)-engine_contour_in(1,1))*0.1;
    x_limits=[engine_contour_in(1,1)-edge_offset,engine_contour_in(1,end)+edge_offset];
    height=0.5*(engine_contour_in(1,end)-engine_contour_in(1,1));
    y_limits=[-height-edge_offset,height+edge_offset];
    xlim(x_limits);ylim(y_limits)
    % Title Block
    x_center=(x_limits(2)+x_limits(1))/2;
    text(x_center,height,engine_name,'FontSize',24,'HorizontalAlignment','center','VerticalAlignment','middle','Color',att_color,'FontName','Courier New');
    text(x_center,height-edge_offset/2,team_name,'FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle','Color',att_color,'FontName','Courier New');
    
    text(x_center,-height+edge_offset/8,engine_type_str,'FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle','Color',att_color,'FontName','Courier New');
    time=string(datetime);
    text(x_center,-height-edge_offset/2,sprintf('%s',time),'FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle','Color',att_color,'FontName','Courier New');
    %%  Create Patch Object
    [offset_x,offset_y]=curve_offset(engine_contour_in(1,:),engine_contour_in(2,:),delta_in);
    right_edge=[offset_x(end);offset_y(end)];
    top_surface=flip([offset_x;offset_y],2);
    left_edge=[offset_x(1) engine_contour_in(1,1);offset_y(1) engine_contour_in(2,1)];
    wall_obj=[engine_contour_in right_edge top_surface left_edge];
%%  Engine Exterior Shading
    hold on
    patch('XData',wall_obj(1,:),'YData',wall_obj(2,:), ...
        'FaceColor',bg_color,'EdgeColor',att_color,'LineStyle',"-");
    patch('XData',wall_obj(1,:),'YData',-wall_obj(2,:), ...
        'FaceColor',bg_color,'EdgeColor',att_color,'LineStyle',"-");
%%  Exterior Hatch
    % Uses open-source function found on mathworks.com to generate
    % hash-lines to emulate a section-view as is seen in Fusion360.
    % Source Files: https://www.mathworks.com/matlabcentral/fileexchange/30733-hatchfill
    engine_wall_patch_pos=patch('XData',wall_obj(1,:),'YData',wall_obj(2,:),'EdgeColor',att_color);
    engine_wall_patch_neg=patch('XData',wall_obj(1,:),'YData',-wall_obj(2,:),'EdgeColor',att_color);
    curve_hatchfill(att_color,engine_wall_patch_pos,'single',70,4,bg_color);
    curve_hatchfill(att_color,engine_wall_patch_neg,'single',-70,4,bg_color);
%% Quiver Annotations
    % Exit Diameter Annotation
    p1=[engine_contour_in(1,end) 0];
    p2=[engine_contour_in(1,end) engine_contour_in(2,end)];
    dist = p1-p2;
    middle_offset=dist(2)/8;
    p1(2)=middle_offset;
    hold on
    quiver(p1(1),p1(2),dist(1),dist(2)-middle_offset,0,att_color,'LineWidth',1,'MaxHeadSize',0.5,'AutoScaleFactor',1);
    p1(2)=-middle_offset;
    hold on
    quiver(p1(1),p1(2),dist(1),-dist(2)+middle_offset,0,att_color,'LineWidth',1,'MaxHeadSize',0.5,'AutoScaleFactor',1);
    D_exit=abs(dist(2)*2);
    text(p1(1),0,sprintf('%.3f [in]',D_exit), ...
        'HorizontalAlignment','center','VerticalAlignment','middle',FontWeight='bold',FontSize=10,Color=att_color);
    % Throat Diameter Annotation
    [~,ind]=min(abs(engine_contour_in(2,:)));
    p1=[engine_contour_in(1,ind) 0];
    p2=[engine_contour_in(1,ind) engine_contour_in(2,ind)];
    dist=p1-p2;
    p1(2)=middle_offset;
    hold on
    quiver(p1(1),p1(2),dist(1),dist(2)-middle_offset,0,att_color,'LineWidth',1,'MaxHeadSize',3,'AutoScaleFactor',1);
    p1(2)=-middle_offset;
    hold on
    quiver(p1(1),p1(2),dist(1),-dist(2)+middle_offset,0,att_color,'LineWidth',1,'MaxHeadSize',3,'AutoScaleFactor',1);
    text(p1(1),0,sprintf('%.3f [in]',D_th_in), ...
        'HorizontalAlignment','center','VerticalAlignment','middle',FontWeight='bold',FontSize=10,Color=att_color);
    % Chamber Diameter Annotation
    p1=[engine_contour_in(1,1) 0];
    p2=[engine_contour_in(1,1) engine_contour_in(2,1)];
    dist=p1-p2;
    p1(2)=middle_offset;
    hold on
    quiver(p1(1),p1(2),dist(1),dist(2)-middle_offset,0,att_color,'LineWidth',1,'MaxHeadSize',0.5,'AutoScaleFactor',1);
    p1(2)=-middle_offset;
    hold on
    quiver(p1(1),p1(2),dist(1),-dist(2)+middle_offset,0,att_color,'LineWidth',1,'MaxHeadSize',0.5,'AutoScaleFactor',1);
    D_chamber=abs(dist(2)*2);
    text(p1(1),0,sprintf('%.3f [in]',D_chamber(1)), ...
        'HorizontalAlignment','center','VerticalAlignment','middle',FontWeight='bold',FontSize=10,Color=att_color); 
    % Wall Thickness Annotation
    p1=[(x_tangent_in/39.37+engine_contour_in(1,1))/2.0 engine_contour_in(2,1)];
    p2=[(x_tangent_in/39.37+engine_contour_in(1,1))/2.0 engine_contour_in(2,1)+delta_in];

    hold on
    quiver(p1(1),p1(2),0,delta_in,0,att_color,'LineWidth',1,'MaxHeadSize',1,'AutoScaleFactor',1);
    quiver(p2(1),p2(2),0,-delta_in,0,att_color,'LineWidth',1,'MaxHeadSize',1,'AutoScaleFactor',1);
    text(p2(1),p2(2)+1.1*delta_in,sprintf('%.3f [in]',delta_in), ...
        'HorizontalAlignment','center','VerticalAlignment','middle',FontWeight='bold',FontSize=10,Color=att_color);
    % Chamber Length Annotation
    if(engine_contour_in(2,1)<engine_contour_in(2,end))
        p1=[engine_contour_in(1,1)+L_ch_in -2*engine_contour_in(2,end)];
        p2=[engine_contour_in(1,1) -2*engine_contour_in(2,end)];
    else
        p1=[engine_contour_in(1,1)+L_ch_in -2*engine_contour_in(2,1)];
        p2=[engine_contour_in(1,1) -2*engine_contour_in(2,1)];
    end
    midpoint=[(p2(1)+p1(1))/2 p1(1)];
    dist=(p1-midpoint);
    middle_offset=dist(1)/3;
    hold on
    p1(1)=midpoint(1)+middle_offset;
    hold on
    quiver(p1(1),p1(2),abs(dist(1))-middle_offset,0,0,att_color,'LineWidth',1,'MaxHeadSize',0.3,'AutoScaleFactor',1);
    p1(1)=midpoint(1)-middle_offset;
    hold on
    quiver(p1(1),p2(2),-abs(dist(1))+middle_offset,0,0,att_color,'LineWidth',1,'MaxHeadSize',0.3,'AutoScaleFactor',1);
    text(engine_contour_in(1,1)+L_ch_in/2.0,p1(2),sprintf('%.3f [in]',L_ch_in), ...
        'HorizontalAlignment','center','VerticalAlignment','middle',FontWeight='bold',FontSize=10,Color=att_color);
    % Engine Length Annotation
    if(engine_contour_in(2,1)<engine_contour_in(2,end))
        p1=[engine_contour_in(1,1) -1.5*engine_contour_in(2,end)];
        p2=[engine_contour_in(1,end) -1.5*engine_contour_in(2,end)];
    else
        p1=[engine_contour_in(1,1) -1.5*engine_contour_in(2,1)];
        p2=[engine_contour_in(1,end) -1.5*engine_contour_in(2,1)];
    end
    midpoint=[(p2(1)+p1(1))/2 p1(2)];
    dist=(p1-midpoint);
    hold on
    p1(1)=midpoint(1)+middle_offset;
    quiver(p1(1),p1(2),abs(dist(1))-middle_offset,0,0,att_color,'LineWidth',1,'MaxHeadSize',0.175,'AutoScaleFactor',1);
    hold on
    p1(1)=midpoint(1)-middle_offset;
    quiver(p1(1),p2(2),-abs(dist(1))+middle_offset,0,0,att_color,'LineWidth',1,'MaxHeadSize',0.175,'AutoScaleFactor',1);
    L_engine=abs(engine_contour_in(1,1)-engine_contour_in(1,end));
    text(engine_contour_in(1,1)+L_engine/2.0,p1(2),sprintf('%.3f [in]',L_engine), ...
        'HorizontalAlignment','center','VerticalAlignment','middle',FontWeight='bold',FontSize=10,Color=att_color);
    hold off
end
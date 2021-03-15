function fPublish_MeasurementProcess( conf_num,...
                                      conf_num_hot,...
                                      measu_file,...
                                      map_env,...
                                      hotspot_current,...
                                      plannedConf_io,...
                                      executedConf_io,...
                                      dir_ )
%



color_hotspots = [1,1,0;...
                  0,1,0;...
                  1,0,1];

% conf_num_hot
col_num = mod(conf_num_hot,3)+(~mod(conf_num_hot,3)*3);
              
% plannedConf_ind  = selectedConfFusion_sioh(conf_num,2);
% plannedConf_orn  = selectedConfFusion_sioh(conf_num,3);
% executedConf_ind = selectedConfReplanned_sioh(conf_num,2);
% executedConf_orn = selectedConfReplanned_sioh(conf_num,3);

plannedConf_ind  = plannedConf_io(1);
plannedConf_orn  = plannedConf_io(2);
executedConf_ind = executedConf_io(1);
executedConf_orn = executedConf_io(2);


[plannedConf_r,plannedConf_c] = ind2sub(size(map_env),plannedConf_ind);
[executedConf_r,executedConf_c] = ind2sub(size(map_env),executedConf_ind);


% /////////////////////////////////////////
% -- optical beams
% /////////////////////////////////////////

%{
M = dlmread([dir_.TomographyLogs,measu_file]);
start_x__all = M(:,3);
start_y__all = M(:,4);
end___x__all = M(:,5);
end___y__all = M(:,6);
weight___all = M(:,8)/max(M(:,8));

beam_colormap = flipud(bone(512));
for i = 1:size(start_x__all,1)
    plot([start_x__all(i),end___x__all(i)]-0.5,...
         [start_y__all(i),end___y__all(i)]-0.5,...
         '-','Color',beam_colormap(round(255*weight___all(i))+100,:));
end
%}

% /////////////////////////////////////////
% -- publish current hotspot
% /////////////////////////////////////////

plot(hotspot_current(:,1),hotspot_current(:,2),...
    'ok','MarkerSize',5,'MarkerFaceColor',color_hotspots(col_num(1),:));


% /////////////////////////////////////////
% -- planned conf
% /////////////////////////////////////////

% -- sequence numbers
conf_num_dist      = 0.75;
font_size_conf_num = 14;
arrow_length       = 3;
color_planned      = [0.4,0.4,0.4];

quiver(plannedConf_r,plannedConf_c,...
    arrow_length*cosd(plannedConf_orn),...
    arrow_length*sind(plannedConf_orn),...
    'LineWidth',1,...
    'color',color_planned,...
    'MaxHeadSize',8,...
    'Marker','o',...
    'MarkerSize',2);

% conf number.       
x3 = plannedConf_r-(conf_num_dist*cosd(plannedConf_orn));
y3 = plannedConf_c-(conf_num_dist*sind(plannedConf_orn));
%{
                    horizontal
                (right,center,left)        
        <-------------------------------|
                                        |
                                        | vertical
                                        | (bottom,
                                        |  middle,
                                        |  top)
                                        µ 
%}

if     plannedConf_orn >= 315 || plannedConf_orn < 045
    HoriAlign = 'right';
    VertAlign = 'middle';
elseif plannedConf_orn >= 045 && plannedConf_orn < 135
    HoriAlign = 'center';
    VertAlign = 'top'; 
elseif plannedConf_orn >= 135 && plannedConf_orn < 255
    HoriAlign = 'left';
    VertAlign = 'middle';
elseif plannedConf_orn >= 255 && plannedConf_orn < 315
    HoriAlign = 'center';
    VertAlign = 'bottom';
else
    error('WTF');
end

text(x3,y3,num2str(conf_num),...
    'HorizontalAlignment',HoriAlign,...
    'VerticalAlignment',VertAlign,...
    'FontWeight','demi',...
    'FontSize',font_size_conf_num,...
    'Color',color_planned,...
    'Interpreter','latex',...
    'EdgeColor',color_planned,...
    'Margin',1,...
    'FontWeight','bold');


% /////////////////////////////////////////
% -- executed conf
% /////////////////////////////////////////

% -- sequence numbers
conf_num_dist      = 0.5;
arrow_length       = 2.75;

quiver(executedConf_r,executedConf_c,...
    arrow_length*cosd(executedConf_orn),...
    arrow_length*sind(executedConf_orn),...
    'LineWidth',0.75,...
    'color','b',...
    'MaxHeadSize',8,...
    'Marker','o',...
    'MarkerSize',2);

% conf number.       
x3 = executedConf_r-(conf_num_dist*cosd(executedConf_orn));
y3 = executedConf_c-(conf_num_dist*sind(executedConf_orn));
%{
                    horizontal
                (right,center,left)        
        <-------------------------------|
                                        |
                                        | vertical
                                        | (bottom,
                                        |  middle,
                                        |  top)
                                        µ 
%}

if     executedConf_orn >= 315 || executedConf_orn < 045
    HoriAlign = 'right';
    VertAlign = 'middle';
elseif executedConf_orn >= 045 && executedConf_orn < 135
    HoriAlign = 'center';
    VertAlign = 'top'; 
elseif executedConf_orn >= 135 && executedConf_orn < 255
    HoriAlign = 'left';
    VertAlign = 'middle';
elseif executedConf_orn >= 255 && executedConf_orn < 315
    HoriAlign = 'center';
    VertAlign = 'bottom';
else
    error('WTF');
end

text(x3,y3,num2str(conf_num),...
    'HorizontalAlignment',HoriAlign,...
    'VerticalAlignment',VertAlign,...
    'FontWeight','demi',...
    'FontSize',font_size_conf_num,...
    'Color','b',...
    'Interpreter','latex',...
    'EdgeColor',[0.25,0.25,1.0],...
    'Margin',1,...
    'FontWeight','bold');







end

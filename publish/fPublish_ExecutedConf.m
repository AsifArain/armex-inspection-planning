function [ ppmm ] = fPublish_ExecutedConf( conf_sequence_num,...
                                           sensing_conf,...
                                           cellsScanned,...
                                           SensorRange,...
                                           map_env,...
                                           confFor,...
                                           FoV,...
                                           para_,...
                                           dir_ )
%
% PLOT LAST EXECUTED CONFIGURATION


fprintf('--> Plot conf #%02d.\n',conf_sequence_num);


if     strcmp(confFor,'detection')
    %conf_color = [0,0.65,0];    
    conf_color = [0.0 0.0 0.6];
    
elseif strcmp(confFor,'tomography')
    %conf_color = 'r';
    conf_color = [0.0 0.0 0.6];
    
end


map_filename      = sprintf('%s_map_raycast.dat',para_.environment);
origin_filename   = sprintf('%s_origin_raycast.dat',para_.environment);
cellsize_filename = sprintf('%s_cellsize_raycast.dat',para_.environment);


map_raycast      = load([dir_.env,map_filename]);
origin_raycast   = load([dir_.env,origin_filename]);
cellsize_raycast = load([dir_.env,cellsize_filename]);


% robot_origin_in_map_x = origin_raycast(1);
% robot_origin_in_map_y = origin_raycast(2);
x_origin_cell = origin_raycast(1);
y_origin_cell = origin_raycast(2);

cellsize_recon = para_.ReconTomographyCellSizeM;


switch para_.SensingSystem
        
        %=======================================
    case 'robot-sampling'
        %=======================================
        
        
        %--- robot logs to M for this conf.
        %====================================================
        %{
        fprintf('---> Generating measurement matrix (M) for the last conf.\n');
        measure_file = sprintf('measurements_conf%d.dat',conf_sequence_num);
        
        visualize = 0;
        [ M_m,M_cell ] = fRobotLogs2M( measure_file,...
                                       para_,...
                                       dir_,...
                                       visualize );

        Mfile = sprintf('M_conf%02d.mat',conf_sequence_num);
        save([dir_.MeasurementLogs,Mfile],'M_m','M_cell');
        %}
        
        %--- M matrix for all conf.
        %====================================================
        

        Mfile = sprintf('M_conf%02d.mat',conf_sequence_num);
        M_current = load([dir_.MeasurementLogs,Mfile]);

        M = M_current.M_m;

        %origin_raycast
        %cellsize_raycast
        origin_recon = origin_raycast.*(cellsize_raycast/cellsize_recon);
        
        M(:,1:4)   = M(:,1:4)./cellsize_recon;
        % M(:,1:4)   = M(:,1:4).*(cellsize_raycast/desired_resolution_m);
        %M(1,1:4)
        M(:,[1,3]) = M(:,[1,3])+origin_recon(1);
        M(:,[2,4]) = M(:,[2,4])+origin_recon(2);
        %M(1,1:4)
        
        start_x__all = M(:,1);
        start_y__all = M(:,2);
        end___x__all = M(:,3);
        end___y__all = M(:,4);
        ppmm         = M(:,5);
        
        
        %=======================================
    case 'simulated-sampling'
        %=======================================
                
        % -- retrive integral measurments
        % ----------------------------
        measu_file = sprintf('measurement_conf%02d.dat',conf_sequence_num);
        M = dlmread([dir_.logs,measu_file]);    
        start_x__all = M(:,3);
        start_y__all = M(:,4);
        end___x__all = M(:,5);
        end___y__all = M(:,6);
        ppmm         = M(:,8);
        
end



%**************
%   PLOT
%**************

if para_.Publish_1tARMEx_CoverageArea == 1

    % -- coverage area
    % ----------------------------
    %plot([start_x__all-0.5,end___x__all-0.5],[start_y__all-0.5,end___y__all-0.5],'-r')
    x_pts = [start_x__all(1),(end___x__all)',start_x__all(1)];
    y_pts = [start_y__all(1),(end___y__all)',start_y__all(1)];

    hCoverageArea = plot(x_pts-0.5,y_pts-0.5,'-','LineWidth',1.5,'color',conf_color);
    hPosition = plot(start_x__all(1)-0.5,start_y__all(1)-0.5,...
        'o','color',conf_color,'MarkerFaceColor',conf_color);

end


if para_.Publish_1tARMEx_OpticalBeams == 1
    % -- optical beams
    % ----------------------------
    %beam_colormap = dlmread('beam_colormap.dat');    
    %for i = 1:size(start_x__all,1)        
    %    plot([start_x__all(i),end___x__all(i)],...
    %         [start_y__all(i),end___y__all(i)],...
    %         '-b');
    %end
    % for i = 1:2:size(start_x__all,1)    
    %     color_this_beam = [min([5000,ppmm(i)])/5000,0,0];
    %     plot([start_x__all(i)-0.5,end___x__all(i)-0.5],...
    %         [start_y__all(i)-0.5,end___y__all(i)-0.5],...
    %         '-','color',color_this_beam);
    % end


    beam_colormap = flipud(hot(512));
    max_color_num = 400;
    max_color_con = 5000;
    for i = 1:3:size(start_x__all,1)
        this_colo_ind=round((min([max_color_con,ppmm(i)])/max_color_con)*max_color_num)+1;
        h_beam = plot([start_x__all(i)-0.5,end___x__all(i)-0.5],...
              [start_y__all(i)-0.5,end___y__all(i)-0.5],...
              '-','color',beam_colormap(this_colo_ind,:));
        h_beam.Color(4) = 0.5;
    end
end


%-- mark scanned cells
%-----------------------------------
if para_.Publish_1tARMEx_ScannedArea == 1
    [scanned_r,scanned_c] = ind2sub(size(map_env),cellsScanned);
    hCoverageCells = plot(scanned_r,scanned_c,'xb','MarkerSize',2);
end


%-- mark orientation arrow
%-------------------------------------
if para_.Publish_1tARMEx_OrientationArrow == 1
    arrow_length = 3;
    line_width_arrow = 1; %0.5; %0.25
    arrow_head_size = 30;
    marker_size_arrow = 2; %2
    quiver(sensing_conf(1)-0.5,sensing_conf(2)-0.5,...
           arrow_length*cosd(sensing_conf(3)),...
           arrow_length*sind(sensing_conf(3)),...
           'LineWidth',line_width_arrow,...
           'color',conf_color,...
           'MaxHeadSize',arrow_head_size,... %'Marker','o',...
           'MarkerSize',marker_size_arrow);
end


% -- symbolic field of view
%---------------------------
if para_.Publish_1tARMEx_SymbolicFOV == 1
    
    start_angle  = sensing_conf(3)-FoV/2;
    end___angle  = sensing_conf(3)+FoV/2;
    r            = 2; %SensorRange.cell/6;
    x1           = r*cosd(start_angle);
    y1           = r*sind(start_angle);
    x2           = r*cosd(end___angle);
    y2           = r*sind(end___angle);

    %-- first line
    %-----------------------
    hLine1 = plot([sensing_conf(1)-0.5,x1+sensing_conf(1)-0.5],...
         [sensing_conf(2)-0.5,y1+sensing_conf(2)-0.5],...
          'color',conf_color,'LineWidth',1); %3     
    %-- second line
    %-----------------------
    hLine2 = plot([sensing_conf(1)-0.5,x2+sensing_conf(1)-0.5],...
         [sensing_conf(2)-0.5,y2+sensing_conf(2)-0.5],...
          'color',conf_color,'LineWidth',1); %3
end
  

%-- position
%-----------------------
if para_.Publish_1tARMEx_ConfPosition == 1

    plot(sensing_conf(1)-0.5,sensing_conf(2)-0.5,...
        'o','color',conf_color,'MarkerFacecolor',conf_color,'MarkerSize',2); %8
end

%-- mark conf num
%------------------------------------
if para_.Publish_1tARMEx_ConfSequenceNum == 1
    font_size = 10; %15;
    % text(start_x__all(1)-0.5,start_y__all(1)-0.5,num2str(conf_sequence_num),...
    %         'FontSize',font_size,...
    %         'Color','k',...
    %         'Interpreter','latex',... %        'EdgeColor',conf_color,...
    %         'Margin',1,...
    %         'FontWeight','bold');

    text_dist = 1;
    % sensing_conf(3)
    % sensing_conf(3)-180
    text(start_x__all(1)-0.5 + text_dist*cosd(sensing_conf(3)-180),...
         start_y__all(1)-0.5 + text_dist*sind(sensing_conf(3)-180),...
         num2str(conf_sequence_num),... %['$c_{',num2str(conf_sequence_num),'}$'],...
         'FontSize',font_size,...
         'Color','k',...
         'Interpreter','latex',... %        'EdgeColor',conf_color,...
         'Margin',1,...
         'FontWeight','bold');
end


if strcmp(para_.SensingSystem,'simulation')
    %pause
end    
% pause(3)
% delete(hPosition)

if max(ppmm)<para_.ConcentrationThreshould4TomographyPPM && strcmp(confFor,'detection')
    delete(hCoverageArea);
end

if para_.Publish_1tARMEx_ScannedArea == 1
    delete(hCoverageCells);
end

pause(1)


file_name = sprintf('fig/phd-final-seminar/%s-1t-armex-step-conf%02d',para_.ExperimentTitle,conf_sequence_num);
export_fig([file_name], '-pdf','-r500')
    
end


% function [ ppmm ] = fPlotExecutedConf( conf_sequence_num,...
%                                        sensing_conf,...
%                                        cellsScanned,...
%                                        map_env,...
%                                        confFor,...
%                                        para_,...
%                                        dir_ )
% %
% % PLOT LAST EXECUTED CONFIGURATION
% 
% 
% 
% if     strcmp(confFor,'detection')
%     conf_color = [0,0.65,0];
% elseif strcmp(confFor,'tomography')
%     conf_color = 'r';
% end
% 
% 
% 
% 
% 
% switch para_.SensingSystem
%         
%         %=======================================
%     case 'robot'
%         %=======================================
%         
%         
%         %--- robot logs to M for this conf.
%         %====================================================
%         fprintf('---> Generating measurement matrix (M) for the last conf.\n');
%         measure_file = sprintf('measurement_conf%d.dat',conf_sequence_num);
%         
%         % visualize = 1;
%         [ M_m,M_cell ] = fRobotLogs2M( measure_file,...
%                                        para_,...
%                                        dir_,...
%                                        visualize );
% 
%         Mfile = sprintf('M_conf%02d.mat',conf_sequence_num);
%         save([dir_.MeasurementLogs,Mfile],'M_m','M_cell');
% 
%         %--- M matrix for all conf.
%         %====================================================
% 
% 
%         Mfile = sprintf('M_conf%02d.mat',conf_sequence_num);
%         M_current = load([dir_.MeasurementLogs,Mfile]);
% 
%         M = M_current.M_m;
% 
%         origin_recon = origin_raycast.*(cellsize_raycast/cellsize_recon);
%         
%         M(:,1:4)   = M(:,1:4)./cellsize_recon;
%         % M(:,1:4)   = M(:,1:4).*(cellsize_raycast/desired_resolution_m);
%         M(:,[1,3]) = M(:,[1,3])+origin_recon(1);
%         M(:,[2,4]) = M(:,[2,4])+origin_recon(2);
%         
%         
%         start_x__all = M(:,1);
%         start_y__all = M(:,2);
%         end___x__all = M(:,3);
%         end___y__all = M(:,4);
%         ppmm         = M(:,5);
%         
%         
%         %=======================================
%     case 'simulation'
%         %=======================================
%                 
%         % -- retrive integral measurments
%         % ----------------------------
%         measu_file = sprintf('measurement_conf%02d.dat',conf_sequence_num);
%         M = dlmread([dir_.logs,measu_file]);    
%         start_x__all = M(:,3);
%         start_y__all = M(:,4);
%         end___x__all = M(:,5);
%         end___y__all = M(:,6);
%         ppmm         = M(:,8);
%         
% end
% 
% 
% 
% 
% 
% % -- plot 
% % ----------------------------
% %plot([start_x__all-0.5,end___x__all-0.5],[start_y__all-0.5,end___y__all-0.5],'-r')
% x_pts = [start_x__all(1),(end___x__all)',start_x__all(1)];
% y_pts = [start_y__all(1),(end___y__all)',start_y__all(1)];
% plot(x_pts-0.5,y_pts-0.5,'-','LineWidth',1.5,'color',conf_color);
% plot(start_x__all(1)-0.5,start_y__all(1)-0.5,...
%     'o','color',conf_color,'MarkerFaceColor',conf_color)
% 
% %
% %beam_colormap = dlmread('beam_colormap.dat');    
% %for i = 1:size(start_x__all,1)        
% %    plot([start_x__all(i),end___x__all(i)],...
% %         [start_y__all(i),end___y__all(i)],...
% %         '-b');
% %end
% 
% 
% %-- mark scanned cells
% %-----------------------------------
% % [scanned_r,scanned_c] = ind2sub(size(map_env),cellsScanned);
% % plot(scanned_r,scanned_c,'xb','MarkerSize',2);
% 
% 
% %-- mark orientation arrow
% %-------------------------------------
% arrow_length = 3;
% line_width_arrow = 0.5; %0.25
% arrow_head_size = 20;
% marker_size_arrow = 1; %2
% quiver(sensing_conf(1)-0.5,sensing_conf(2)-0.5,...
%        arrow_length*cosd(sensing_conf(3)),...
%        arrow_length*sind(sensing_conf(3)),...
%        'LineWidth',line_width_arrow,...
%        'color',conf_color,...
%        'MaxHeadSize',arrow_head_size,... %'Marker','o',...
%        'MarkerSize',marker_size_arrow);
% 
% 
% 
% 
% %-- mark conf num
% %------------------------------------
% text(start_x__all(1)-0.5,start_y__all(1)-0.5,num2str(conf_sequence_num),...        
%         'FontSize',15,...
%         'Color','k',...
%         'Interpreter','latex',... %        'EdgeColor',conf_color,...
%         'Margin',1,...
%         'FontWeight','bold');
%     
% 
% end

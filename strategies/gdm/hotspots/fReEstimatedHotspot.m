function [ hotspot_estimated,...
           cell_coords_x,...
           cell_coords_y,...
           mean_map ] = fReEstimatedHotspot( conf_num_hot,...
                                             hotspot_current,...                                             
                                             executedConf_num_thisHot,...
                                             map_env,...
                                             para_,...
                                             dir_,...
                                             visualize )

% fReEstimatedHotspot re-estimate hotspot based on the new measurements.
%
% Author:  Asif Arain
% Project: SPP-REM (JFR-2016)

% ____________________________________________________________________ %
%                     integral measurements                            %
% ____________________________________________________________________ %

% -- measurement matrix
% 1 conf_num,...              % configuration number
% 2 sensor___timestamp,...    % time stamp
% 3 beam_initial_x,...        % beam start position along x-axis
% 4 beam_initial_y,...        % beam start position along y-axis
% 5 beam_voxel_end_xy(1),...  % beam final position along x-axis
% 6 beam_voxel_end_xy(2),...  % beam final position along y-axis
% 7 beam____length,...        % beam length
% 8 ppm_m,...                 % integral concentrtion value (ppm.m)
% 9 weight

% dir_.env
% visualize = 0

%**********************************************************************
%
%                           MAP:
%
%**********************************************************************

cellsize_recon = para_.ReconTomographyCellSizeM;

file__map_env      = sprintf('%s_map_coverage.dat',para_.environment);
file__origin_env   = sprintf('%s_origin_coverage.dat',para_.environment);
file__cellsize_env = sprintf('%s_cellsize_coverage.dat',para_.environment);

map_env      = load([dir_.env,file__map_env]);
origin_env   = load([dir_.env,file__origin_env]);
cellsize_env = load([dir_.env,file__cellsize_env]);


file__map_raycast      = sprintf('%s_map_raycast.dat',para_.environment);
file__origin_raycast   = sprintf('%s_origin_raycast.dat',para_.environment);
file__cellsize_raycast = sprintf('%s_cellsize_raycast.dat',para_.environment);

map_raycast      = load([dir_.env,file__map_raycast]);
origin_raycast   = load([dir_.env,file__origin_raycast]);
cellsize_raycast = load([dir_.env,file__cellsize_raycast]);


%**********************************************************************
%
%                   INTEGRAL MEASUREMENTS
%
%**********************************************************************
switch para_.SensingSystem
        
        %=======================================
    case 'robot'
        %=======================================
        
        %--- robot logs to M for this conf.
        %====================================================
        fprintf('---> Generating measurement matrix (M) for the last conf.\n');
        measure_file = sprintf('measurements_conf%d.dat',conf_sequence_num);
        % visualize = 1;
        [ M_m,M_cell ] = fRobotLogs2M( measure_file,...
                                       para_,...
                                       dir_,...
                                       visualize );

        Mfile = sprintf('M_conf%02d.mat',conf_sequence_num);
        save([dir_.MeasurementLogs,Mfile],'M_m','M_cell');

        
        
        %--- M matrix for all conf.
        %====================================================

        % -- measurement matrix
        fprintf('---> Combining all measurements (M).\n');

        M = [];
        for num = 1:conf_sequence_num

            Mfile = sprintf('M_conf%02d.mat',num);
            M_current = load([dir_.MeasurementLogs,Mfile]);

            M = [M; M_current.M_m];
            %M = [M; M_current.M_cell];
        end

        %{
        origin_recon = origin_raycast.*(cellsize_raycast/cellsize_recon);
        M(:,1:4)   = M(:,1:4)./cellsize_recon;
        % M(:,1:4)   = M(:,1:4).*(cellsize_raycast/desired_resolution_m);
        M(:,[1,3]) = M(:,[1,3])+origin_recon(1);
        M(:,[2,4]) = M(:,[2,4])+origin_recon(2);
        %}
        
        M(:,1:4)   = M(:,1:4)./cellsize_env;
        M(:,[1,3]) = M(:,[1,3])+origin_env(1);
        M(:,[2,4]) = M(:,[2,4])+origin_env(2);

        if visualize == 1

            %--- environment map
            figure('name','Optical beams'); 
            imshow(map_env','InitialMagnification',800); hold on;
            set(gca,'YDir','normal'); colormap('gray'); caxis([-2 0.5])
            plot(origin_env(1),origin_env(2),'ob');

            %{
            %--- optical beams (with color proportional to concentrations)
            beam_colormap = flipud(bone(512));
            beam_weight = M(:,5)/max(M(:,5));
            for i = 1:size(M,1)        
                plot([M(i,1),M(i,3)],[M(i,2),M(i,4)],...
                    '-','Color',beam_colormap(round(255*beam_weight(i))+100,:));
            end
            %--- coverage area
            x_pts = [M(1,1);M(:,3);M(1,1)];
            y_pts = [M(1,2);M(:,4);M(1,2)];
            plot(x_pts-0.0,y_pts-0.0,'-','LineWidth',1.5,'color','b');
            %}


            M_env = M(:,1:4).*(cellsize_recon/cellsize_env);

            %--- optical beams (with color proportional to concentrations)
            beam_colormap = flipud(bone(512));
            beam_weight = M(:,5)/max(M(:,5));
            for i = 1:size(M_env,1)
                plot([M_env(i,1),M_env(i,3)],[M_env(i,2),M_env(i,4)],...
                    '-','Color',beam_colormap(round(255*beam_weight(i))+100,:));
            end
            %--- coverage area
            x_pts = [M_env(1,1);M_env(:,3);M_env(1,1)];
            y_pts = [M_env(1,2);M_env(:,4);M_env(1,2)];
            plot(x_pts-0.0,y_pts-0.0,'-','LineWidth',1.5,'color','b');

        end


        

        %=======================================
    case 'simulation'
        %=======================================
        
        
        M      = [];
        beam_x = [];
        beam_y = [];
        for i = 1:numel(executedConf_num_thisHot)%conf_num_hot-1%n_c

            measu_file = sprintf('measurement_conf%02d.dat',...
                executedConf_num_thisHot(i));

            M_current = dlmread([dir_.logs,measu_file]);

            %M_current = dlmread([dir_.meas_capsule,measu_file{i}]);
            % - update list of all beams
            for j = 1:size(M_current,1)
                beam_x = [beam_x,M_current(j,3),M_current(j,5),NaN];
                beam_y = [beam_y,M_current(j,4),M_current(j,6),NaN];
            end
            M = [M;M_current];
        end


        % -- splite M matrix
        % {
        start_x__all = M(:,3);
        start_y__all = M(:,4);
        end___x__all = M(:,5);
        end___y__all = M(:,6);
        %M(:,8)
        weight___all = M(:,8)/max(M(:,8));
        %}
        
        %-- final M matrix (new)
        M = M(:,[3,4,5,6,8]);
end










% CELL_SIZE = 1;

% -- line model for gdm
fprintf('- Line model for 2d gdm.\n');

recon_size = 'measurement-beams';
% [ A,...
%   ppmm,...
%   cell_coords_x,...
%   cell_coords_y,...
%   num_cells_x,...
%   num_cells_y ] = fLineModel( M,...
%                               map_env,...
%                               recon_size,...
%                               cellsize_recon,...
%                               para_,...
%                               dir_ );
[ A,...
  ppmm,...
  cell_coords_x,...
  cell_coords_y,...
  num_cells_x,...
  num_cells_y ] = fLineModelJFR( M,...
                                 recon_size,...
                                 cellsize_recon,...
                                 para_,...
                                 dir_ );
A = full(A);
% dir_.Recon

%--- save cell_coords_x
x_coord_file = sprintf('x_coord_conf%02d.dat',conf_sequence_num);
fileID = fopen([dir_.Recon,x_coord_file],'wt'); fclose(fileID);
dlmwrite([dir_.Recon,x_coord_file],cell_coords_x,'-append','delimiter',' ');
fileID = fopen([dir_.Recon,x_coord_file],'a'); fclose(fileID);

%--- save cell_coords_x
x_coord_file = sprintf('x_coord.dat');
fileID = fopen([dir_.Recon,x_coord_file],'wt'); fclose(fileID);
dlmwrite([dir_.Recon,x_coord_file],cell_coords_x,'-append','delimiter',' ');
fileID = fopen([dir_.Recon,x_coord_file],'a'); fclose(fileID);

%--- save cell_coords_y
y_coord_file = sprintf('y_coord_conf%02d.dat',conf_sequence_num);
fileID = fopen([dir_.Recon,y_coord_file],'wt'); fclose(fileID);
dlmwrite([dir_.Recon,y_coord_file],cell_coords_y,'-append','delimiter',' ');
fileID = fopen([dir_.Recon,y_coord_file],'a'); fclose(fileID);

%--- save cell_coords_y
y_coord_file = sprintf('y_coord.dat');
fileID = fopen([dir_.Recon,y_coord_file],'wt'); fclose(fileID);
dlmwrite([dir_.Recon,y_coord_file],cell_coords_y,'-append','delimiter',' ');
fileID = fopen([dir_.Recon,y_coord_file],'a'); fclose(fileID);
% ____________________________________________________________________ %
%                2d gas distribution map (line model)                  %
% ____________________________________________________________________ %

fprintf('- 2d gas distribution map.\n');
[ mean_map ] = fGDM( A,...
                     ppmm,...
                     cell_coords_x,...
                     cell_coords_y,...
                     num_cells_x,...
                     num_cells_y,...
                     visualize );
                
                
%999999999999999999999999999999999999999999999999999999999999999999999999999
map_recon = mean_map';

%--- save mean map
recon_file = sprintf('reconstruction_conf%02d.dat',conf_sequence_list_hc(end));
fileID = fopen([dir_.Recon,recon_file],'wt'); fclose(fileID);
dlmwrite([dir_.Recon,recon_file],map_recon,'-append','delimiter',' ');
fileID = fopen([dir_.Recon,recon_file],'a'); fclose(fileID);

%--- save mean map all
recon_file = sprintf('reconstruction.dat');
fileID = fopen([dir_.Recon,recon_file],'wt'); fclose(fileID);
dlmwrite([dir_.Recon,recon_file],map_recon,'-append','delimiter',' ');
fileID = fopen([dir_.Recon,recon_file],'a'); fclose(fileID);

map_recon_col = mean_map;
gdm_colormap = flipud(autumn(512)); % for main fig. 
max_concentration = max(map_recon_col(:)); % main fig.
delta = max_concentration/(size(gdm_colormap,1)-1);
linear_color_number = round(map_recon_col(:)./delta)+1;

colR = round(gdm_colormap(linear_color_number,1)*255);
colG = round(gdm_colormap(linear_color_number,2)*255);
colB = round(gdm_colormap(linear_color_number,3)*255);

thisfile = sprintf('reconstructionColorR.dat');
fileID = fopen([dir_.Recon,thisfile],'wt'); fclose(fileID);
dlmwrite([dir_.Recon,thisfile],colR,'-append','delimiter',' ');
fileID = fopen([dir_.Recon,thisfile],'a'); fclose(fileID);

thisfile = sprintf('reconstructionColorR_conf%02d.dat',conf_sequence_list_hc(end));
fileID = fopen([dir_.Recon,thisfile],'wt'); fclose(fileID);
dlmwrite([dir_.Recon,thisfile],colR,'-append','delimiter',' ');
fileID = fopen([dir_.Recon,thisfile],'a'); fclose(fileID);

thisfile = sprintf('reconstructionColorG.dat');
fileID = fopen([dir_.Recon,thisfile],'wt'); fclose(fileID);
dlmwrite([dir_.Recon,thisfile],colG,'-append','delimiter',' ');
fileID = fopen([dir_.Recon,thisfile],'a'); fclose(fileID);

thisfile = sprintf('reconstructionColorG_conf%02d.dat',conf_sequence_list_hc(end));
fileID = fopen([dir_.Recon,thisfile],'wt'); fclose(fileID);
dlmwrite([dir_.Recon,thisfile],colG,'-append','delimiter',' ');
fileID = fopen([dir_.Recon,thisfile],'a'); fclose(fileID);

thisfile = sprintf('reconstructionColorB.dat');
fileID = fopen([dir_.Recon,thisfile],'wt'); fclose(fileID);
dlmwrite([dir_.Recon,thisfile],colB,'-append','delimiter',' ');
fileID = fopen([dir_.Recon,thisfile],'a'); fclose(fileID);

thisfile = sprintf('reconstructionColorB_conf%02d.dat',conf_sequence_list_hc(end));
fileID = fopen([dir_.Recon,thisfile],'wt'); fclose(fileID);
dlmwrite([dir_.Recon,thisfile],colB,'-append','delimiter',' ');
fileID = fopen([dir_.Recon,thisfile],'a'); fclose(fileID);


%-- copy to the ROS folder for publication
filesFrom = dir_.Recon;
filesTo = sprintf('%stwo-step-exploration-tomography-adaptive/',dir_.ROSLogs);

%-- copy command
syncCommand = sprintf('rsync -a --ignore-existing %s %s',filesFrom,filesTo);

%-- copy log files
system(syncCommand);
%999999999999999999999999999999999999999999999999999999999999999999999999999




% mean_map
if sum(mean_map(:)) == 0
    mean_map = 0.1*rand(size(mean_map));
end
W_recon = mean_map./(max(mean_map(:)));

[map_row,map_col] = meshgrid(1:size(mean_map,1),1:size(mean_map,2));


% max(mean_map(:))

% W_recon(:)

map_row_Wmean = sum(sum(W_recon.*map_row'))/sum(W_recon(:));
map_col_Wmean = sum(sum(W_recon.*map_col'))/sum(W_recon(:));

% -- matrix operations to find weighted mean
%{
mean_row = ([1:size(W,2)]*(sum(W,1))')/sum(W(:)) % or
mean_row = sum(sum(W.*row'))/sum(W(:))

mean_col = ([1:size(W,2)]*(sum(W,1))')/sum(W(:)) % or
mean_col = sum(sum(W.*col'))/sum(W(:))
%}

% plot(map_col_Wmean,map_row_Wmean,'sg')
% plot(map_row_Wmean,map_col_Wmean,'sg')

% round(map_row_Wmean)
map_row_Wmean_absolute = cell_coords_y(round(map_row_Wmean));
map_col_Wmean_absolute = cell_coords_x(round(map_col_Wmean));

if visualize == 1
plot(map_col_Wmean_absolute,map_row_Wmean_absolute,'sr')
end

% hotspot_newborn = [map_row_Wmean_absolute,map_col_Wmean_absolute]
hotspot_estimated = [map_col_Wmean_absolute,map_row_Wmean_absolute];


% -- estimation for single conf
if (conf_num_hot-1) == 1
    
    %
    % -- configuration placements
    conf_positions = unique([start_x__all,start_y__all],'rows');
    conf_________x = conf_positions(:,1);
    conf_________y = conf_positions(:,2);
    
    dist2prevHotspot = pdist2([conf_________x,conf_________y],hotspot_current,...
                              'euclidean');
     
    ang2estimatedHotspot = rad2deg(angle2Points([conf_________x,conf_________y],...
                                     hotspot_estimated));
    
    hotspot_estimated_x = conf_________x + dist2prevHotspot*cosd(ang2estimatedHotspot);
    hotspot_estimated_y = conf_________y + dist2prevHotspot*sind(ang2estimatedHotspot);
    
    hotspot_estimated = [hotspot_estimated_x,hotspot_estimated_y];
    
end

% ____________________________________________________________________ %
%                             visualization                            %
% ____________________________________________________________________ %

if visualize == 1
    
    %figure; 
    hold on;
    
    % -- all beams
    %plot(beam_x,beam_y,'-','color',[0.9,0.9,0.9])
    
    % -- save beam color map
    %{
    %beam_colormap = colormap(flipud(colormap(autumn(256))));
    %beam_colormap = colormap(flipud(colormap(copper(512))));
    %beam_colormap = colormap(cool(256));
    beam_colormap = colormap(flipud(colormap(bone(512))));
    fileID = fopen('beam_colormap.dat','a'); fclose(fileID);
    dlmwrite('beam_colormap.dat',beam_colormap,'precision','%.4f');
    fileID = fopen('beam_colormap.dat','a'); fclose(fileID);
    %}
    beam_colormap = dlmread('beam_colormap.dat');
    
    for i = 1:size(start_x__all,1)        
        plot([start_x__all(i),end___x__all(i)],...
             [start_y__all(i),end___y__all(i)],...
             '-','Color',beam_colormap(round(255*weight___all(i))+100,:));
    end
        
    %plot(Xpts(:,1),Xpts(:,2),'xr')
    
    % -- high concentration beams
    %plot(beam_x_con,beam_y_con,'-','color',[1.0,0.1,0.1])
    
    % -- intersection points
    %plot(beam_intersections(1,:),beam_intersections(2,:),'og')
    
    
    % -- clusters after identified
    %{
    colors = lines(num_of_clusters);
    for i = 1:num_of_clusters
        %{
        plot(clusters{i}(:,1),clusters{i}(:,2),...
            'o','color',colors(i,:),'MarkerFaceColor',colors(i,:))
        plot(clusters{i}(:,1),clusters{i}(:,2),...
            'o','color',colors(i,:),'MarkerFaceColor',colors(i,:))
        %}
        
        plot(Xpts(Xclusters==i,1),...
             Xpts(Xclusters==i,2),...
             'o','color',colors(i,:),'MarkerFaceColor',colors(i,:),...
             'MarkerSize',1.65);
    end
    %}

    %axis equal
    
    % -- hotspot of all intersection points
    %%plot(hotspot_intersections(1),hotspot_intersections(2),'or')

    % -- hotspots of clusters
    %plot(hotspot_clusters(:,1),hotspot_clusters(:,2),...
    %    'dk','MarkerSize',10,'MarkerFaceColor','g')
    
    % -- hotspot previous
    %plot(hotspot_current(1),hotspot_current(2),...
    %    'ok','MarkerSize',10,'MarkerFaceColor','y')
    
    % -- hotspot final
    %plot(hotspot_newborn(1),hotspot_newborn(2),...
    %    'ok','MarkerSize',10,'MarkerFaceColor','g')
    
end

% pause

end
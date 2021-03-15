function [ hc_estimated,...
           cell_coords_x,...
           cell_coords_y,...
           map_recon_local ] = fReEstimatedHcs1t( conf_sequence_list_hc,...
                                                  SensorRange,...
                                                  OBSenv,...
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


tComput_i = tic; % initialize pre computational time.

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
    case 'robot-sampling'
        %=======================================
        
        %--- robot logs to M for this conf.
        %====================================================
        %{
        %fprintf('---> Generating measurement matrix (M) for the last conf.\n');
        fprintf('---> Generating M matrix for conf# %02d.\n',conf_sequence_list_hc(end));
        
        measure_file = sprintf('measurements_conf%d.dat',conf_sequence_list_hc(end));
        % visualize = 0;
        [ M_m,M_cell ] = fRobotLogs2M( measure_file,...
                                       para_,...
                                       dir_,...
                                       visualize );

        Mfile = sprintf('M_conf%02d.mat',conf_sequence_list_hc(end));
        save([dir_.MeasurementLogs,Mfile],'M_m','M_cell');
        %}
        
        %--- M matrix for all conf.
        %====================================================
        
        % -- measurement matrix
        %fprintf('---> Combining all measurements (M).\n');
        disp(['---> Combining measurements (M) of confs: ',...
            num2str(conf_sequence_list_hc)]);
        
        M = [];
        %for num = 1:conf_sequence_list_hc
        for num = conf_sequence_list_hc
            
            Mfile = sprintf('M_conf%02d.mat',num);
            M_current = load([dir_.logs,Mfile]);
            
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
        
        %origin_recon = origin_raycast.*(cellsize_raycast/cellsize_recon);        
        M(:,1:4)   = M(:,1:4)./cellsize_env;
        M(:,[1,3]) = M(:,[1,3])+origin_env(1);
        M(:,[2,4]) = M(:,[2,4])+origin_env(2);
        
        %visualize = 1
        
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
        
        %pause
        %{
        M = [];
        for i = conf_sequence_list_hc
            measu_file = sprintf('measurements_conf%d.dat',i);
            [ M_m,M_cell ] = fRobotLogs2M( measure_file,...
                                           para_,...
                                           dir_,...
                                           visualize );
            [dir_.logs,measu_file]
            M_current = dlmread([dir_.logs,measu_file])
            pause
            M = [M;M_current];
        end
        
        cellsize_raycast
        cellsize_recon
        
        origin_recon = origin_raycast.*(cellsize_raycast/cellsize_recon)
        M
        M(:,1:4)   = M(:,1:4)./cellsize_recon;
        % M(:,1:4)   = M(:,1:4).*(cellsize_raycast/desired_resolution_m);
        M(:,[1,3]) = M(:,[1,3])+origin_recon(1);
        M(:,[2,4]) = M(:,[2,4])+origin_recon(2);
        %}
        

        %=======================================
    case 'simulated-sampling'
        %=======================================
        
        disp(['---> Combining measurements (M) of confs: ',...
            num2str(conf_sequence_list_hc)]);
        
        M      = [];
        beam_x = [];
        beam_y = [];
        for i = conf_sequence_list_hc

            measu_file = sprintf('measurement_conf%02d.dat',i);
            M_current = dlmread([dir_.logs,measu_file]);

            %-- DEBUG
            %h_confmark(i) = plot(M_current(1,3),M_current(1,4),'ob','MarkerSize',8);
            % - update list of all beams
            for j = 1:size(M_current,1)
                beam_x = [beam_x,M_current(j,3),M_current(j,5),NaN];
                beam_y = [beam_y,M_current(j,4),M_current(j,6),NaN];
            end
            M = [M;M_current];
        end

        %-- DEBUG
        % disp('Confs used to estimate the hotspot are marked with blue circle');
        % pause
        % delete(h_confmark)
        % pause

        % -- splite M matrix
        %{
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



%**********************************************************************
%
%                    (LOCAL) RECONSTRUCTION MAP
%
%**********************************************************************

% CELL_SIZE = 1;

% -- line model for gdm
%fprintf('- Line model for 2d gdm.\n');
fprintf('---> Line model.\n');

recon_size = 'measurement-beams';

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
x_coord_file = sprintf('x_coord_conf%02d.dat',conf_sequence_list_hc(end));
fileID = fopen([dir_.Recon,x_coord_file],'wt'); fclose(fileID);
dlmwrite([dir_.Recon,x_coord_file],cell_coords_x,'-append','delimiter',' ');
fileID = fopen([dir_.Recon,x_coord_file],'a'); fclose(fileID);

%--- save cell_coords_x
x_coord_file = sprintf('x_coord.dat');
fileID = fopen([dir_.Recon,x_coord_file],'wt'); fclose(fileID);
dlmwrite([dir_.Recon,x_coord_file],cell_coords_x,'-append','delimiter',' ');
fileID = fopen([dir_.Recon,x_coord_file],'a'); fclose(fileID);

%--- save cell_coords_y
y_coord_file = sprintf('y_coord_conf%02d.dat',conf_sequence_list_hc(end));
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

fprintf('---> Local reconstruction.\n');

[ mean_map ] = fGDM( A,...
                     ppmm,...
                     cell_coords_x,...
                     cell_coords_y,...
                     num_cells_x,...
                     num_cells_y,...
                     visualize );
                 


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
switch para_.SensingSystem
    case 'robot-sampling'
        filesFrom = dir_.Recon;
        filesTo = sprintf('%sone-step-exploration/',dir_.ROSLogs);

        %-- copy command
        syncCommand = sprintf('rsync -a %s %s',filesFrom,filesTo);

        %-- copy log files
        system(syncCommand);
        
        fprintf('---- The reconstruction is synchronized to display.\n');

        
    case 'simulated-sampling'
        %        
end



% if visualize ==1
%     % ----------------    
%     % -- splite M matrix
%     start_x__all = M(:,3);
%     start_y__all = M(:,4);
%     end___x__all = M(:,5);
%     end___y__all = M(:,6);
%     plot([start_x__all';end___x__all'],[start_y__all';end___y__all'],'-g')
%     %plot([start_y__all';end___y__all'],[start_x__all';end___x__all'],'-g')
%     axis equal
%     
%     %figure; imshow(mean_map','InitialMagnification',800); hold on;
%     %set(gca,'YDir','normal')
%     %colormap('hot'); % not go back to gray scale.
%     %%plot([start_x__all';end___x__all'],[start_y__all';end___y__all'],'-g')
%     %plot([start_y__all';end___y__all'],[start_x__all';end___x__all'],'-g')    
%     %pause(1)
%     
% end

% if visualize ==1
%     
%   
%     % - ground truth
%     map_conc_______gt = dlmread([dir_.maps,gt_filename]);
%     
%     figure; hold on;
%     image(map_conc_______gt,'CDataMapping','scaled')
%     colormap hot; 
%     %colormap gray; 
%     colormap(flipud(colormap))
%     colorbar;
%     
%     % ----------------    
%     % -- splite M matrix
%     start_x__all = M(:,3);
%     start_y__all = M(:,4);
%     end___x__all = M(:,5);
%     end___y__all = M(:,6);
%     plot([start_x__all';end___x__all'],[start_y__all';end___y__all'],'-g')
%     %plot([start_y__all';end___y__all'],[start_x__all';end___x__all'],'-g')
%     axis equal
%     
%     conf_x = unique(start_x__all)
%     conf_y = unique(start_y__all)
% end




%---- in case reconstruction map is empty
if sum(mean_map(:)) == 0
    mean_map = 0.1*rand(size(mean_map));
end

map_recon_local = mean_map';

%**********************************************************************
%
%                   HOTSPOT ESTIMATION
%
%**********************************************************************

switch para_.HotspotReEstimationMethod
    
    case 'ClusteringBased'
        
        fprintf('---> Cluster based estimation of hcs.\n');
        
        %visualize = 1
        [ hc_estimated ] = fClustering4Hotspots( map_env,...
                                                 map_recon_local,...
                                                 cell_coords_x,...
                                                 cell_coords_y,...
                                                 SensorRange,...
                                                 OBSenv,...
                                                 cellsize_env,...
                                                 para_,...
                                                 dir_,...
                                                 visualize );
    case 'SimpleWeightedMean'
        
        fprintf('---> Simple weighted mean based estimation of hcs.\n');
        
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

        %visualize = 1
        if visualize == 1
        plot(map_col_Wmean_absolute,map_row_Wmean_absolute,'sr')
        end

        % hotspot_newborn = [map_row_Wmean_absolute,map_col_Wmean_absolute]
        hc_estimated = [map_col_Wmean_absolute,map_row_Wmean_absolute];
        
end





% -- estimation for single conf
% if (conf_num_hot-1) == 1
%     
%     %
%     % -- configuration placements
%     conf_positions = unique([start_x__all,start_y__all],'rows');
%     conf_________x = conf_positions(:,1);
%     conf_________y = conf_positions(:,2);
%     
%     dist2prevHotspot = pdist2([conf_________x,conf_________y],hotspot_current,...
%                               'euclidean');
%      
%     ang2estimatedHotspot = rad2deg(angle2Points([conf_________x,conf_________y],...
%                                      hotspot_estimated));
%     
%     hotspot_estimated_x = conf_________x + dist2prevHotspot*cosd(ang2estimatedHotspot);
%     hotspot_estimated_y = conf_________y + dist2prevHotspot*sind(ang2estimatedHotspot);
%     
%     hotspot_estimated = [hotspot_estimated_x,hotspot_estimated_y];
%     
% end

% ____________________________________________________________________ %
%                             visualization                            %
% ____________________________________________________________________ %
% visualize = 1
% if visualize == 1
%     
%     
%     start_x__all = M(:,1);
%     start_y__all = M(:,2);
%     end___x__all = M(:,3);
%     end___y__all = M(:,4);
%     weight___all = M(:,5)/max(M(:,5));
%         
%     %figure; 
%     hold on;
%         
%     beam_colormap = flipud(bone(512));
%     for i = 1:size(start_x__all,1)        
%         plot([start_x__all(i),end___x__all(i)],...
%              [start_y__all(i),end___y__all(i)],...
%              '-','Color',beam_colormap(round(255*weight___all(i))+100,:));
%     end 
% end
%}
        


%**********************************************************************
% Total computation time.
%**********************************************************************
tComput = 1e-4*(round(toc(tComput_i)*1e4));
fprintf(1,'---- Computation time: %0.4f sec \n',tComput);


end


function [ hc_estimated ] = fClustering4Hotspots( map_env,...
                                                  map_recon_local,...
                                                  cell_coords_x,...
                                                  cell_coords_y,...
                                                  SensorRange,...
                                                  OBSenv,...
                                                  cellsize_env,...
                                                  para_,...
                                                  dir_,...
                                                  visualize )
%
%


% map_hotspot = zeros(size(map_recon_local));

map_hotspot = zeros(size(map_env));
for i = 1:numel(cell_coords_x)
    for j = 1:numel(cell_coords_y)
        map_hotspot(cell_coords_x(i),cell_coords_y(j)) = map_recon_local(i,j);
    end
end
       


% ind = find(gd_map_hotspot(:))>500;
% gd_map_hotspot(ind) = 1;
% map_hotspot((map_recon_local(:)>para_.highConcentrationThreshold_PPM)) = 1;
% map_hotspot((map_hotspot(:)>para_.highConcentrationThreshold_PPM)) = 1;
% map_hotspot((map_hotspot(:)>500)) = 1;


map_positive = zeros(size(map_env));
map_positive((map_hotspot(:)>para_.highConcentrationThreshold_PPM)) = 1;

% test_map_resized2_obs = find(test_map_resized2(:)==0);
% gd_map_coverage(test_map_resized2_obs) = 0;


% -- hotspots based on clustering

% - high concentration cells
% [high_r,high_c] = find(map_hotspot);
[high_r,high_c] = find(map_positive);
Hcells = [high_r,high_c];



% plot
% visualize = 1
if visualize==1
figure('name','Positive concentrations'); 
imshow(map_env','InitialMagnification',350); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
color_positive_con = [255,140,0]/255;
pause(1)
for i = 1:size(Hcells,1)
    plot(Hcells(i,1)+0.5,Hcells(i,2)+0.5,...
         's',...
         'color',color_positive_con,...
         'MarkerFaceColor',color_positive_con,...
         'MarkerSize',2.75);
end
pause(2)
filename = sprintf('%s%s-hotspots-positive-concentrations',...
    dir_.hotspots,para_.ExperimentTitle);
export_fig(filename, '-pdf','-r500')

pause(1)
end


%*******************************************************************
% 
%               CLUSTERING
% 
%*******************************************************************

%-- linkage cutoff for clusters
cut_off = 2.5; %3.5; %median([Z(end-2,3) Z(end-1,3)])


%--- filter outliers
%---------------------------------------------
% if size(Hcells,1)>1
%     
%     
%     %-- cluster linkage tree
%     Z = linkage(Hcells,'single');
%     %figure; dendrogram(Z,'ColorThreshold',cut_off)
% 
%     % - high concentration clusters
%     clusters2Filter = cluster(Z,'cutoff',cut_off,'criterion','distance');
%         
%     % - number of clusters
%     num_of_clusters = numel(unique(clusters2Filter));
%     
%     %-- weighted map
%     Wmap = map_hotspot./sum(map_hotspot(:)); % total weight is 1    
%     
%     %-- filter
%     filter_ind = [];
%     for i = 1:num_of_clusters
%         
%         %-- this cluster
%         thisCluster = Hcells(clusters2Filter==i,:);
%         
%         %-- if number of observation are less than a number
%         if size(thisCluster,1) < 5
%             
%             %-- concentration weight of this cluster
%             thisWeight = 0;
%             for j = 1:size(thisCluster,1)
%                 thisWeight = thisWeight + Wmap(thisCluster(j,1),thisCluster(j,2));
%             end
%             %thisWeight
%             %-- if concentration weight of this cluster is less than a
%             %number
%             if thisWeight < 0.01
%                 %-- filter out the observations (cells)
%                 %Hcells_filtered(clusters2Filter==i,:) = [];
%                 
%                 filter_ind = [filter_ind;find(clusters2Filter==i)];
%                 
%             end
%         end
%     end
%     
%     Hcells(filter_ind,:) = [];
%     
%     
%     
%     % plot
%     if visualize==1
%         figure('name','Filtered observations (cells)'); 
%         imshow(map_env','InitialMagnification',350); hold on;
%         set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
%         pause(1)
%         for i = 1:size(Hcells,1)
%             plot(Hcells(i,1),Hcells(i,2),...
%                  's',...
%                  'color','k',...%color_positive_con,...
%                  'MarkerFaceColor',color_positive_con,...
%                  'MarkerSize',3.75); %2.75
%         end
%         pause(2)
%         filename = 'hotspots-filtered-cells';
%         % print('-painters','-dpdf','-r500',filename); % pdf    
%         export_fig([dir_.TomographyMaps,filename], '-pdf','-r500')
%         pause(1)
%     end
%         
% end


if size(Hcells,1)>1

    %*******************************************
    %               FIRST LAYER
    %*******************************************
    
    %*******************************************
    % Step-2: CLUSTERING
    % clusters based on distance
    %*******************************************


    % - cluster linkage tree
    Z = linkage(Hcells,'single');
    % figure; dendrogram(Z);

    % - high concentration clusters
    Hclusters1 = cluster(Z,'cutoff',cut_off,'criterion','distance');

    % Hclusters = cluster(Z,'cutoff',1);
    % Hclusters = cluster(Z,'cutoff',2.5,'Depth',2);
    % Hclusters = cluster(Z,'cutoff',1,'MaxClust',5);

    % - number of clusters
    num_of_clusters1 = numel(unique(Hclusters1));



    % plot
    % visualize = 1
    if visualize==1
        figure('name','clusters - 1st layer'); 
        imshow(map_env','InitialMagnification',350); hold on;
        set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
        pause(1)
        % colors = lines(num_of_clusters1);
        %colors = lines(9);
        % colors = linspecer(20);
        % colors = linspecer(30);
        colors = linspecer(num_of_clusters1);
        for i = 1:num_of_clusters1    
            plot(Hcells(Hclusters1==i,1)+0.5,Hcells(Hclusters1==i,2)+0.5,...
                 's',...
                 'color',colors(i,:),...
                 'MarkerFaceColor',colors(i,:),...
                 'MarkerSize',2.75);
        end
        pause(2)
        filename = sprintf('%s%s-hotspots-clusters1-agglomerative',...
            dir_.hotspots,para_.ExperimentTitle);
        %print('-painters','-dpdf','-r500',filename); % pdf
        export_fig(filename, '-pdf','-r500')
        pause(1)
    end

    %*******************************************
    %           SECOND LAYER
    %*******************************************
    
    %*******************************************
    % Step-3: CLUSTERING
    % sub-divide lengthy clusters
    %*******************************************

    hot_clusters = zeros(size(Hcells,1),1);
    clust_num = 1;

    for i = 1:num_of_clusters1

        cluster1_this = [Hcells(Hclusters1==i,1),Hcells(Hclusters1==i,2)];

        dist_this = zeros(size(cluster1_this,1));
        for j = 1:size(dist_this,1)
            dist_this(j,:) = pdist2(cluster1_this(j,:),cluster1_this,'euclidean');
        end

        num_of_desired_clusters2_this = ...
            round(max(dist_this(:))/(SensorRange.cell/2));
        if num_of_desired_clusters2_this<1
            num_of_desired_clusters2_this = 1;
        end


        if num_of_desired_clusters2_this>1

            %if num_of_clusters2_this<1
            %    num_of_clusters2_this = 1;
            %end
                   
            distance_type = 'cityblock';
            %distance_type = 'cosine';
            %distance_type = 'correlation';
            %distance_type = 'hamming'; % beykar            
            [Hclusters2_this,~] = kmeans( cluster1_this,...
                                          num_of_desired_clusters2_this,...
                                          'Distance',distance_type );

            % - number of clusters
            num_of_desired_clusters2_this = numel(unique(Hclusters2_this));
            hot_clusters(Hclusters1==i) = Hclusters2_this+clust_num-1;
            clust_num = clust_num+num_of_desired_clusters2_this;

        else

            hot_clusters(Hclusters1==i) = clust_num;
            clust_num = clust_num+1;
        end
    end


    % - number of clusters
    num_of_clusters = numel(unique(hot_clusters));



    % plot
    if visualize==1
    figure('name','clusters - 2nd layer'); 
    imshow(map_env','InitialMagnification',350); hold on;
    set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
    pause(1)
    % colors = lines(num_of_clusters);
    colors = linspecer(num_of_clusters+num_of_clusters1);
    for i = 1:num_of_clusters
        plot(Hcells(hot_clusters==i,1)+0.5,Hcells(hot_clusters==i,2)+0.5,...
             's',...
             'color',colors(num_of_clusters1+i,:),...
             'MarkerFaceColor',colors(num_of_clusters1+i,:),...
             'MarkerSize',2.75);
    end
    pause(2)
    filename = sprintf('%s%s-hotspots-clusters2-kmeans',...
        dir_.hotspots,para_.ExperimentTitle);
    %print('-painters','-dpdf','-r500',filename); % pdf
    export_fig(filename, '-pdf','-r500');

    pause(1)

    end 

    
    %************************************************
    % Step-4: ESTIMATE HOTSPOT CENTERS
    % hotspot centers are weighted mean in each cluster
    %************************************************
    
    % -- hotspots are weighted mean
    Wmap = map_hotspot./max(map_hotspot(:));
    hc_estimated = zeros(num_of_clusters,2);
    hcs_weights = zeros(num_of_clusters,1);
    hcs_ppmm = zeros(num_of_clusters,1);
    for i = 1:num_of_clusters

        Hcluster = [Hcells(hot_clusters==i,1),Hcells(hot_clusters==i,2)];

        Hcluster_wmean_row = 0;
        Hcluster_wmean_col = 0;
        Hcluster_wmean_wgt = 0;
        Hcluster_ppmm      = 0;
        for j = 1:size(Hcluster,1)
            Hcluster_wmean_row = Hcluster_wmean_row + ...
                Wmap(Hcluster(j,1),Hcluster(j,2))*Hcluster(j,1);
            Hcluster_wmean_col = Hcluster_wmean_col + ...
                Wmap(Hcluster(j,1),Hcluster(j,2))*Hcluster(j,2);
            Hcluster_wmean_wgt = Hcluster_wmean_wgt + ...
                Wmap(Hcluster(j,1),Hcluster(j,2));
            Hcluster_ppmm = Hcluster_ppmm +...
                map_hotspot(Hcluster(j,1),Hcluster(j,2));
        end


        Hcluster_wmean = [Hcluster_wmean_row/Hcluster_wmean_wgt,...
                          Hcluster_wmean_col/Hcluster_wmean_wgt];        
        hc_estimated(i,:) = round(Hcluster_wmean);
        hcs_weights(i) = Hcluster_wmean_wgt;       
        hcs_ppmm(i) = Hcluster_ppmm;
        
    end
    
    %hcw_before = [hc_estimated,hcs_weights];
    %************************************************
    % Step-4: ESTIMATE HOTSPOT CENTERS
    % move estimated sources to unoccupied space
    %************************************************
    [ hc_estimated ] = fMovingOccupiedHcsToFreeSpace( hc_estimated,...
                                                      OBSenv,...
                                                      map_env,...
                                                      cellsize_env,...
                                                      para_ );
    
    %hcw_after = [hc_estimated,hcs_weights];
    
    %find(hcw_after(:)~=hcw_before(:))
    %pause
    
    if visualize==1
        
        figure('name','clusters-2 with hotspots'); 
        imshow(map_env','InitialMagnification',350); hold on;
        set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
        pause(1)
        
        colors = lines(num_of_clusters);
        for i = 1:num_of_clusters
            plot(Hcells(hot_clusters==i,1)+0.5,Hcells(hot_clusters==i,2)+0.5,...
                 's',...
                 'color',colors(i,:),...
                 'MarkerFaceColor',colors(i,:),...
                 'MarkerSize',4); %5
        end
        
        % -- hotspots
        for i = 1:size(hc_estimated,1)
            plot(hc_estimated(i,1),hc_estimated(i,2),...
                 'o','color','k','MarkerFaceColor','r','MarkerSize',5); %7
        end
    end
    
    
    
    %************************************************
    % Step-5: ESTIMATE CENTERS
    % nearby centers are combined
    %************************************************   
    
    %-- merg all hcs close enough
    %dist_ind_kaamkay = 1;
    
    % -- temp fix
    %hc_estimated    
    %size(hc_estimated,1)
    if size(hc_estimated,1) < 2
        dist_ind_kaamkay = [];
    else 
        dist_ind_kaamkay = 1;
    end
    
    testind = 0;
    while numel(dist_ind_kaamkay)>=1

        testind = testind+1;
        %hcs_dist = fPathDist2(hotspots,hotspots,map_env);
        
        %-- combinations and pairwise distance
        %-----------------------------------------------
        scs_comb_ind = nchoosek(1:size(hc_estimated,1),2);
        scs_comb_dist = inf(size(scs_comb_ind,1),1);
        for i = 1:size(scs_comb_dist,1)
            %-- from -- to
            sc_from = hc_estimated(scs_comb_ind(i,1),:);
            sc_to   = hc_estimated(scs_comb_ind(i,2),:);
            %-- traveling dist
            %dist_this = fPathDist2(sc_from,sc_to,map_env);
            dist_this = pdist2(sc_from,sc_to,'euclidean');
            scs_comb_dist(i) = dist_this;
            %if dist_this<=para_.SCsCombiningDist_cells
            %    break;
            %end
        end
        
        [dist_val,dist_ind] = sort(scs_comb_dist);
        dist_ind_kaamkay = dist_ind(dist_val<=para_.SCsCombiningDist_cells);
        
        if numel(dist_ind_kaamkay)>=1
            
            sc_first = hc_estimated(scs_comb_ind(dist_ind_kaamkay(1),1),:);
            w_first  = hcs_weights(scs_comb_ind(dist_ind_kaamkay(1),1));
            ppm_first = hcs_ppmm(scs_comb_ind(dist_ind_kaamkay(1),1));
            
            sc_second = hc_estimated(scs_comb_ind(dist_ind_kaamkay(1),2),:);
            w_second  = hcs_weights(scs_comb_ind(dist_ind_kaamkay(1),2));
            ppm_second = hcs_ppmm(scs_comb_ind(dist_ind_kaamkay(1),2));
            
            w_mirged = w_first + w_second;
            
            sc_mirged_r = ((w_first*sc_first(1,1)) + (w_second*sc_second(1,1))) / w_mirged;
            sc_mirged_c = ((w_first*sc_first(1,2)) + (w_second*sc_second(1,2))) / w_mirged;
            
            sc_mirged = [sc_mirged_r,sc_mirged_c];
            
            ppm_mirged = ppm_first + ppm_second;
            
            hc_estimated(scs_comb_ind(dist_ind_kaamkay(1),1),:) = sc_mirged;
            hc_estimated(scs_comb_ind(dist_ind_kaamkay(1),2),:) = [];
            
            hcs_weights(scs_comb_ind(dist_ind_kaamkay(1),1)) = w_mirged;
            hcs_weights(scs_comb_ind(dist_ind_kaamkay(1),2)) = [];
            
            hcs_ppmm(scs_comb_ind(dist_ind_kaamkay(1),1)) = ppm_mirged;
            hcs_ppmm(scs_comb_ind(dist_ind_kaamkay(1),2)) = [];
            
        end
        
        %-- tmp fix
        %disp('temporary fix - line 1023');
        warning('This is temporary fix, fix it...')
        %hc_estimated
        if size(hc_estimated,1) < 2
            dist_ind_kaamkay = []
        end
        
    end
    
    
    %*******************************************
    %           nth LAYER
    %*******************************************
    
    
    hcs_prevCluster = hc_estimated;
    ws_prevCluster  = hcs_weights;
    ppmm_prevCluster  = hcs_ppmm;
    
    
    %{
    cond_nClustering = size(hcs_prevCluster,1)>1;    
    
    while cond_nClustering
    
        % - cluster linkage tree
        Z = linkage(hcs_prevCluster,'single');
        %figure; dendrogram(Z);

        % - high concentration clusters
        clusters_this = cluster(Z,'cutoff',8.5,'criterion','distance');

        % - number of clusters
        num_of_clusters_this = numel(unique(clusters_this));      
        
        % -- hotspots are weighted mean
        hcs_thisCluster = zeros(num_of_clusters_this,2);
        ws_thisCluster  = zeros(num_of_clusters_this,1);
        ppmm_thisCluster  = zeros(num_of_clusters_this,1);
        
        for i = 1:num_of_clusters_this

            candidateCells = [hcs_prevCluster(clusters_this==i,1),...
                              hcs_prevCluster(clusters_this==i,2)];
            candiateWeights = ws_prevCluster(clusters_this==i);
            candiatePPMM = ppmm_prevCluster(clusters_this==i);
            
            thisHc_r = 0;
            thisHc_c = 0;
            thisHc_w = 0;
            thisHc_ppmm = 0;
            for j = 1:size(candidateCells,1)
                %
                thisHc_r = thisHc_r + candiateWeights(j)*candidateCells(j,1);
                thisHc_c = thisHc_c + candiateWeights(j)*candidateCells(j,2);
                thisHc_w = thisHc_w + candiateWeights(j);
                thisHc_ppmm = thisHc_ppmm + candiatePPMM(j);
            end
            
            thisHc = [thisHc_r/thisHc_w,thisHc_c/thisHc_w];

            hcs_thisCluster(i,:) = round(thisHc);
            ws_thisCluster(i) = thisHc_w;
            ppmm_thisCluster(i) = thisHc_ppmm;
            
        end
        
        if visualize==1
        
            figure('name','clusters-n with hotspots'); 
            imshow(map_env','InitialMagnification',350); hold on;
            set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
            pause(1)

            colors = lines(num_of_clusters_this);
            for i = 1:num_of_clusters_this
                plot(hcs_prevCluster(clusters_this==i,1),...
                     hcs_prevCluster(clusters_this==i,2),...
                     's',...
                     'color','k',...
                     'MarkerFaceColor',colors(i,:),...
                     'MarkerSize',4); %5
            end

            % -- hotspots
            for i = 1:size(hcs_thisCluster,1)
                plot(hcs_thisCluster(i,1),hcs_thisCluster(i,2),...
                     'o','color','k','MarkerFaceColor','r','MarkerSize',5); %7
            end
        end
        
        cond_nClustering = size(hcs_thisCluster,1)>1 & ...
                           (size(hcs_prevCluster,1)-size(hcs_thisCluster,1))>0;
        
        hcs_prevCluster = hcs_thisCluster;
        ws_prevCluster = ws_thisCluster;
        ppmm_prevCluster = ppmm_thisCluster;
        
        %pause
        
    end
    
    if visualize==1
    pause(2)
    filename = 'hotspots-clusters-n-hcs';
    % print('-painters','-dpdf','-r500',filename); % pdf    
    export_fig([dir_.Solutions,filename], '-pdf','-r500')
    pause(1)
    end
    %}
    
    %hc_estimated = round(hcs_prevCluster);
    hc_estimated = round(hcs_prevCluster*100)/100;
    
    
    %*************************************************
    %-- high pass filtering of hcs based on weights
    %*************************************************
    total_weight = sum(ws_prevCluster);
    hcs_weights = ws_prevCluster;
    hcs_ppmm = ppmm_prevCluster;
    num_of_hcs = size(hc_estimated,1);
    
    %qualification_threshould = (total_weight/num_of_hcs*0.5)
    qualification_w = (total_weight/num_of_hcs*para_.HcsCutoffWeightFactor);    
        
    %hotspots = hotspots(hcs_weights>qualification_threshould,:)
    hc_estimated = hc_estimated( hcs_weights>qualification_w |...
                                 hcs_ppmm>para_.HcsCutoffPPMM,: );
    
    
    
    %-- moving hcs from occupied cells to the nearest unoccupied cells
    [ hc_estimated ] = fMovingOccupiedHcsToFreeSpace( hc_estimated,...
                                                      OBSenv,...
                                                      map_env,...
                                                      cellsize_env,...
                                                      para_ );
    hc_estimated = unique(hc_estimated,'rows');

    
    %{
    %-- merg all hcs close enough
    toMergHcs = [0,0];    
    while numel(toMergHcs)>=2 
        hcs_dist = fPathDist2(hc_estimated,hc_estimated,map_env);

        for i = 1:size(hc_estimated,1)
            hcs_dist(i,i) = inf;
        end
        [r,c] = ind2sub(size(hcs_dist), find(hcs_dist<10));

        toMergHcs = unique(sort([r,c],2),'rows');

        for i = 1:size(toMergHcs,1)
            hc_estimated(toMergHcs(i,1),1) = (hc_estimated(toMergHcs(i,1),1)+hc_estimated(toMergHcs(i,2),1))/2;
            hc_estimated(toMergHcs(i,1),2) = (hc_estimated(toMergHcs(i,1),2)+hc_estimated(toMergHcs(i,2),2))/2;

            hc_estimated(toMergHcs(i,2),1) = hc_estimated(toMergHcs(i,1),1);
            hc_estimated(toMergHcs(i,2),2) = hc_estimated(toMergHcs(i,1),2);

        end
        hc_estimated = round (hc_estimated);
        hc_estimated = unique(hc_estimated,'rows');
    end
    %}
    
else
    
    hc_estimated = Hcells;
    
end



%///////////////////////////////////////////// 
% HOTSPOTS OVER CLUSTERS
%/////////////////////////////////////////////
if visualize==1
figure('name','clusters with hotspots'); 
imshow(map_env','InitialMagnification',350); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
pause(1)
% color_clusters = lines(num_of_clusters);
for i = 1:num_of_clusters
    plot(Hcells(hot_clusters==i,1)+0.5,Hcells(hot_clusters==i,2)+0.5,...
         's',...
         'color',colors(num_of_clusters1+i,:),...
         'MarkerFaceColor',colors(num_of_clusters1+i,:),...
         'MarkerSize',2.75); %5
end
% -- hotspots
for i = 1:size(hc_estimated,1)
    plot(hc_estimated(i,1),hc_estimated(i,2),...
         'o','color','k','MarkerFaceColor','r','MarkerSize',3); %7
end

pause(2)
filename = sprintf('%s%s-hotspots-clusters',...
    dir_.hotspots,para_.ExperimentTitle);
%print('-painters','-dpdf','-r500',filename); % pdf
export_fig(filename, '-pdf','-r500');
pause(1)
end



%*******************************************************************
% HOTSPOTS OVER COARSE MAP
%*******************************************************************
if visualize==1
    
figure('name','coarse map with hotspots'); 
imshow(map_env','InitialMagnification',350); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])

% ------------- select color map for gdm --------------
gdm_colormap = flipud(hot(512)); % for main fig. 
max_concentration = max(map_hotspot(:)); % main fig.
delta = max_concentration/(size(gdm_colormap,1)-1);

[gd_row,gd_col] = size(map_hotspot);

for i = 1:gd_row
    for j = 1:gd_col
        % if positive concentration

        if map_hotspot(i,j)>=0
            
            if map_env(i,j) > 0

            % ---------------------- useful trash -------------------
            %col = gdm_colormap(round(map_recon(i,j)/delta)+1,:);

            % --- if concentration is greater than 1000 ppm, its still 1000 ppm
            %col = gdm_colormap( round(min(map_recon(i,j),3000)/delta)+1, :);
            %col = gdm_colormap( round(min(map_recon(i,j),1.78e+04)/delta)+1, :);
            %col = gdm_colormap( round(min(map_recon(i,j),inf)/delta)+1, :);
            % -------------------------------------------------------

            % ----------------------------------------------------
            % linear color code
            % ----------------------------------------------------
            %if strcmp(color_scale,'linear')
            % {
            linear_color_number = round( min(map_hotspot(i,j),max_concentration)/delta)+1;
            col = gdm_colormap( linear_color_number,:);                
            %}
            %end
            % ----------------------------------------------------

            %----- to avoid black borders of earch cell -----
            % {
            if sum(col) == 3
                col = col-(1e-10);
            end
            %}
            % ----- plot ------------
            %plot(j,i,'s','Color',col,'MarkerFaceColor',col,'MarkerSize',4);
            %plot(i,j,'s','Color',col,'MarkerFaceColor',col,'MarkerSize',4);
            %plot(cell_coords_x(i),cell_coords_y(j),'s','Color',col,'MarkerFaceColor',col,'MarkerSize',4);

            %plot(cell_coords_y(i)+0.5,cell_coords_x(j)+0.5,...
            %    's','Color',col,'MarkerFaceColor',col,'MarkerSize',5);

            %plot(cell_coords_x(j)+0.5,cell_coords_y(i)+0.5,...
            %    's','Color',col,'MarkerFaceColor',col,'MarkerSize',4.75);
            plot(i,j,...
                's',...
                'MarkerEdgeColor',col,...
                'MarkerFaceColor',col,...
                'MarkerSize',2.75); %4.85
            
            end
        end
    end
end   

pause(1)

for i = 1:size(hc_estimated,1)    
    plot(hc_estimated(i,1),hc_estimated(i,2),...
         'o','color','k','MarkerFaceColor','r','MarkerSize',3); %7
end

pause(2)
filename = sprintf('%s%s-hotspots-coarse-map',...
    dir_.hotspots,para_.ExperimentTitle);
%print('-painters','-dpdf','-r500',filename); % pdf
export_fig(filename, '-pdf','-r500');
pause(1)
end
 

% -----------------------------------------------------------------------
% hotspots on occupied cells
% -----------------------------------------------------------------------
[ hc_estimated ] = fMovingOccupiedHcsToFreeSpace( hc_estimated,...
                                                  OBSenv,...
                                                  map_env,...
                                                  cellsize_env,...
                                                  para_ );


if visualize==1
pause(1)
for i = 1:size(hc_estimated,1)
    plot(hc_estimated(i,1),hc_estimated(i,2),...
         'o','color','b','MarkerSize',7);
end
pause(1)
end




end


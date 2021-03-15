function [ M,...
           mean_map,...
           cell_coords_x,...
           cell_coords_y ] = fReconstruction( meas_files,...
                                              map_env,...
                                              para_,...
                                              dir_,...
                                              visualize )


% fReconstruction performs reconstruction from the integral concentration measurements.
%
% Author:  Asif Arain
% Project: Simulator for the Evaluation of Sensing Geometries (JFR-2016)


% CELL_SIZE = para_.ReconTomographyCellSizeM;

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



%*****************************************************************************************
% 
%                                   MEASUREMENTS
% 
%*****************************************************************************************

% -- measurement matrix
fprintf('- measurement matrix for line model.\n');

% meas_files = dir([dir_.logs,'measurement*_timestmp*_hotspot*_adaptive*']);
files_names = dir([dir_.logs,meas_files])

M = [];
for i = 1:size(files_names,1)
    files_names(i).name;
    
    M_current = dlmread([dir_.logs,files_names(i).name]);
    M = [M;M_current];
end



%*****************************************************************************************
% 
%                                   LINE MODEL
% 
%*****************************************************************************************
fprintf('- Line model for 2d gdm.\n');


recon_size = 'environment-map';

% [ A,...
%   ppmm,...
%   cell_coords_x,...
%   cell_coords_y,...
%   num_cells_x,...
%   num_cells_y ] = fLineModel( M,...
%                               map_env,...
%                               recon_size,...
%                               CELL_SIZE,...
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

% --- save cell_coords_x
% fileID = fopen([dir_.recon_cellx,x_coord_file],'wt'); fclose(fileID);
% dlmwrite([dir_.recon_cellx,x_coord_file],cell_coords_x,'-append','delimiter',' ');
% fileID = fopen([dir_.recon_cellx,x_coord_file],'a'); fclose(fileID);

% --- save cell_coords_y
% fileID = fopen([dir_.recon_celly,y_coord_file],'wt'); fclose(fileID);
% dlmwrite([dir_.recon_celly,y_coord_file],cell_coords_y,'-append','delimiter',' ');
% fileID = fopen([dir_.recon_celly,y_coord_file],'a'); fclose(fileID);


%*****************************************************************************************
% 
%                              GAS DISTRIBUTION MAP
% 
%*****************************************************************************************

fprintf('- 2d gas distribution map.\n');

[ mean_map ] = fGDM( A,...
                     ppmm,...
                     cell_coords_x,...
                     cell_coords_y,...
                     num_cells_x,...
                     num_cells_y,...
                     visualize );


                
%*****************************************************************************************
% 
%                              VISUALIZATION
% 
%*****************************************************************************************
                
                
if visualize ==1
    % ----------------    
    % -- splite M matrix
    start_x__all = M(:,3);
    start_y__all = M(:,4);
    end___x__all = M(:,5);
    end___y__all = M(:,6);
    plot([start_x__all';end___x__all'],[start_y__all';end___y__all'],'-g')
    %plot([start_y__all';end___y__all'],[start_x__all';end___x__all'],'-g')
    axis equal
    
    %figure; imshow(mean_map','InitialMagnification',800); hold on;
    %set(gca,'YDir','normal')
    %colormap('hot'); % not go back to gray scale.
    %%plot([start_x__all';end___x__all'],[start_y__all';end___y__all'],'-g')
    %plot([start_y__all';end___y__all'],[start_x__all';end___x__all'],'-g')    
    %pause(1)
    
end

if visualize ==1
    
  
    % - ground truth    
    gt_filename = sprintf('map_conc_%d.dat',para_.GroundTruthTimeStamp);
    map_conc_______gt = dlmread([dir_.GroundTruthLogs,gt_filename]);
    
    
    figure; hold on;
    image(map_conc_______gt,'CDataMapping','scaled')
    colormap hot; 
    %colormap gray; 
    colormap(flipud(colormap))
    colorbar;
    
    % ----------------    
    % -- splite M matrix
    start_x__all = M(:,3);
    start_y__all = M(:,4);
    end___x__all = M(:,5);
    end___y__all = M(:,6);
    plot([start_x__all';end___x__all'],[start_y__all';end___y__all'],'-g')
    %plot([start_y__all';end___y__all'],[start_x__all';end___x__all'],'-g')
    axis equal
    
    conf_x = unique(start_x__all)
    conf_y = unique(start_y__all)
end



end


function [ map_recon ] = fMapRecon( map_env,...
                                    OBS,...
                                    origin_env,...
                                    cell_coords_x,...
                                    cell_coords_y,...
                                    para_,...
                                    dir_,...
                                    visualize )
                                
%


% initialize gdm
if para_.ReconDetectionArtificial
    map_recon = zeros(size(map_env));
    
    % --- low concentration
    % generated a random concentrations
    lower_limit = 0;
    upper_limit = 50;
    %map_recon = (upper_limit-lower_limit).*rand(size(map_env)) + lower_limit;
    low_concentration = (upper_limit-lower_limit).*rand(numel(map_env)-numel(OBS.ind),1) + lower_limit;
    free_cells = find(map_env);
    map_recon(free_cells) = low_concentration; % low concentration
    
    % ---- high concentration
    lower_limit = 500;
    upper_limit = 5000;
    %map_recon = (upper_limit-lower_limit).*rand(size(map_env)) + lower_limit;
    %high_concentration = (upper_limit-lower_limit).*rand(15,1) + lower_limit;
    high_concentration = (upper_limit-lower_limit).*rand(para_.num_of_artificial_hotspots,1) + lower_limit;
          
    num_of_hotspots = para_.num_of_artificial_hotspots;%15;
    rnd_ind = randi([0 numel(free_cells)],1,num_of_hotspots);
        
    map_recon(free_cells(rnd_ind)) = high_concentration; % low concentration
    
    if visualize==1
    figure; imshow(map_recon)
    end
    
    
else
    
    load([dir_.gdplan_recon,'coarse_gdm.mat'],...
        'mean_map','cell_coords_x','cell_coords_y');
       
    %{
    % ------------- select color map for gdm --------------
    gdm_colormap = flipud(colormap('hot')); % for main fig. 
    % -----------------------------------------------------
    
    
    % ------ color step size --------------------------------
    max_concentration = max(mean_map(:)) % main fig.
    delta = max_concentration/(size(gdm_colormap,1)-1);
    
    
    mean_gd_plot = mean_map';
    [con_row,con_col] = size(mean_gd_plot)
    size(cell_coords_x)
    size(cell_coords_y)
    
    
    for i = 1:con_row%numel(cell_coords_x)%map_row
        for j = 1:con_col%numel(cell_coords_y)%map_col
            % if positive concentration
            
            if mean_gd_plot(i,j)>=0
                                
                linear_color_number = round( min(mean_gd_plot(i,j),max_concentration)/delta)+1;
                col = gdm_colormap( linear_color_number, :);
                
                %plot(cell_coords_x(i)+0.5,cell_coords_y(j)+0.5,...
                %    's','Color',col,'MarkerFaceColor',col,'MarkerSize',4.5);
                plot(i,j,'s','Color',col,'MarkerFaceColor',col,'MarkerSize',4.5);
                
                
            end
        end
    end
    %}
    
end

map_recon = mean_map';

%****************************************************************************
% 
%                 Adjust GDM according to environment map:
% 
%****************************************************************************

warning('NEED TO BE ADJUSTED');

if     strcmp(para_.TheCoarseMapIsFrom,'SimulatedMeasurement')
    
    cell_coords_x = cell_coords_x - origin_env(1);
    cell_coords_y = cell_coords_y - origin_env(2);
    
elseif strcmp(para_.TheCoarseMapIsFrom,'RobotMeasurement')
    %
    %    
else
    error('The sampling source for coarse map is unknown.');
end


% indices of robot origin in environment map.
% env_map_origin_ind_row = round((robot_origin_x-mapCut(2))/pix_cell);
% env_map_origin_ind_col = round((robot_origin_y-mapCut(1))/pix_cell);
env_map_origin_ind_row = origin_env(1);
env_map_origin_ind_col = origin_env(2);
env_map_origin2end_row = size(map_env,1)-env_map_origin_ind_row;
env_map_origin2end_col = size(map_env,2)-env_map_origin_ind_col;

% indices of mean map origin
% cell_coords_x
% cell_coords_y

% pause
% [~,mean_map_origin_ind_row] = min(abs(fliplr(cell_coords_y)))
[~,mean_map_origin_ind_row] = min(abs((cell_coords_x)));
[~,mean_map_origin_ind_col] = min(abs(cell_coords_y));

mean_map_origin2end_row = numel(cell_coords_x)-mean_map_origin_ind_row;
mean_map_origin2end_col = numel(cell_coords_y)-mean_map_origin_ind_col;

% indices differences between environment map and mean map.
deltaC_left   = env_map_origin_ind_col-mean_map_origin_ind_col;
deltaC_right  = env_map_origin2end_col-mean_map_origin2end_col;
deltaR_bottom = env_map_origin2end_row-mean_map_origin2end_row;
deltaR_top    = env_map_origin_ind_row-mean_map_origin_ind_row;


% add/delete columns from left of the map
if deltaC_left > 0
    map_recon = [zeros(size(map_recon,1),deltaC_left),map_recon];
elseif deltaC_left < 0
    map_recon(:,1:abs(deltaC_left)) = [];
    warning('First %d columns of gas concentration map are deleted',deltaC_left)
end

% add/delete coulmns from right of the map
if deltaC_right > 0
    map_recon = [map_recon,zeros(size(map_recon,1),deltaC_right)];
elseif deltaC_right < 0
    map_recon(:,end-abs(deltaC_right)+1:end) = [];
    warning('Last %d columns of gas concentration map are deleted',deltaC_right)
end

% add/delete rows from top of the map
if deltaR_top > 0
    map_recon = [zeros(deltaR_top,size(map_recon,2));map_recon];
elseif deltaR_top < 0
    map_recon(1:abs(deltaR_top),:) = [];
    warning('First %d rows of gas concentration map are deleted',deltaR_top)
end

% add/delete rows form the bottom of the map
if deltaR_bottom > 0
    map_recon = [map_recon;zeros(deltaR_bottom,size(map_recon,2))];
elseif deltaR_bottom < 0
    map_recon(end-abs(deltaR_bottom)+1:end,:) = [];
    warning('Last %d rows of gas concentration map are deleted',deltaR_bottom)
end

% size(map_recon)

if visualize==1
figure('name','gas distribution map'); 
imshow(map_recon','InitialMagnification',750); hold on;
set(gca,'YDir','normal'); colormap('gray'); %caxis([-2 2])
end



end


function [ A,...
           ppmm,...
           cell_coords_x,...
           cell_coords_y,...
           num_cells_x,...
           num_cells_y ] = fLineModelJFR( M,...
                                          recon_size,...
                                          CELL_SIZE,...
                                          para_,...
                                          dir_ )

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



% fprintf('\nLoading: %s', measuremets_file);
 

start_x = M(:,1);
start_y = M(:,2);
end___x = M(:,3);
end___y = M(:,4);
ppmm    = M(:,5);



% conf_num,...              % configuration number
% sensor___timestamp,...    % time stamp
% beam_initial_x,...        % beam start position along x-axis
% beam_initial_y,...        % beam start position along y-axis
% beam_end_____x,...        % beam final position along x-axis
% beam_end_____y,...        % beam final position along y-axis
% beam____length,...        % beam length
% ppm_m

map_filename = sprintf('%s_map_coverage.dat',para_.environment);
% map___environment = dlmread([dir_.env,para_.EnvironmentFigFile]);
map_env = dlmread([dir_.env,map_filename]);

map_env = 1-map_env;
[env_row,env_col] = size(map_env);
% pause




if strcmp(recon_size,'measurement-beams')

      
    
    switch para_.SensingSystem

            %=======================================
        case 'robot'
            %=======================================
            
            %{
            file__cellsize_env = sprintf('%s_cellsize_coverage.dat',para_.environment);   
            cellsize_env = load([dir_.env,file__cellsize_env]);
            %max_x = env_row * (cellsize_env/CELL_SIZE)
            
            min_x = max([round(min([min(start_x),min(end___x)])-CELL_SIZE-15),1]) * cellsize_env;
            max_x = min([round(max([max(start_x),max(end___x)])+CELL_SIZE+15),env_row]) * cellsize_env;
            min_y = max([round(min([min(start_y),min(end___y)])-CELL_SIZE-15),1]) * cellsize_env;
            max_y = min([round(max([max(start_y),max(end___y)])+CELL_SIZE+15),env_col]) * cellsize_env;
            %}
            
            min_x = max([round(min([min(start_x),min(end___x)])-15),1]);
            max_x = min([round(max([max(start_x),max(end___x)])+15),env_row]);
            min_y = max([round(min([min(start_y),min(end___y)])-15),1]);
            max_y = min([round(max([max(start_y),max(end___y)])+15),env_col]);
            
    
            %=======================================
        case 'simulated-sampling'
            %=======================================
            
            %{
            min_x = max([round(min([min(start_x),min(end___x)])-CELL_SIZE-15),1]);
            max_x = min([round(max([max(start_x),max(end___x)])+CELL_SIZE+15),env_row]);
            min_y = max([round(min([min(start_y),min(end___y)])-CELL_SIZE-15),1]);
            max_y = min([round(max([max(start_y),max(end___y)])+CELL_SIZE+15),env_col]);
            %}
            
            min_x = max([round(min([min(start_x),min(end___x)])-15),1]);
            max_x = min([round(max([max(start_x),max(end___x)])+15),env_row]);
            min_y = max([round(min([min(start_y),min(end___y)])-15),1]);
            max_y = min([round(max([max(start_y),max(end___y)])+15),env_col]);
            
            

    end
    

elseif strcmp(recon_size,'environment-map')
    
    %{
    min_x = 0;
    max_x = env_row; %100
    min_y = 0;
    max_y = env_col; %100
    %}
        
    %min_x = 0;
    %min_y = 0;
    min_x = 1;
    min_y = 1;
    
    %----------------------------------------
    
    switch para_.SensingSystem

            %=======================================
        case 'robot'
            %=======================================
            
            %{
            file__cellsize_env = sprintf('%s_cellsize_coverage.dat',para_.environment);   
            cellsize_env = load([dir_.env,file__cellsize_env]);
            %max_x = env_row * (cellsize_env/CELL_SIZE)
            max_x = env_row * cellsize_env;
            max_y = env_col * cellsize_env;
            
            num_cells_x = ceil((max_x - min_x)./CELL_SIZE);
            num_cells_y = ceil((max_y - min_y)./CELL_SIZE);

            cell_coords_x = (min_x/cellsize_env) + CELL_SIZE*(1:num_cells_x);
            cell_coords_y = (min_y/cellsize_env) + CELL_SIZE*(1:num_cells_y);
            %}
            
            max_x = env_row;
            max_y = env_col;

            %=======================================
        case 'simulated-sampling'
            %=======================================
            max_x = env_row;
            max_y = env_col;
               
            %{
            num_cells_x = ceil((max_x - min_x)./CELL_SIZE);
            num_cells_y = ceil((max_y - min_y)./CELL_SIZE);

            cell_coords_x = min_x + CELL_SIZE*(1:num_cells_x);
            cell_coords_y = min_y + CELL_SIZE*(1:num_cells_y);
            %}
    end
        
end



% num_cells_x = ceil((max_x - min_x)./CELL_SIZE)
% num_cells_y = ceil((max_y - min_y)./CELL_SIZE)
num_cells_x = ceil((max_x - min_x)./(cellsize_recon/cellsize_env))+1;
num_cells_y = ceil((max_y - min_y)./(cellsize_recon/cellsize_env))+1;

% cell_coords_x = min_x + CELL_SIZE*(1:num_cells_x)
% cell_coords_y = min_y + CELL_SIZE*(1:num_cells_y)
cell_coords_x = linspace(min_x,max_x,num_cells_x);
cell_coords_y = linspace(min_y,max_y,num_cells_y);

% size(cell_coords_x)
% size(cell_coords_y)

tot_num_cells = num_cells_x * num_cells_y;
tot_num_measurements = size(ppmm,1);
A = sparse(tot_num_measurements,tot_num_cells);

% size(A)
%build the design matrix

fprintf('---- Num of measurements to process: %d \n',tot_num_measurements);
for m = 1:tot_num_measurements
    %m
    %A(m,:) = GetMeasurementPath(start_x(m),end_x(m),start_y(m),end_y(m),cell_coords_x,cell_coords_y);
    
    %start_x(m)
    %end___x(m)
    %start_y(m)
    %end___y(m)
    %cell_coords_x
    %cell_coords_y
    something = GetMeasurementPath_asif(start_x(m),end___x(m),start_y(m),end___y(m),cell_coords_x,cell_coords_y);
    %something = GetMeasurementPath(start_x(m),end_x(m),start_y(m),end_y(m),cell_coords_x,cell_coords_y);
    %size(something)
    A(m,:) = something;
    
    %{
    if (mod(m,100) == 0)
        fprintf('%d measurements processed! \n',m);
    end
    %}
end
fprintf('---- %d measurements are processed.\n',m);

% name = strcat(strcat(strtok(INPUT_FILE_NAME,'.'),'_line_'),strrep(sprintf('%1.2f',CELL_SIZE),'.',''),'.mat');
% save([INPUT_FILE_LOCATION,name]);
% fprintf('---> Processing completed!, data saved as %s\n', name);

% save([M_LINE_LOC,M_LINE_FILE]);
% fprintf('---> Processing completed!, data saved as %s\n', M_LINE_FILE);

end

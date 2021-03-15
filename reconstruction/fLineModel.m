function [ A,...
           ppmm,...
           cell_coords_x,...
           cell_coords_y,...
           num_cells_x,...
           num_cells_y ] = fLineModel( M,...
                                       map_env,...
                                       recon_size,...
                                       CELL_SIZE,...
                                       para_,...
                                       dir_ )




 
% fprintf('\nLoading: %s', measuremets_file);
 

start_x = M(:,3);
start_y = M(:,4);
end___x = M(:,5);
end___y = M(:,6);
ppmm    = M(:,8);
% ppmm    = M(:,8)/1e+5;


% conf_num,...              % configuration number
% sensor___timestamp,...    % time stamp
% beam_initial_x,...        % beam start position along x-axis
% beam_initial_y,...        % beam start position along y-axis
% beam_end_____x,...        % beam final position along x-axis
% beam_end_____y,...        % beam final position along y-axis
% beam____length,...        % beam length
% ppm_m

% calculate boundaries
% min_x = min(vertcat(min(start_x),min(end___x)))-CELL_SIZE
% max_x = max(vertcat(max(start_x),max(end___x)))+CELL_SIZE
% min_y = min(vertcat(min(start_y),min(end___y)))-CELL_SIZE
% max_y = max(vertcat(max(start_y),max(end___y)))+CELL_SIZE

% min_x = min(vertcat(min(start_x),min(end_x)));
% max_x = max(vertcat(max(start_x),max(end_x)));
% min_y = min(vertcat(min(start_y),min(end_y)));
% max_y = max(vertcat(max(start_y),max(end_y)));

% varing x and y lengths. (working)
% min_x = min(vertcat(min(start_x),min(end___x)))-CELL_SIZE/2
% max_x = max(vertcat(max(start_x),max(end___x)))+CELL_SIZE/2
% min_y = min(vertcat(min(start_y),min(end___y)))-CELL_SIZE/2
% max_y = max(vertcat(max(start_y),max(end___y)))+CELL_SIZE/2
% pause
% min_x = 1
% max_x = 100
% min_y = 1
% max_y = 37



% min_x = max([round(min([min(start_x),min(end___x)])-CELL_SIZE/2),1])
% max_x = min([round(max([max(start_x),max(end___x)])+CELL_SIZE/2),env_row])
% min_y = max([round(min([min(start_y),min(end___y)])-CELL_SIZE/2),1])
% max_y = min([round(max([max(start_y),max(end___y)])+CELL_SIZE/2),env_col])


% --- tmp commented 20161120
% min_x = max([round(min([min(start_x),min(end___x)])-CELL_SIZE-3),1]);
% max_x = min([round(max([max(start_x),max(end___x)])+CELL_SIZE+3),env_row]);
% min_y = max([round(min([min(start_y),min(end___y)])-CELL_SIZE-3),1]);
% max_y = min([round(max([max(start_y),max(end___y)])+CELL_SIZE+3),env_col]);

% ---

% map___environment = dlmread([dir_.env,para_.EnvironmentFigFile]);
% map___environment = 1-map___environment;
% [env_row,env_col] = size(map___environment);


[env_row,env_col] = size(map_env);

recon_size

if strcmp(recon_size,'measurement-beams')

    min_x = max([round(min([min(start_x),min(end___x)])-CELL_SIZE-15),1]);
    max_x = min([round(max([max(start_x),max(end___x)])+CELL_SIZE+15),env_row]);
    min_y = max([round(min([min(start_y),min(end___y)])-CELL_SIZE-15),1]);
    max_y = min([round(max([max(start_y),max(end___y)])+CELL_SIZE+15),env_col]);

elseif strcmp(recon_size,'environment-map')
    
    min_x = 0;
    max_x = env_row; %100
    min_y = 0;
    max_y = env_col; %100
    
end



% fixed as derscribed in the environment model
% load(['trinchetto-gas_tomography_simulator/Models/EnvModel_',EnvModel])
% load(['trinchetto-gas_tomography_simulator/Models/EnvModel/',EnvModel])

% load([ENV_LOC,ENV_FILE])
% min_x = round(min(end_x))-7 %min(environment_model.cell_coords_x);
% max_x = round(max(end_x))+7 %max(environment_model.cell_coords_x);
% min_y = round(min(end_y))-7 %min(environment_model.cell_coords_y);
% max_y = round(max(end_x))+7 %max(environment_model.cell_coords_y);
% roundingthreshold = 0.001;
% min_x = round(min(min(start_x),min(end_x))-roundingthreshold+0.5);
% max_x = round(max(max(start_x),max(end_x))-roundingthreshold+0.5);

% disp('--- to be corrected [min_x,max_x,min_y,max_y], look rounding off!')
% min_x = round(min( min(start_x) , min(end_x) ))
% max_x = round(max( max(start_x) , max(end_x) ))
% min_y = round(min( min(start_y) , min(end_y) ))
% max_y = round(max( max(start_y) , max(end_y) ))
% Note: setting a rounding threshold can be set as,
% round(x-roundingtreshold+0.5)
% However, this does not work for negative numbers.



% min_x = -8-CELL_SIZE/2
% max_x = 8+CELL_SIZE/2
% min_y = -8-CELL_SIZE/2
% max_y = 8+CELL_SIZE/2

% min_x = 35; %round(min(end___x))-7 %min(environment_model.cell_coords_x);
% max_x = 66; %round(max(end___x))+7 %max(environment_model.cell_coords_x);
% min_y = 32; %round(min(end___y))-7 %min(environment_model.cell_coords_y);
% max_y = 66; %round(max(end___x))+7 %max(environment_model.cell_coords_y);

% min_x = round(min(end___x))-7 %min(environment_model.cell_coords_x);
% max_x = round(max(end___x))+7 %max(environment_model.cell_coords_x);
% min_y = round(min(end___y))-7 %min(environment_model.cell_coords_y);
% max_y = round(max(end___x))+7 %max(environment_model.cell_coords_y);


num_cells_x = ceil((max_x - min_x)./CELL_SIZE);
num_cells_y = ceil((max_y - min_y)./CELL_SIZE);

% fprintf('\n ---> Min, Max X coordinates (%f,%f)',min_x,max_x);
% fprintf('\n ---> Min, Max Y coordinates (%f,%f)',min_y,max_y);
% fprintf('\n ---> Total cells %d \n',num_cells_x*num_cells_y);

 
% cell_coords_x = min_x + CELL_SIZE*(0:num_cells_x-1);
% cell_coords_y = min_y + CELL_SIZE*(0:num_cells_y-1);

cell_coords_x = min_x + CELL_SIZE*(1:num_cells_x);
cell_coords_y = min_y + CELL_SIZE*(1:num_cells_y);

% cell_coords_x = min_x + CELL_SIZE*(0:num_cells_x)
% cell_coords_y = min_y + CELL_SIZE*(0:num_cells_y)



tot_num_cells = num_cells_x * num_cells_y;
tot_num_measurements = size(ppmm,1);
A = sparse(tot_num_measurements,tot_num_cells);

% size(A)
%build the design matrix
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
    
    if (mod(m,100) == 0)
        fprintf('%d measurements processed! \n',m);
    end
end

% name = strcat(strcat(strtok(INPUT_FILE_NAME,'.'),'_line_'),strrep(sprintf('%1.2f',CELL_SIZE),'.',''),'.mat');
% save([INPUT_FILE_LOCATION,name]);
% fprintf('---> Processing completed!, data saved as %s\n', name);

% save([M_LINE_LOC,M_LINE_FILE]);
% fprintf('---> Processing completed!, data saved as %s\n', M_LINE_FILE);

end

function [ mean_map ] = fGDM( A,...
                              ppmm,...
                              cell_coords_x,...
                              cell_coords_y,...
                              num_cells_x,...
                              num_cells_y,...
                              visualize )
% f2dMapLineModel generates 2d GDM with line model.



LAMBDA = 1e-10;
% LAMBDA = 1e-3;
% LAMBDA = 1e-1;
% LAMBDA = 0;


% A_FILE = (['M_',num2str(SIMULATION_FILE),'_line_',strrep(sprintf('%1.2f',CELL_SIZE),'.',''),'.mat']);
% A_FILE = strrep(A_FILE,'SimCase_','');
% load([INPUT_FILE_LOCATION,A_FILE]);

% load([M_LINE_LOC,M_LINE_FILE])

% fprintf('\nLoading A matrix file: %s \n',A_FILE);
% fprintf('\nLoading A matrix file: %s \n',M_LINE_FILE);

% [mean_map_vector, variance_map_vector] = CalculateGDM(A,ppmm,LAMBDA);
%[mean_map_vector,~] = CalculateGDM(A,ppmm,LAMBDA);
[mean_map_vector] = CalculateGDM(A,ppmm,LAMBDA);
%[X,Y]               = meshgrid(cell_coords_x,cell_coords_y); % -- tempo
mean_map          = reshape(mean_map_vector,num_cells_y,num_cells_x);

%mean_map            = mean_map';

% variance_map = reshape(variance_map_vector,num_cells_y,num_cells_x);

% cell_coords_x
% cell_coords_y

% max(mean_map(:))
if visualize == 1
    figure; hold on;
    %pcolor(X,Y,mean_map);
    image(cell_coords_x,cell_coords_y,mean_map,'CDataMapping','scaled')
    % title(sprintf('Mean Map, file %s', A_FILE));
    colormap hot; 
    %colormap gray; 
    colormap(flipud(colormap))
    colorbar;
    
    %pause
end



end
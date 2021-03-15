function [ M,...
           mean_map,...
           cell_coords_x,...
           cell_coords_y ] = fReconstructionFixed( para_,...
                                                   dir_,...
                                                   visualize )


% fReconstruction performs reconstruction from the integral concentration measurements.
%
% Author:  Asif Arain
% Project: Simulator for the Evaluation of Sensing Geometries (JFR-2016)



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







switch para_.SensingSystem
        
        %=======================================
    case 'robot'
        %=======================================
        
        
        %--- robot logs to M matrix for all conf.
        %====================================================
        fprintf('---> Generating measurement matrix (M) for all the conf.\n');
        
        meas_files = 'measurements_conf*.dat*';
        files_names = dir([dir_.logs,meas_files]);
        
        
        % {
        for i = 1:size(files_names,1)
            measure_file = sprintf('measurements_conf%d.dat',i);
            % visualize = 1;
            [ M_m,M_cell ] = fRobotLogs2M( measure_file,...
                                           para_,...
                                           dir_,...
                                           visualize );
        
            Mfile = sprintf('M_conf%02d.mat',i);
            save([dir_.MeasurementLogs,Mfile],'M_m','M_cell');
        end
        %}


        %--- M matrix for all conf.
        %====================================================

        % -- measurement matrix
        fprintf('---> Combining all measurements (M).\n');

        M = [];
        for num = 1:size(files_names,1)

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
    case 'simulated-sampling'
        %=======================================
        
        
        % -- measurement matrix
        fprintf('- measurement matrix for line model.\n');

        % meas_files = dir([dir_.logs,'measurement*_timestmp*_hotspot*_adaptive*']);
        meas_files = 'measurement_conf*.dat*';
        files_names = dir([dir_.logs,meas_files])

        M = [];
        for i = 1:size(files_names,1)
            files_names(i).name;

            M_current = dlmread([dir_.logs,files_names(i).name]);
            M = [M;M_current];
        end
        
        %-- final M matrix (new)
        M = M(:,[3,4,5,6,8]);
end



%	LINE MODEL 
%====================================================
fprintf('---> Generating line model.\n');

recon_size = 'environment-map';
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



%--- save cell_coords_x
x_coord_file = sprintf('x_coord.dat');
fileID = fopen([dir_.Recon,x_coord_file],'wt'); fclose(fileID);
dlmwrite([dir_.Recon,x_coord_file],cell_coords_x,'-append','delimiter',' ');
fileID = fopen([dir_.Recon,x_coord_file],'a'); fclose(fileID);

%--- save cell_coords_y
y_coord_file = sprintf('y_coord.dat');
fileID = fopen([dir_.Recon,y_coord_file],'wt'); fclose(fileID);
dlmwrite([dir_.Recon,y_coord_file],cell_coords_y,'-append','delimiter',' ');
fileID = fopen([dir_.Recon,y_coord_file],'a'); fclose(fileID);


 
%----  GAS DISTRIBUTION MAP
%====================================================

fprintf('---> Building gas distribution map.\n');

[ mean_map ] = fGDM( A,...
                     ppmm,...
                     cell_coords_x,...
                     cell_coords_y,...
                     num_cells_x,...
                     num_cells_y,...
                     visualize );


%--- save cell_coords_y
% recon_file = sprintf('reconstruction.dat');
% fileID = fopen([dir_.Recon,recon_file],'wt'); fclose(fileID);
% dlmwrite([dir_.Recon,recon_file],mean_map,'-append','delimiter',' ');
% fileID = fopen([dir_.Recon,recon_file],'a'); fclose(fileID);                 
                



map_recon = mean_map';

%--- save mean map all
recon_file = sprintf('reconstruction.dat');
fileID = fopen([dir_.Recon,recon_file],'wt'); fclose(fileID);
dlmwrite([dir_.Recon,recon_file],map_recon,'-append','delimiter',' ');
fileID = fopen([dir_.Recon,recon_file],'a'); fclose(fileID);

map_recon_col = mean_map;
gdm_colormap = flipud(autumn(512)); % for main fig.
max(map_recon_col(:))
max_concentration = min( [max(map_recon_col(:)),para_.PublishUpperPPM] ) % main fig.
delta = max_concentration/(size(gdm_colormap,1)-1);
% linear_color_number = round(map_recon_col(:)./delta)+1;
linear_color_number = round(min(map_recon_col(:),max_concentration)./delta)+1;
 

colR = round(gdm_colormap(linear_color_number,1)*255);
colG = round(gdm_colormap(linear_color_number,2)*255);
colB = round(gdm_colormap(linear_color_number,3)*255);

thisfile = sprintf('reconstructionColorR.dat');
fileID = fopen([dir_.Recon,thisfile],'wt'); fclose(fileID);
dlmwrite([dir_.Recon,thisfile],colR,'-append','delimiter',' ');
fileID = fopen([dir_.Recon,thisfile],'a'); fclose(fileID);

thisfile = sprintf('reconstructionColorG.dat');
fileID = fopen([dir_.Recon,thisfile],'wt'); fclose(fileID);
dlmwrite([dir_.Recon,thisfile],colG,'-append','delimiter',' ');
fileID = fopen([dir_.Recon,thisfile],'a'); fclose(fileID);

thisfile = sprintf('reconstructionColorB.dat');
fileID = fopen([dir_.Recon,thisfile],'wt'); fclose(fileID);
dlmwrite([dir_.Recon,thisfile],colB,'-append','delimiter',' ');
fileID = fopen([dir_.Recon,thisfile],'a'); fclose(fileID);



% %-- copy to the ROS folder for publication
% filesFrom = dir_.Recon;
% filesTo = sprintf('%sone-step-exploration/',dir_.ROSLogs);
% 
% %-- copy command
% syncCommand = sprintf('rsync -a --ignore-existing %s %s',filesFrom,filesTo);
% 
% %-- copy log files
% system(syncCommand);



%-- copy to the ROS folder for publication
switch para_.SensingSystem
    case 'robot'
        
        filesFrom = dir_.Recon;
        filesTo = dir_.ROSLogsThis; %sprintf('%sone-step-exploration/',dir_.ROSLogs);
        
        %-- copy command
        syncCommand = sprintf('rsync -a %s %s',filesFrom,filesTo);

        %-- copy log files
        system(syncCommand);
        
        fprintf('---- The reconstruction is synchronized to display.\n');
        
    case 'simulated-sampling'
        %
end



%---- VISUALIZATION
%====================================================
                
if visualize == 1
    
    %map_env = map_env;
    
    figure('name','reconstrucion on environment map'); 
    imshow(map_env','InitialMagnification',800); hold on;
    set(gca,'YDir','normal');
    colormap('gray'); caxis([-1 1.25]);
           
    map_recon = mean_map';
    
    
    % ------------- select color map for gdm --------------
    gdm_colormap = flipud(hot(512)); % for main fig. 
           
    % ------ color step size --------------------------------    
    max_concentration = max(map_recon(:)); % main fig.
    delta = max_concentration/(size(gdm_colormap,1)-1);
    
    
    
    for i = 1:numel(cell_coords_x)%1:numel(cell_coords_y)%map_row
        for j = 1:numel(cell_coords_y)%1:numel(cell_coords_x)%map_col
            % if positive concentration
            
            if map_recon(i,j)>=0
                
                if map_env(i,j)>0
                % ---------------------- useful trash -------------------
                %col = gdm_colormap(round(map_gd(i,j)/delta)+1,:);
                
                % --- if concentration is greater than 1000 ppm, its still 1000 ppm
                %col = gdm_colormap( round(min(map_gd(i,j),3000)/delta)+1, :);
                %col = gdm_colormap( round(min(map_gd(i,j),1.78e+04)/delta)+1, :);
                %col = gdm_colormap( round(min(map_gd(i,j),inf)/delta)+1, :);
                % -------------------------------------------------------
                
                % ----------------------------------------------------
                % linear color code
                % ----------------------------------------------------
                
                linear_color_number = round( min(map_recon(i,j),max_concentration)/delta)+1;
                col = gdm_colormap( linear_color_number,:);                
                
                
                %----- to avoid black borders of earch cell -----
                % {
                if sum(col) == 3
                    col = col-(1e-10);
                end
                %}
                % ----- plot ------------
                plot(cell_coords_x(i)+0.0,cell_coords_y(j)+0.0,...
                    's',...
                    'MarkerEdgeColor',col,...
                    'MarkerFaceColor',col,...
                    'MarkerSize',5);
                end
            end
        end
    end

end


                
if visualize == 1
    %pause
    % ----------------    
    % -- splite M matrix
    start_x__all = M(:,1);
    start_y__all = M(:,2);
    end___x__all = M(:,3);
    end___y__all = M(:,4);
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

% if visualize == 1
%     
%   
%     % - ground truth    
%     gt_filename = sprintf('map_conc_%d.dat',para_.GroundTruthTimeStamp);
%     map_conc_______gt = dlmread([dir_.GroundTruthLogs,gt_filename]);
%     
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



end
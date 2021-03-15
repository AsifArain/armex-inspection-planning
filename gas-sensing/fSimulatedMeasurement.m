function fSimulatedMeasurement( sensing_conf,...
                                measu_file,...
                                map_env,...
                                para_,...
                                dir_,...
                                visualize )

% fMeasurements collects integral concentrations.
%
% Author:  Asif Arain
% Project: Simulator for the Evaluation of Sensing Geometries (JFR-2016)

% visualize = 1;
% ____________________________________________________________________ %
%                                                                      %
%                           environment map                            %
% ____________________________________________________________________ %

% map___environment = dlmread([dir_.env,para_.EnvironmentFigFile]);
% map___environment = 1-map___environment;




if visualize == 1
    % figure; imshow(map___environment,'InitialsensorMagnification',700); hold on;

    figure; hold on;
    [map_row,map_col] = size(map_env);
    for i = 1:map_row
        for j = 1:map_col
            % if obstacle/free
            if ~map_env(i,j) % occupied
                col = [0.4 0.4 0.4];
                faceCol = [0.4 0.4 0.4];
            elseif map_env(i,j) % free
                col = [0.9 0.9 0.9];
                faceCol = [1 1 1];
                %col = [1 1 1];
            end        
            plot(i+0.5,j+0.5,'S','Color',col,...
                             'MarkerFaceColor',faceCol,...
                             'MarkerSize',20);
        end
    end

    set(gca,'xtick',1:1:map_row);
    set(gca,'ytick',1:1:map_col);
    grid on; axis equal;
end


% ____________________________________________________________________ %
%                                                                      %
%     sensing configuration - collecting integral measurements         %
%                                                                      %
% ____________________________________________________________________ %

% --- write to a file
% initialize and clear file 
fileID = fopen([dir_.logs,measu_file],'wt'); fclose(fileID); % writing to follow


for conf_num = 1:size(sensing_conf,1)
    
    % -- sensing parameters
    sensor___________x = sensing_conf(conf_num,1);
    sensor___________y = sensing_conf(conf_num,2);
    sensor_orientation = sensing_conf(conf_num,3);
    sensor_________fov = sensing_conf(conf_num,4);
    sensor_______range = sensing_conf(conf_num,5);
    sensor___samp_freq = sensing_conf(conf_num,6);
    sensor___timestamp = sensing_conf(conf_num,7);
    
    % -- beam directions
    sensor__directions = linspace(sensor_orientation-sensor_________fov/2,...
                                  sensor_orientation+sensor_________fov/2,...
                                  sensor___samp_freq);
    % adjust directions >= 360 deg
    sensor__directions(sensor__directions>=360) =...
        sensor__directions(sensor__directions>=360)-360;
    
    % adjust directions < 0 deg
    sensor__directions(sensor__directions<0) =...
        sensor__directions(sensor__directions<0)+360;
    
    % gas concentration map at desired time stamp
    % ____________________________________________
    %{
    if sensor___timestamp ~= 0 % if time stamp is not 0 (0 is for avg. map)
        
        %concentrationMapName = sprintf('map_conc_stamp_%02d.dat',sensor___timestamp);
        %concentrationMapName = sprintf('map_conc_stmp_%03d.dat',sensor___timestamp);
        %concentrationMapName = sprintf('map_conc_%d.dat',sensor                      
        %[dir_.MeasurementLogsLoc,concentrationMapName]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             ___timestamp)
        %map_concentration = dlmread([dir_.MeasurementLogsLoc,concentrationMapName]);
        
        gt_filename = sprintf('map_conc_%d.dat',para_.GroundTruthTimeStamp);
        map_concentration = dlmread([dir_.GroundTruthLogs,gt_filename]);
        
    elseif sensor___timestamp == 0 % use avg concentration map
        
        %map_concentration = dlmread([dir_.maps,'map_conc_______gt.dat']);
        %map_concentration = dlmread([dir_.maps,'map_conc_000to000.dat']);
        map_concentration = dlmread([dir_.GroundTruthLogs,'map_conc_0.dat']);
        
    end
    %}
    map_concentration = dlmread([dir_.gt,'map_con.dat']);
    
    if visualize == 1 && conf_num == 1
        fVisualizeCon( map_concentration );
    end
    
    
    % -- visualize gas concentration map
    %{
    if visualize == 1 && conf_num == 1
        
        [map_row,map_col] = size(map_concentration);

        % pl select color scale
        color_scale = 'linear'; % insets
        %color_scale = 'nonlinear'; % main fig.
        
        % ------------- select color map for gdm --------------
        % colormap hot; 
        % colormap gray; 
        % colormap(flipud(colormap('hot')));
        gdm_colormap = colormap(flipud(colormap('hot'))); % for main fig. 
        %gdm_colormap = colormap(flipud(colormap('bone')));
        %gdm_colormap = colormap(flipud(colormap('copper')));
        %gdm_colormap = colormap(flipud(colormap('pink')));
        %gdm_colormap = colormap(flipud(colormap('cool')));
        %gdm_colormap = colormap(flipud(colormap('spring')));
        %gdm_colormap = colormap(flipud(colormap('summer'))); % for insets
        %gdm_colormap = colormap(flipud(colormap('autumn')));
        %gdm_colormap = colormap('summer');
        %gdm_colormap = colormap('autumn');
        colormap('gray'); % not go back to gray scale.
        %caxis([-5 2]) % for resulting gdm
        caxis([-12 4]) % for sensor placement solution
        %caxis([-10 1.5]) % for gdm insets
        % load gdm_colormap.mat
        % delta = max(map_gd(:))/(size(gdm_colormap,1)-1);
        %load publish/gdm_colormap.mat
        % -----------------------------------------------------
        
        % ------ color step size --------------------------------
        %max_concentration = 8000; % inset gdm
        max_concentration = max(map_concentration(:)); % main fig.
        delta = max_concentration/(size(gdm_colormap,1)-1);    
        %delta = max(map_gd(:))/(size(gdm_colormap,1)-1);
        %delta = 500/(size(gdm_colormap,1)-1); % lets limit to 500ppm

        % lets plot for the scale of 500 ppm.
        %delta = 3000/(size(gdm_colormap,1)-1);
        %max(map_concentration(:))
        %pause
        %delta = 1.78e+04/(size(gdm_colormap,1)-1)
        
        % ------------------- nonlinear color code ---------------------  
        if strcmp(color_scale,'nonlinear')
            %b = (exp(256/1))/255;
            %a = 256 /(log(b*256));

            x = 1:1:256;
            %epsilon = 1/(exp(1)-1); % standard
            %epsilon = 1/(exp(0.05)-1); % slower
            epsilon = 1/(exp(0.10)-1); % slower
            %epsilon = 1/(exp(0.15)-1); % slower
            %epsilon = 1/(exp(0.20)-1); % slower
            %epsilon = 1/(exp(0.90)-1); % slower
            %epsilon = 1/(exp(1.5)-1); % standard
            normlogsum_penalty = log(1+(abs(x)/epsilon));
            %normlogsum_penalty = 1.05.^(1+(abs(x)/epsilon));
            normlogsum_penalty = normlogsum_penalty/(normlogsum_penalty(x==max(x))); % normalized
            normlogsum_penalty = round(256*normlogsum_penalty);
            normlogsum_penalty(1)
            normlogsum_penalty(end)
            %figure; plot(x,normlogsum_penalty); % to test

            % publish nonlinear steps
            linear_percentage    = 10:10:100;
            linear_percentage_ind  = round((linear_percentage/100)*numel(x));
            nonlinear_percentage_ind = normlogsum_penalty(linear_percentage_ind);
            nonlinear_percentage = (nonlinear_percentage_ind/numel(x))*100;
            nonlinear_percentage = round(nonlinear_percentage);
            %disp('nonlinear_percentage = ',num2str(nonlinear_percentage))
            disp_stng = ['nonlinear color code scale is: ',num2str(nonlinear_percentage)];
            %disp(disp_stng)
            disp(disp_stng)
        end
        % ----------------------------------------------------------------
        
        
        
        for i = 1:map_row
            for j = 1:map_col
                % if positive concentration
                if map_concentration(i,j)>0
                    
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
                    if strcmp(color_scale,'linear')
                    linear_color_number = round( min(map_concentration(i,j),max_concentration)/delta)+1;
                    col = gdm_colormap( linear_color_number, :);
                    end
                    % ----------------------------------------------------

                    % ----------------------------------------------------
                    % nonlinear color code
                    % ----------------------------------------------------
                    if strcmp(color_scale,'nonlinear')
                    linear_color_number = round( min(map_concentration(i,j),max_concentration)/delta)+1;
                    %nonlinear_color_number = round( a*log(b*linear_color_number) );
                    %nonlinear_color_number = linear_color_number;
                    nonlinear_color_number = normlogsum_penalty(linear_color_number);
                    %if nonlinear_color_number <1
                    %    nonlinear_color_number = 1;
                    %end
                    col = gdm_colormap( nonlinear_color_number, :);                
                    end
                    % ----------------------------------------------------

                    % ----- plot ------------
                    %plot(j,i,'s','Color',col,'MarkerFaceColor',col,'MarkerSize',4);
                    %plot(i,j,'s','Color',col,'MarkerFaceColor',col,'MarkerSize',4);
                    %plot(i-0.5,j-0.5,'s','Color',col,'MarkerFaceColor',col,'MarkerSize',4);
                    plot(i+0.5,j+0.5,'s','Color',col,'MarkerFaceColor',col,'MarkerSize',15);
                end
            end
        end       
    end
    %}

    
    % -- integral measurements along beams
    for num = 1:sensor___samp_freq
        
        % -- beam parameters
        beam_initial_x = sensor___________x;
        beam_initial_y = sensor___________y;
        beam_direction = sensor__directions(num);
        beam_____range = sensor_______range;
        
        % -- beam on a simple grid
        [ beam_end_____x,...
          beam_end_____y,...
          beam_voxel___list,...
          beam_voxel_points,...
          beam_voxel_length] = fBeamOnSimpleGrid( beam_initial_x,...
                                                  beam_initial_y,...
                                                  beam_direction,...
                                                  beam_____range,...
                                                  map_env );
        
        % -- integral concentations
        [ ppm_m,...
          beam_voxel____end,...
          beam_voxel_end_xy,...
          beam_voxel___PPMM] = fBeamPPMM( beam_voxel___list,...
                                          beam_voxel_length,...
                                          beam_voxel_points,...
                                          map_env,...
                                          map_concentration);
                        
        % -- beam length in occluded map (beam range is the max length
        % without obstacles)
        %beam_voxel___list
        %beam_voxel_end_xy
        %beam_initial_x
        %beam_initial_y
        beam____length = sqrt( (beam_voxel_end_xy(1)-beam_initial_x)^2 +...
                               (beam_voxel_end_xy(2)-beam_initial_y)^2 );
        
        % -- write integral concentration data to the file
        dlmwrite([dir_.logs,measu_file],...
                 [conf_num,...              % configuration number
                  sensor___timestamp,...    % time stamp
                  beam_initial_x,...        % beam start position along x-axis
                  beam_initial_y,...        % beam start position along y-axis
                  beam_voxel_end_xy(1),...  % beam final position along x-axis
                  beam_voxel_end_xy(2),...  % beam final position along y-axis
                  beam____length,...        % beam length
                  ppm_m],...                % integral concentrtion value (ppm.m)
                  '-append','delimiter',' ');
        fileID = fopen([dir_.logs,measu_file],'a'); fclose(fileID);
                
        
        % -- visualize beams
        if visualize == 1
            
            % -- voxels from start to end in occuleded map
            h_voxelPPMM = plot(beam_voxel___PPMM(1,:)+0.5,...
                               beam_voxel___PPMM(2,:)+0.5,'sr');
            
            % -- beam full range
            %plot([beam_initial_x,beam_end_____x],[beam_initial_y,beam_end_____y],...
            %    '-r','MarkerFaceColor','b','MarkerSize',5);
            
            % -- beam start position to end position in occluded map            
            plot([beam_initial_x,beam_voxel_end_xy(1)],...
                 [beam_initial_y,beam_voxel_end_xy(2)],...
                 '-','Color',[176,196,222]/256,'MarkerSize',5);
            
            % - beam end positions in occluded map
            %plot(beam_voxel_end_xy(1),beam_voxel_end_xy(2),...
            %    'xb','MarkerFaceColor','b','MarkerSize',5);
            
            delete(h_voxelPPMM)
        end
        
        
        
    end
    %fprintf('- Integral measurements collected for conf #%d.\n',conf_num);
    
    if visualize == 1
        % -- mark configuration
        plot(sensor___________x,sensor___________y,...
             'ob','MarkerFaceColor','b','MarkerSize',5);
    end
    
    
        
end

   
end

function fVisualizeCon( map_concentration )

        [map_row,map_col] = size(map_concentration);

        % pl select color scale
        color_scale = 'linear'; % insets
        %color_scale = 'nonlinear'; % main fig.
        
        % ------------- select color map for gdm --------------
        % colormap hot; 
        % colormap gray; 
        % colormap(flipud(colormap('hot')));
        gdm_colormap = colormap(flipud(colormap('hot'))); % for main fig. 
        %gdm_colormap = colormap(flipud(colormap('bone')));
        %gdm_colormap = colormap(flipud(colormap('copper')));
        %gdm_colormap = colormap(flipud(colormap('pink')));
        %gdm_colormap = colormap(flipud(colormap('cool')));
        %gdm_colormap = colormap(flipud(colormap('spring')));
        %gdm_colormap = colormap(flipud(colormap('summer'))); % for insets
        %gdm_colormap = colormap(flipud(colormap('autumn')));
        %gdm_colormap = colormap('summer');
        %gdm_colormap = colormap('autumn');
        colormap('gray'); % not go back to gray scale.
        %caxis([-5 2]) % for resulting gdm
        caxis([-12 4]) % for sensor placement solution
        %caxis([-10 1.5]) % for gdm insets
        % load gdm_colormap.mat
        % delta = max(map_gd(:))/(size(gdm_colormap,1)-1);
        %load publish/gdm_colormap.mat
        % -----------------------------------------------------
        
        % ------ color step size --------------------------------
        %max_concentration = 8000; % inset gdm
        max_concentration = max(map_concentration(:)); % main fig.
        delta = max_concentration/(size(gdm_colormap,1)-1);    
        %delta = max(map_gd(:))/(size(gdm_colormap,1)-1);
        %delta = 500/(size(gdm_colormap,1)-1); % lets limit to 500ppm

        % lets plot for the scale of 500 ppm.
        %delta = 3000/(size(gdm_colormap,1)-1);
        %max(map_concentration(:))
        %pause
        %delta = 1.78e+04/(size(gdm_colormap,1)-1)
        
        % ------------------- nonlinear color code ---------------------  
        if strcmp(color_scale,'nonlinear')
            %b = (exp(256/1))/255;
            %a = 256 /(log(b*256));

            x = 1:1:256;
            %epsilon = 1/(exp(1)-1); % standard
            %epsilon = 1/(exp(0.05)-1); % slower
            epsilon = 1/(exp(0.10)-1); % slower
            %epsilon = 1/(exp(0.15)-1); % slower
            %epsilon = 1/(exp(0.20)-1); % slower
            %epsilon = 1/(exp(0.90)-1); % slower
            %epsilon = 1/(exp(1.5)-1); % standard
            normlogsum_penalty = log(1+(abs(x)/epsilon));
            %normlogsum_penalty = 1.05.^(1+(abs(x)/epsilon));
            normlogsum_penalty = normlogsum_penalty/(normlogsum_penalty(x==max(x))); % normalized
            normlogsum_penalty = round(256*normlogsum_penalty);
            normlogsum_penalty(1)
            normlogsum_penalty(end)
            %figure; plot(x,normlogsum_penalty); % to test

            % publish nonlinear steps
            linear_percentage    = 10:10:100;
            linear_percentage_ind  = round((linear_percentage/100)*numel(x));
            nonlinear_percentage_ind = normlogsum_penalty(linear_percentage_ind);
            nonlinear_percentage = (nonlinear_percentage_ind/numel(x))*100;
            nonlinear_percentage = round(nonlinear_percentage);
            %disp('nonlinear_percentage = ',num2str(nonlinear_percentage))
            disp_stng = ['nonlinear color code scale is: ',num2str(nonlinear_percentage)];
            %disp(disp_stng)
            disp(disp_stng)
        end
        % ----------------------------------------------------------------

        for i = 1:map_row
            for j = 1:map_col
                % if positive concentration
                if map_concentration(i,j)>0
                    
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
                    if strcmp(color_scale,'linear')
                    linear_color_number = round( min(map_concentration(i,j),max_concentration)/delta)+1;
                    col = gdm_colormap( linear_color_number, :);
                    end
                    % ----------------------------------------------------

                    % ----------------------------------------------------
                    % nonlinear color code
                    % ----------------------------------------------------
                    if strcmp(color_scale,'nonlinear')
                    linear_color_number = round( min(map_concentration(i,j),max_concentration)/delta)+1;
                    %nonlinear_color_number = round( a*log(b*linear_color_number) );
                    %nonlinear_color_number = linear_color_number;
                    nonlinear_color_number = normlogsum_penalty(linear_color_number);
                    %if nonlinear_color_number <1
                    %    nonlinear_color_number = 1;
                    %end
                    col = gdm_colormap( nonlinear_color_number, :);                
                    end
                    % ----------------------------------------------------

                    % ----- plot ------------
                    %plot(j,i,'s','Color',col,'MarkerFaceColor',col,'MarkerSize',4);
                    %plot(i,j,'s','Color',col,'MarkerFaceColor',col,'MarkerSize',4);
                    %plot(i-0.5,j-0.5,'s','Color',col,'MarkerFaceColor',col,'MarkerSize',4);
                    plot(i+0.5,j+0.5,'s','Color',col,'MarkerFaceColor',col,'MarkerSize',15);
                end
            end
        end  

end
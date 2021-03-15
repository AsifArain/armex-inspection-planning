function fHumanExpert( initialConfNum,dir_,para_ )
%=====================================================================
%   AN ASSISTANCE FOR HUMAN EXPERT TO EXPLORE THE ENVIRONMENT
%=====================================================================
% 

% conf_sequence_num = 18
% fPublishHumanReconstruction( conf_sequence_num,para_ )
% pause

%initialConfNum = 23; %13;%13;


para_.MissionStrategy = 'h-expert';
para_.FilePrefix = sprintf('%s-h',para_.FilePrefix2);

%==============================
%---- ROS
%=============================
% ros_service_conf = rossvcclient('/executed_conf_num');
ros_service_conf = rossvcclient(para_.ROSServiceExecutedConfNum);



% initialize and clear conf to execute file 
%dir_.exeConf   = '/media/a/sft/ros-catkin-ws/jfr-experiments/src/exploration_planning/plan_execution/logs/prismaforum5-04/human-exploration/';
% para_.confFile = 'plan_onestepexploration.txt';
%para_.confFile = 'robot_plan.dat';

%dir_.exeConf
%para_.confFile
%fileID = fopen([dir_.exeConf,para_.confFile],'wt'); fclose(fileID); % writing to follow

% para_.SamplingType = 'adaptive';

% %=============================
% %---- directories
% %=============================
% [ dir_ ] = fDirHumanJFR( para_ );
% dir_.Recon = dir_.Solutions;
% dir_.logs = dir_.MeasurementLogs;
% 
% 
% %-- plan type
% para_.PlanType = 'human-exploration';
% dir_.ROSLogsThis = sprintf('%s%s/',dir_.ROSLogs,para_.PlanType);


%=============================
%---- YAML file data
%=============================
% dataYAML = ReadYamlRaw([dir_.env,para_.EnvironmentYamlFile]);
% dataYAML_origin = cell2mat(dataYAML.origin);
% robot_origin_x = abs(dataYAML_origin(1)/dataYAML.resolution)
% robot_origin_y = abs(dataYAML_origin(2)/dataYAML.resolution)


% %============================
% %   DIARY
% %============================
% DateTimeNow = datestr(now,'yyyy-mm-dd-HH:MM:SS');
% diaryfilename = sprintf('diary_%s',DateTimeNow);
% diary([dir_.Solutions,diaryfilename]);



%////////////////////////////////////////////
% 
%               EXPLORATION
% 
%////////////////////////////////////////////
visualize = 0;
conf_sequence_num = initialConfNum;

while true
    % {
    %-- wait for the current conf to execute
    %------------------------------------------    
    fprintf('--> Conf# %02d is being executed...\n',conf_sequence_num);
    tSampling_i = tic;    
    conf_num_response = call(ros_service_conf);
    
    while conf_num_response.Sum ~= conf_sequence_num
        %disp('Waiting for conf to execute...')
        %fprintf('--> Last executed conf was %02d.\n',conf_num_response.Sum);
        %fprintf('--> Conf %02d is being executed..\n',conf_sequence_num);
        conf_num_response = call(ros_service_conf);
        pause(1)
    end    
    fprintf('--> Conf# %02d is executed....\n',conf_sequence_num);

    %-- sampling time
    tSampling = 1e-4*(round(toc(tSampling_i)*1e4));
    fprintf(1,'---- Sampling time: %0.4f sec \n',tSampling);

    
    %-- updated (copy) log files in the local directory
    %----------------------------------------------------
    %system('sh copyLogsFromRobot.sh');
    %system('sh copyLogsFromThisComputerLab0120171121.sh');
    %system('sh copyLogsFromThisComputerPrismaForum0120171201.sh');
    %system('sh copyLogsFromThisComputerPrismaForum0120171214.sh');    
    %filesFrom = sprint(['/media/a/sft/ros-catkin-ws/jfr-experiments/src/'...
    %                    'exploration_planning/plan_execution/logs/%s/'...
    %                    'human-exploration/'],para_.ExperimentTitle);
    %filesFrom = sprintf('%shuman-exploration/',dir_.ROSLogs);
    filesFrom = dir_.ROSLogsThis;
    filesTo = dir_.MeasurementLogs;
    %-- copy command
    syncCommand = sprintf('rsync -raz %s %s',filesFrom,filesTo);      
    %-- copy log files
    fprintf('--> Measurement files are being syncronised..\n');
    system(syncCommand);
    
        
    %-- perform reconstruction
    %--------------------------------    
    [ M,...
      mean_map,...
      cell_coords_x,...
      cell_coords_y ] = fReconstructionAdaptive( conf_sequence_num,...
                                                 para_,...
                                                 dir_,...
                                                 visualize );
                                        
    %-- publish this reconstruction
    %-----------------------------------
    %{
    fPublishReconstruction( mean_map,...
                            cell_coords_x,...
                            cell_coords_y,...
                            M,...
                            para_,...
                            dir_ );
    %}                    
    %-- updated (copy) reconstruction files in husky
    %--------------------------------
    %system('sh copyReconstructionsToRobot.sh');
    %system('sh copyReconstructionsToPublisher.sh');
    %system('sh copyReconstructionsToPublisherPrismaForum0120171214.sh');    
    %filesFrom = sprint(['/media/a/sft/ros-catkin-ws/jfr-experiments/src/'...
    %                    'exploration_planning/plan_execution/logs/%s/'...
    %                    'human-exploration/'],para_.ExperimentTitle);                    
    %filesFrom = dir_.Solutions;    
    filesFrom = dir_.Solutions;
    filesTo = dir_.ROSLogsThis;
    %filesTo = sprintf('%shuman-exploration/',dir_.ROSLogs);    
    %-- copy command
    syncCommand = sprintf('rsync -raz %s %s',filesFrom,filesTo);    
    %-- copy log files
    fprintf('--> Reconstruction is being syncronised..\n');
    system(syncCommand);
                            
    %-- update conf. sequence number
    %-----------------------------------
    conf_sequence_num = conf_sequence_num+1;
    
    %pause
    
end

% diary save
% diary off


end


function fPublishReconstruction( mean_map,...
                                 cell_coords_x,...
                                 cell_coords_y,...
                                 M,...
                                 para_,...
                                 dir_ )
% 
% 



% --- environment map data
% ================================================

file__map_env       = sprintf('%s_map_%s.dat',para_.environment,para_.EnvMapToPublish);
file__origin_env    = sprintf('%s_origin_%s.dat',para_.environment,para_.EnvMapToPublish);
file__cellsize_env  = sprintf('%s_cellsize_%s.dat',para_.environment,para_.EnvMapToPublish);

map_env      = load([dir_.env,file__map_env]);
origin_env   = load([dir_.env,file__origin_env]);
cellsize_env = load([dir_.env,file__cellsize_env]);


% --- plot environment
% ================================================

figure('name','reconstrucion on environment map'); 
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal');
colormap('gray'); caxis([-1 0.15]);

plot(origin_env(1),origin_env(2),'*r');

% --- plot reconstruction
% ================================================

map_recon = mean_map';

% ------------- select color map for gdm --------------
gdm_colormap = flipud(hot(512)); % for main fig. 

% ------ color step size --------------------------------    
max_concentration = max(map_recon(:)); % main fig.
delta = max_concentration/(size(gdm_colormap,1)-1);



for i = 1:numel(cell_coords_x)%1:numel(cell_coords_y)%map_row
    for j = 1:numel(cell_coords_y)%1:numel(cell_coords_x)%map_col
        % if positive concentration

        if map_recon(i,j) >0 %>=0

            %if map_env(i,j)>0
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
            
            
            %plot(cell_coords_x(i)+0.0,cell_coords_y(j)+0.0,...
            plot( (cell_coords_x(i)*(para_.ReconTomographyCellSizeM/cellsize_env) ) +0.0,...
                  (cell_coords_y(j)*(para_.ReconTomographyCellSizeM/cellsize_env) ) +0.0,...
                 's',...
                 'MarkerEdgeColor','k',... 'MarkerEdgeColor',col,...
                 'MarkerFaceColor',col,...
                 'MarkerSize',5); %5
            %end
        end
    end
end



%--- sensing positions
%----------------------------------
sensingPositions = [unique(M(:,1)),unique(M(:,2))];
plot(sensingPositions(:,1),sensingPositions(:,2),'og');
% plot(sensingPositions(end,1),sensingPositions(end,2),'ok','MarkerFaceColor','k');
plot(sensingPositions(end,1),sensingPositions(end,2),'*g');







end


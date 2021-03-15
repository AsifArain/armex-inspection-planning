function fRobotMeasurement( sensing_confs,...
                            conf_sequence_num,...
                            dir_,...
                            para_ )
%
% EXECUTE THE CONFIGURATION USING A REAL PLATFORM


ros_service_conf = rossvcclient(para_.ROSServiceExecutedConfNum);


%-- create plan file, if not already exists
if exist([dir_.ROSLogsThis,para_.RobotGlobalPlanFileName],'file') == 0    
    fileID = fopen([dir_.ROSLogsThis,para_.RobotGlobalPlanFileName],'wt'); 
    fclose(fileID);
end

%///////////////////////////////////////////////////////////////
%
%   Write sensing conf to the file for robot
%
%///////////////////////////////////////////////////////////////

switch para_.PlanType
    
        %============================================================================
    case {'one-step-exploration','two-step-exploration-tomography-adaptive'}
        %============================================================================
        
        conf_num_response = call(ros_service_conf);
        if conf_sequence_num > conf_num_response.Sum

            dlmwrite([dir_.ROSLogsThis,para_.RobotGlobalPlanFileName],...
                     [sensing_confs(:,1),... % x-position
                      sensing_confs(:,2),... % y-position
                      sensing_confs(:,3),... % orientation
                      ],...               
                      '-append','delimiter',' ');
            fileID = fopen([dir_.ROSLogsThis,para_.RobotGlobalPlanFileName],'a'); 
            fclose(fileID);
        end
        
        %============================================================================
    case {'two-step-exploration-detection','two-step-exploration-tomography-fixed'}
        %============================================================================
        
        dlmwrite([dir_.ROSLogsThis,para_.RobotGlobalPlanFileName],...
                 [sensing_confs(:,1),... % x-position
                  sensing_confs(:,2),... % y-position
                  sensing_confs(:,3),... % orientation
                  ],...               
                  'delimiter',' ');
        fileID = fopen([dir_.ROSLogsThis,para_.RobotGlobalPlanFileName],'a'); 
        fclose(fileID);
        
end
    

%-- copy command
syncCommand = sprintf('rsync -raz %s %s',dir_.ROSLogsThis,dir_.logs);


% fprintf('--> Conf# %02d is being executed...\n',conf_sequence_num);
fprintf('--> Till conf #%02d to execute...\n',conf_sequence_num);

conf_num_response = call(ros_service_conf);
fprintf('--> Last conf #%02d was executed....\n',conf_num_response.Sum);

prevExecutedConf = conf_num_response.Sum;

tSampling_i = tic;

conf_num_response = call(ros_service_conf);
% while conf_num_response.Sum ~= conf_sequence_num
while conf_num_response.Sum < conf_sequence_num
    
    %disp('Waiting for conf to execute...')
    %fprintf('--> Last executed conf was %02d.\n',conf_num_response.Sum);
    %fprintf('--> Conf %02d is being executed..\n',conf_sequence_num);
    conf_num_response = call(ros_service_conf);
    %conf_num_response.Sum
    %conf_sequence_num
    pause(1)
    
    %-- copy log files
    %system(syncCommand);
    
    if conf_num_response.Sum > prevExecutedConf
        fprintf('--> Now conf #%02d is executed....\n',conf_num_response.Sum);
        prevExecutedConf = conf_num_response.Sum;
    end
    
end
% fprintf('--> Conf# %02d is executed....\n',conf_sequence_num);

%-- sampling time
tSampling = 1e-4*(round(toc(tSampling_i)*1e4));
fprintf(1,'---- Sampling time: %0.4f sec \n',tSampling);


%-- copy log files
fprintf('--> Measurement files are being syncronised..\n');
system(syncCommand);
    


%-------------
%--- robot logs to M for this conf.
%====================================================
% % fprintf('---> Generating measurement matrix (M) for the last conf.\n');
% fprintf('---> Generating M matrix for conf# %02d.\n',conf_sequence_num);
% measure_file = sprintf('measurements_conf%d.dat',conf_sequence_num);
% 
% visualize = 0;
% [ M_m,M_cell ] = fRobotLogs2M( measure_file,...
%                                para_,...
%                                dir_,...
%                                visualize );
% 
% Mfile = sprintf('M_conf%02d.mat',conf_sequence_num);
% save([dir_.MeasurementLogs,Mfile],'M_m','M_cell');



%///////////////////////////////////////////////////////////////
%
%   Robot logs to M for the executed confs.
%
%//////////////////////////////////////////////////////////////
switch para_.PlanType
    
        %============================================================================
    case {'one-step-exploration','two-step-exploration-tomography-adaptive'}
        %============================================================================
        
        fprintf('---> Generating M matrix for conf #%02d.\n',conf_sequence_num);
        measure_file = sprintf('measurements_conf%d.dat',conf_sequence_num);

        visualize = 0;
        [ M_m,M_cell ] = fRobotLogs2M( measure_file,...
                                       para_,...
                                       dir_,...
                                       visualize );

        Mfile = sprintf('M_conf%02d.mat',conf_sequence_num);
        %save([dir_.MeasurementLogs,Mfile],'M_m','M_cell');
        save([dir_.logs,Mfile],'M_m','M_cell');
        
        %============================================================================
    case {'two-step-exploration-detection','two-step-exploration-tomography-fixed'}
        %============================================================================
        
        fprintf('---> Generating M matrix for all confs.\n');
        
        for this_num = 1:conf_sequence_num
            
            fprintf('---- M matrix for conf #%02d.\n',this_num);
            measure_file = sprintf('measurements_conf%d.dat',this_num);

            visualize = 0;
            [ M_m,M_cell ] = fRobotLogs2M( measure_file,...
                                           para_,...
                                           dir_,...
                                           visualize );

            Mfile = sprintf('M_conf%02d.mat',this_num);
            %save([dir_.MeasurementLogs,Mfile],'M_m','M_cell');
            save([dir_.logs,Mfile],'M_m','M_cell');
        end
        
end
    


end
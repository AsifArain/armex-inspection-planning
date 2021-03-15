function [ cellsToBeScanned,...
           cellsScanned ] = fExecuteAllConfs( sensing_confs,...
                                              conf_sequence_num,...
                                              map_env,...
                                              FoV,...
                                              SensorRange,...
                                              cellsToBeScanned,...
                                              para_,...
                                              dir_ )
%
% EXECUTE THE CONFIGURATION AND UPDATE THE LIST OF CELLS TO BE SCANNED


fprintf('\n    ****************************************************\n');
fprintf('             SAMPLING FOR CONFIGURATION #%02d\n',conf_sequence_num);
fprintf('    ****************************************************\n');


switch para_.SensingSystem


        %======================================
    case 'robot'
        %======================================

        sensing_confs = sensing_confs(:,1:3);
        
        
        %-- create plan file, if not already exists
        %if exist([dir_.ROSLogsThis,para_.RobotGlobalPlanFileName],'file') == 0    
            fileID = fopen([dir_.ROSLogsThis,para_.RobotGlobalPlanFileName],'wt'); 
            fclose(fileID);
        %end
        
        % -- write sensing conf to the file for robot
        dlmwrite([dir_.ROSLogsThis,para_.RobotGlobalPlanFileName],...
                 [sensing_confs(:,1),... % x-position
                  sensing_confs(:,2),... % y-position
                  sensing_confs(:,3),... % orientation
                  ],...
                  'delimiter',' ');
        fileID = fopen([dir_.ROSLogsThis,para_.RobotGlobalPlanFileName],'a'); 
        fclose(fileID);
        

        % dir_.ROSLogsThis

        %-- copy command
        syncCommand = sprintf('rsync -raz %s %s',dir_.ROSLogsThis,dir_.logs);

        
        
        
        ros_service_conf = rossvcclient(para_.ROSServiceExecutedConfNum);


        fRobotMeasurement( sensing_confs,...
                           conf_sequence_num,...
                           ros_service_conf,...
                           dir_,...
                           para_ );


                       
        % -- cells scanned
        %====================================================
        %TODO
        warning('To be fixed according to the robot scan, and not simulation');
        visualize = 0;
        fprintf('--> Coverage for conf %02d\n',conf_sequence_num);
        [ cellsScanned ] = fSingleConfCoverage( sensing_confs,...
                                                map_env,...
                                                FoV,...
                                                SensorRange,...
                                                para_,...
                                                visualize );

        % -- update the global list of cells to be scanned
        %====================================================
        scannedInd = ismember(cellsToBeScanned,cellsScanned);
        cellsToBeScanned(scannedInd) = [];
                       
                       
        %======================================
    case 'simulation'
        %======================================

        fprintf('--> Measurments for conf %02d\n',conf_sequence_num);
        % -- output measurement file name
        measu_file = sprintf('measurement_conf%02d.dat',conf_sequence_num);
        visualize = 0;

        % -- if output file does not exist
        %fname_out = [dir_.meas{time_stamp_number},measu_file];
        %if exist(fname_out, 'file') == 0
            % -- collect integral measurements
            %fprintf('Remote Sensor:\n');
            %fprintf('Set# %02d, %s:\n',i, conf_file);
            fSimulatedMeasurement( sensing_confs,...
                                   measu_file,...
                                   map_env,...
                                   para_,...
                                   dir_,...
                                   visualize );
        %end
        
        
        %--- save executed conf to file
        %====================================================
        ConfFile = 'executed_confs.txt';
        %-- create plan file, if not already exists
        if exist([dir_.Solutions,ConfFile],'file') == 0
            fileID = fopen([dir_.Solutions,ConfFile],'wt'); 
            fclose(fileID);
        end
        % -- write integral concentration data to the file
        dlmwrite([dir_.Solutions,ConfFile],sensing_confs,'-append','delimiter',' ');
        fileID = fopen([dir_.Solutions,ConfFile],'a'); fclose(fileID);

        
        
                               
        % -- cells scanned
        %====================================================
        fprintf('--> Coverage for conf %02d\n',conf_sequence_num);
        [ cellsScanned ] = fSingleConfCoverage( sensing_confs,...
                                                map_env,...
                                                FoV,...
                                                SensorRange,...
                                                para_,...
                                                visualize );

        % -- update the global list of cells to be scanned
        %====================================================
        scannedInd = ismember(cellsToBeScanned,cellsScanned);
        cellsToBeScanned(scannedInd) = [];


end






end

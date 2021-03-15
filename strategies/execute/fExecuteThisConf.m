function [ cellsToBeScanned,...
           cellsScanned ] = fExecuteThisConf( sensing_conf,...
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

% fprintf(['    <x,y,theta,fov,r,s,t> \n'...
%          '    <%02.2f,%02.2f,%02.2f,%02.2f,%02.2f,%02.2f,%02.2f>\n'],sensing_conf);
fprintf(['    <x,y,theta,fov,r,s,t,for> \n'...
         '    <%02.2f,%02.2f,%02.2f,%02.2f,%02.2f,%02.2f,%02.2f,%02.2f>\n'],sensing_conf);

switch para_.SensingSystem


        %======================================
    case 'robot'
        %======================================

        sensing_conf = sensing_conf(:,1:3);        

        fRobotMeasurement( sensing_conf,...
                           conf_sequence_num,...
                           dir_,...
                           para_ );


                       
        % -- cells scanned
        %====================================================
        %TODO
        %warning('To be fixed according to the robot scan, and not simulation');
        visualize = 0;
        fprintf('--> Coverage for conf %02d\n',conf_sequence_num);
        [ cellsScanned ] = fSingleConfCoverage( sensing_conf,...
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
    case 'simulated-sampling'
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
            fSimulatedMeasurement( sensing_conf,...
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
        if exist([dir_.str,ConfFile],'file') == 0
            fileID = fopen([dir_.str,ConfFile],'wt'); 
            fclose(fileID);
        end
        % -- write integral concentration data to the file
        dlmwrite([dir_.str,ConfFile],sensing_conf,...
            '-append','delimiter',' ','precision','%06.2f');
        fileID = fopen([dir_.str,ConfFile],'a'); fclose(fileID);

        
        
                               
        % -- cells scanned
        %====================================================
        fprintf('--> Coverage for conf %02d\n',conf_sequence_num);
        [ cellsScanned ] = fSingleConfCoverage( sensing_conf,...
                                                map_env,...
                                                FoV,...
                                                SensorRange );

        % -- update the global list of cells to be scanned
        %====================================================
        scannedInd = ismember(cellsToBeScanned,cellsScanned);
        cellsToBeScanned(scannedInd) = [];


end






end

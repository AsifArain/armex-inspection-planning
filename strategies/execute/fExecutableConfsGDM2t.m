function [sensing_confs] = fExecutableConfsGDM2t( selectedConfFusion_ioh,...
                                                  confSequenceNum,...
                                                  SensorRange,...
                                                  map_env,...
                                                  FoV,...
                                                  dir_,...
                                                  para_ )
%
% INITIAL COVERAGE PLAN FOR GAS DETECTION



% -- add sequence number to the fused conf
selectedConfPlanned_sioh = [confSequenceNum,selectedConfFusion_ioh];
% -- sorted by sequence number
selectedConfPlanned_sioh = sortrows(selectedConfPlanned_sioh,1);

% -- sensing conf to execute
[confs_r,confs_c] = ind2sub(size(map_env),selectedConfPlanned_sioh(:,2));
confs_o = selectedConfPlanned_sioh(:,3);


switch para_.SensingSystem
    
    case 'robot-sampling'        
        sensing_confs = [confs_r+0.0,confs_c+0.0,confs_o];
                
        %
        %-- create plan file, if not already exists
        if exist([dir_.ROSLogsThis,para_.RobotLocalPlanFileName],'file') == 0    
            fileID = fopen([dir_.ROSLogsThis,para_.RobotLocalPlanFileName],'wt'); 
            fclose(fileID);
        end
        %-- write local/current plan
        dlmwrite([dir_.ROSLogsThis,para_.RobotLocalPlanFileName],...
                 [sensing_confs(:,1),... % x-position
                  sensing_confs(:,2),... % y-position
                  sensing_confs(:,3),... % orientation
                  ],...               
                  'delimiter',' ');
        fileID = fopen([dir_.ROSLogsThis,para_.RobotLocalPlanFileName],'a'); 
        fclose(fileID);
        
        %-- copy command
        thisFile = ([dir_.ROSLogsThis,para_.RobotLocalPlanFileName]);
        thatFile = ([dir_.logs,para_.RobotLocalPlanFileName]);
        syncCommand = sprintf('rsync -raz %s %s',thisFile,thatFile);
        system(syncCommand);
        
        
        
    case 'simulated-sampling'
        sensing_confs = [confs_r+0.5,confs_c+0.5,confs_o];
end


% -- other sensing parameters
conf_fov = FoV;
conf_r   = SensorRange.cell;
conf_frq = FoV*2;
conf_tim = para_.GroundTruthTimeStamp;

sensing_confs(:,4) = conf_fov*ones(size(sensing_confs,1),1);
sensing_confs(:,5) = conf_r*  ones(size(sensing_confs,1),1);
sensing_confs(:,6) = conf_frq*ones(size(sensing_confs,1),1);
sensing_confs(:,7) = conf_tim*ones(size(sensing_confs,1),1);



end

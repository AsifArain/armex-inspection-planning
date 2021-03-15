function [sensing_confs] = fExecutableConfsGD( FoVfaces,...
                                               SensorRange,...
                                               FoV,...
                                               oV,...
                                               C,...
                                               confSequenceNums,...
                                               confKept,...
                                               map_env,...
                                               dir_,...
                                               para_ )
%
% INITIAL COVERAGE PLAN FOR GAS DETECTION

f = FoVfaces.num;
Conf = zeros(size(oV,2),1);
Conf(confKept) = C;
Conf_ind = find(Conf);                      % IDs of selected conf.
Conf_cell = fix(Conf_ind/f)+(~~mod(Conf_ind,f));   % cell num for each conf num.
[Conf_r,Conf_c] = ind2sub(size(map_env),Conf_cell); % row and col of cell num.
Conf_o = FoVfaces.lgt(mod(Conf_ind,f)+(~mod(Conf_ind,f)*f));  % conf num within a cell.

% -- pose
% sensing_confs = [Conf_r+0.5,Conf_c+0.5,Conf_o];

switch para_.SensingSystem
    
    case 'robot-sampling'        
        sensing_confs = [Conf_r+0.0,Conf_c+0.0,Conf_o];
                
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
        sensing_confs = [Conf_r+0.5,Conf_c+0.5,Conf_o];
end




% -- other parameters
conf_fov = FoV;
conf_r   = SensorRange.cell;
conf_frq = FoV*2;
conf_tim = para_.GroundTruthTimeStamp;

sensing_confs(:,4) = conf_fov*ones(size(sensing_confs,1),1);
sensing_confs(:,5) = conf_r*  ones(size(sensing_confs,1),1);
sensing_confs(:,6) = conf_frq*ones(size(sensing_confs,1),1);
sensing_confs(:,7) = conf_tim*ones(size(sensing_confs,1),1);


% -- conf by sequence
sensing_confs = [confSequenceNums,sensing_confs];
sensing_confs = sortrows(sensing_confs,1);
sensing_confs = sensing_confs(:,2:end);



 




end

close all; clear all; clc;



para_.ExperimentTitle = 'prismaforum5-10-publish';
% para_.Strategy = 'one-step-exploration'; 
% para_.Strategy = 'two-step-exploration-detection'; 
para_.Strategy = 'two-step-exploration-tomography-fixed'; 

dir_.env = 'environments/';
file_map_env = 'prismaforum5_map_coverage.dat';

map_env = load([dir_.env,file_map_env]);


dir_.ROSLogsThis = sprintf(['/media/a/sft/ros-catkin-ws/jfr-experiments/src/'...
    'exploration_planning/plan_execution/logs/%s/%s/'],para_.ExperimentTitle,para_.Strategy);

% confs_filename = 'planned_confs_global.dat';
% path_filename = 'planned_path_global.dat';


% sensing_confs = load([dir_.ROSLogsThis,confs_filename]);
% sensing_confs = sensing_confs(:,1:2);

para_.ConfType = 'Planned';


%[ travelingPath,travelingDist ] = fTravelingPathDist( sensing_confs,map_env );
[ travelingPath,travelingDist ] = fTravelingPathDist( map_env,para_,dir_ );





% fileID = fopen([dir_.ROSLogsThis,path_filename],'wt'); fclose(fileID);
% % -- write integral concentration data to the file
% dlmwrite([dir_.ROSLogsThis,path_filename],travelingPath,'delimiter',' ');
% fileID = fopen([dir_.ROSLogsThis,path_filename],'a'); fclose(fileID);






%-- traveling dist
%------------------------------
% translational_steps = abs(diff(travelingPath));
% steps_xy = sum(translational_steps,2);
% 
% straight_step_size = 1*0.5;
% diagonal_step_size = 1.4142*0.5;
% 
% straight_dist = numel(find(steps_xy==1))*straight_step_size;
% diagonal_dist = numel(find(steps_xy==2))*diagonal_step_size;
% 
% 
% travelingDist = straight_dist+diagonal_dist

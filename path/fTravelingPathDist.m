function [ travelingPath,travelingDist ] = fTravelingPathDist( map_env,para_,dir_ )
%
%*************************************************************************
% Traveling path between a list of sensing configurations.
%*************************************************************************
%
% @AsifArain - 07-Apr-2018


file__cell_env = sprintf('%s_cellsize_coverage.dat',para_.environment);
cellsize_env   = load([dir_.env,file__cell_env]);

%******************************************************
%-- initialize
%******************************************************
time_i = tic;
travelingPath = [];

para_.ConfType

%******************************************************
%   SENSING CONFIGURATIONS
%******************************************************
switch para_.ConfType
    
    case 'Planned'
        
        %-- real experiments
        %{
        confFile = 'planned_confs_global.dat';
        sensing_confs = load([dir_.ROSLogsThis,confFile]);
        %}
        
        %-- simulation experiments
        confFile = 'executed_confs_final.txt';
        sensing_confs = load([dir_.str,confFile]);
        sensing_confs = floor(sensing_confs);
        
        
        sensing_confs = sensing_confs(:,1:2);
        %sensing_confs(end+1,:) = sensing_confs(1,:);
        
    case 'Executed'
        
        confFile = 'executed_confs.dat';
        sensing_confs = load([dir_.ROSLogsThis,confFile]);
        
        file__org_conf = sprintf('%s_origin_conf.dat',para_.environment);
        origin_conf = load([dir_.env,file__org_conf]);
        
        file__cell_conf = sprintf('%s_cellsize_conf.dat',para_.environment);
        cellsize_conf   = load([dir_.env,file__cell_conf]);
        
        sensing_confs(:,1) = (sensing_confs(:,1)/cellsize_conf+origin_conf(1)-0.5+1);
        sensing_confs(:,2) = (sensing_confs(:,2)/cellsize_conf+origin_conf(2)-0.5+1);
        
        sensing_confs = round(sensing_confs);
        sensing_confs = sensing_confs(:,1:2);
        %sensing_confs(end+1,:) = sensing_confs(1,:);
        
end




switch para_.MissionStrategy 
    case {'one-step-exploration-jfr','two-step-exploration-jfr','1t-armex','2t-armex'}
        initial_position = [53,25];
        %initial_position = [26 19] %-- simulation experiment sample-env-paper-04
    case 'human-expert'
        initial_position = [212,20];
end


sensing_confs = [initial_position; sensing_confs; initial_position];


sensing_confs

size(sensing_confs)


%----- FOR INSITU -- temporary
% numel(find(map_env))
% 
% [conf1,conf2] = ind2sub(size(map_env),find(map_env));
% sensing_confs = [[conf1(1) conf2(2)];[conf1(end) conf2(end)]]
%-------------------------------


%******************************************************
%-- find nearst conf interatively
%******************************************************
for i = 1:size(sensing_confs,1)-1
    
    
    %-----------------
    confFrom = sensing_confs(i,:);
    confTo   = sensing_confs(i+1,:);
    
    [ OptimalPath ] = fPath2Confs(confFrom,confTo,map_env );
    
    travelingPath = [travelingPath;flipud(OptimalPath)];
    
end
   
%******************************************************
% prepare final list
%******************************************************
% travelingPath = unique(travelingPath,'rows');
    

%******************************************************
% traveling dist
%******************************************************
translational_steps = abs(diff(travelingPath));
steps_xy = sum(translational_steps,2);

straight_step_size = 1*cellsize_env;
diagonal_step_size = 1.4142*cellsize_env;

straight_dist = numel(find(steps_xy==1))*straight_step_size;
diagonal_dist = numel(find(steps_xy==2))*diagonal_step_size;

travelingDist = straight_dist+diagonal_dist;

fprintf(1,'-- traveling dist: %0.4f m \n',travelingDist)


pause
%******************************************************
% save path and dist (data files) in the log dir
%******************************************************
path_filename = 'planned_path_global.dat';
dist_filename = 'planned_dist_global.dat';

switch para_.SensingSystem
    case 'robot'
        
        fileID = fopen([dir_.ROSLogsThis,path_filename],'wt'); fclose(fileID);
        dlmwrite([dir_.ROSLogsThis,path_filename],travelingPath,'delimiter',' ');
        fileID = fopen([dir_.ROSLogsThis,path_filename],'a'); fclose(fileID);

        fileID = fopen([dir_.ROSLogsThis,dist_filename],'wt'); fclose(fileID);
        dlmwrite([dir_.ROSLogsThis,dist_filename],travelingDist,'delimiter',' ');
        fileID = fopen([dir_.ROSLogsThis,dist_filename],'a'); fclose(fileID);
        
    case 'simulation'
        
        fileID = fopen([dir_.str,path_filename],'wt'); fclose(fileID);
        dlmwrite([dir_.str,path_filename],travelingPath,'delimiter',' ');
        fileID = fopen([dir_.str,path_filename],'a'); fclose(fileID);

        fileID = fopen([dir_.str,dist_filename],'wt'); fclose(fileID);
        dlmwrite([dir_.str,dist_filename],travelingDist,'delimiter',' ');
        fileID = fopen([dir_.str,dist_filename],'a'); fclose(fileID);
end

%******************************************************
% Total computation time.
%******************************************************
tComputation = 1e-4*(round(toc(time_i)*1e4));
fprintf(1,'-- computation time: %0.4f sec \n',tComputation)

end


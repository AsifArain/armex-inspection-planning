function [ dir_ ] = fDir2tARMEx( para_ )
%
%

fprintf(1,['\n\n',...
    '------------------------------------------\n',...
    '  Setting directory parameters - 2t-ARMEx \n',...
    '------------------------------------------\n'])

dir_.env = ['environments/'];
dir_.gt  = ['results/',para_.ExperimentTitle,'/'];
% dir_.str = sprintf('results/%s/%s/',...
%                    para_.ExperimentTitle,...
%                    para_.MissionStrategy);
              
dir_.str = sprintf('results/%s/%s/%s/',...
                   para_.ExperimentTitle,...
                   para_.MissionStrategy,...
                   para_.SensingSystem);
              
dir_.diary = sprintf('%s/diary/',dir_.str);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%   GAS DETECTION PLANNING DIRECTORIES
%
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               
dir_.gdplan       = sprintf('%s/detection/',dir_.str);
dir_.gdplan_spp   = sprintf('%s/spp/',dir_.gdplan);
dir_.gdplan_logs  = sprintf('%s/measurements/%s',dir_.gdplan);
dir_.gdplan_recon = sprintf('%s/recon/%s',dir_.gdplan);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%   GAS DISTRIBUTION MAPPING PLANNING DIRECTORIES
%
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% dir_gdmParameters = sprintf('%s/gdm/alpha%03d_beta%03d_gamma%03d',...
%                    dir_str,...
%                    para_.alpha*100,...
%                    para_.beta*100,...
%                    para_.gamma*100);
% dir_.gdmplan          = sprintf('%s/',dir_gdmParameters);

dir_.gdmplan          = sprintf('%s/gdm/',dir_.str);
% dir_.gdmplan          = sprintf('%s/gdm-withhotspotfusion/',dir_.str);
% dir_.gdmplan          = sprintf('%sgdm-withouthotspotfusion/',dir_.str);
dir_.gdmplan_hotspots = sprintf('%shotspots/',dir_.gdmplan);
dir_.gdmplan_spp      = sprintf('%sspp/',dir_.gdmplan);
dir_.gdmplan_logs     = sprintf('%smeasurements/',dir_.gdmplan);
dir_.gdmplan_recon    = sprintf('%srecon/',dir_.gdmplan);
dir_.gdmplan_evaluation      = sprintf('%sevaluation/',dir_.gdmplan);


ros_pkg = '/dev/cod/ros-catkin-ws/armex/src/mission_strategies/plan_execution/';
dir_.RobotPlan = sprintf('%s',ros_pkg);
dir_.ROSLogs = sprintf('%slogs/%s/',ros_pkg,para_.ExperimentTitle);
dir_.ROSLogsThis = sprintf('%s%s/',dir_.ROSLogs,para_.MissionStrategy);


%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%   CREATE DIRECTORIES
%
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if ~exist(dir_.env,'dir')
    mkdir(dir_.env);
end

if ~exist(dir_.gt,'dir')
    mkdir(dir_.gt);
end

if ~exist(dir_.diary,'dir')
    mkdir(dir_.diary);
end

if ~exist(dir_.gdplan,'dir')
    mkdir(dir_.gdplan);
end

if ~exist(dir_.gdplan_spp,'dir')
    mkdir(dir_.gdplan_spp);
end

if ~exist(dir_.gdplan_logs,'dir')
    mkdir(dir_.gdplan_logs);
end

if ~exist(dir_.gdplan_recon,'dir')
    mkdir(dir_.gdplan_recon);
end

if ~exist(dir_.gdmplan,'dir')
    mkdir(dir_.gdmplan);
end

if ~exist(dir_.gdmplan_hotspots,'dir')
    mkdir(dir_.gdmplan_hotspots);
end

if ~exist(dir_.gdmplan_spp,'dir')
    mkdir(dir_.gdmplan_spp);
end

if ~exist(dir_.gdmplan_logs,'dir')
    mkdir(dir_.gdmplan_logs);
end

if ~exist(dir_.gdmplan_recon,'dir')
    mkdir(dir_.gdmplan_recon);
end

if ~exist(dir_.gdmplan_evaluation,'dir')
    mkdir(dir_.gdmplan_evaluation);
end


end
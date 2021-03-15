function [ dir_ ] = fDir1tARMEx( para_ )
%
%


fprintf(1,['\n\n',...
    '------------------------------------------\n',...
    '  Setting directory parameters - 1t-ARMEx \n',...
    '------------------------------------------\n'])

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%   DIRECTORIES
%
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

dir_.env          = ['environments/'];
dir_.gt           = ['results/',para_.ExperimentTitle,'/'];
% dir_.str          = sprintf('results/%s/%s/',para_.ExperimentTitle,para_.MissionStrategy);
dir_.str          = sprintf('results/%s/%s/%s/',para_.ExperimentTitle,para_.MissionStrategy,para_.SensingSystem);
% dir_.str          = sprintf('results/%s/%s/coverage%03d/',...
%     para_.ExperimentTitle,para_.MissionStrategy,para_.minSensingCoverage1t*100);
% dir_.str          = sprintf('results/%s/%s/ees/',para_.ExperimentTitle,para_.MissionStrategy);
%dir_.str          = sprintf('results/sensing_coverage_comparison/%s/%s/coverage%03d/',...
%    para_.ExperimentTitle,para_.MissionStrategy,para_.minSensingCoverage1t*100);

dir_.logs         = sprintf('%smeasurements/',dir_.str);
dir_.recon        = sprintf('%srecon/',dir_.str);
dir_.hotspots     = sprintf('%shotspots/',dir_.str);
dir_.evaluation   = sprintf('%sevaluation/',dir_.str);
dir_.sppdetection = sprintf('%sspp-detection/',dir_.str);
dir_.sppgdm       = sprintf('%sspp-gdm/',dir_.str);
dir_.diary        = sprintf('%sdiary/',dir_.str);


ros_pkg = '/media/a/sft/ros-catkin-ws/jfr-experiments/src/exploration_planning/plan_execution/';
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

if ~exist(dir_.str,'dir')
    mkdir(dir_.str);
end

if ~exist(dir_.logs,'dir')
    mkdir(dir_.logs);
end

if ~exist(dir_.recon,'dir')
    mkdir(dir_.recon);
end

if ~exist(dir_.hotspots,'dir')
    mkdir(dir_.hotspots);
end

if ~exist(dir_.evaluation,'dir')
    mkdir(dir_.evaluation);
end

if ~exist(dir_.sppdetection,'dir')
    mkdir(dir_.sppdetection);
end

if ~exist(dir_.sppgdm,'dir')
    mkdir(dir_.sppgdm);
end

if ~exist(dir_.diary,'dir')
    mkdir(dir_.diary);    
end

end
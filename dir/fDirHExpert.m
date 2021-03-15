function [ dir_ ] = fDirHumanJFR( para_ )
%
%



dir_.env = 'environments/';
dir_.GroundTruthLogs = ['results/',para_.ExperimentTitle,'/'];



dir_explo = sprintf('results/%s/%s',...
                   para_.ExperimentTitle,...
                   para_.ExplorationStrategy);

dir_tom2 = dir_explo;
% sprintf('%s/alpha%03d_beta%03d_gamma%03d',...
%                    dir_explo,...
%                    alpha*100,...
%                    beta*100,...
%                    gamma*100);

% dir_.MeasurementLogs = sprintf('%s/measurement-logs/',dir_tom2);
dir_.MeasurementLogs = sprintf('%s/measurement-logs/%s/',dir_tom2,para_.SensingSystem);
dir_.Solutions = sprintf('%s/solutions/',dir_tom2);
dir_.Evaluation = sprintf('%s/evaluation/',dir_tom2);


ros_pkg = '/media/a/sft/ros-catkin-ws/jfr-experiments/src/exploration_planning/plan_execution/';
dir_.RobotPlan = sprintf('%s',ros_pkg);
dir_.ROSLogs = sprintf('%slogs/%s/',ros_pkg,para_.ExperimentTitle);



if ~exist(dir_.env,'dir')
    mkdir(dir_.env);
end

if ~exist(dir_.GroundTruthLogs,'dir')
    mkdir(dir_.GroundTruthLogs);
end


if ~exist(dir_.MeasurementLogs,'dir')
    mkdir(dir_.MeasurementLogs);
end


if ~exist(dir_.Solutions,'dir')
    mkdir(dir_.Solutions);
end

if ~exist(dir_.Evaluation,'dir')
    mkdir(dir_.Evaluation);
end

%----------------------------------------------
% add path
% ---------------------------------------------
addpath(genpath([pwd,'/preprocess']))
addpath(genpath([pwd,'/integration']))
addpath(genpath([pwd,'/tsp']))
addpath(genpath([pwd,'/publish']))
addpath(genpath([pwd,'/template_matching']))
addpath(genpath([pwd,'/detection']))
addpath(genpath([pwd,'/robot-plan']))
addpath(genpath([pwd,'/experts-solution']))
addpath(genpath([pwd,'/measurements']))
addpath(genpath([pwd,'/reconstruction']))
addpath(genpath([pwd,'/dispersion']))
addpath(genpath([pwd,'/evaluation']))
addpath(genpath([pwd,'/fusion']))

end
function [ para_ ] = fAlgorithmParams()

%
%


fprintf(1,['\n\n',...
    '------------------------------------------\n',...
    '  Setting algorithm parameters           \n',...
    '------------------------------------------\n'])

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%       EXPERIMENT TYPE (SIMULATION OR REAL)
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% % para_.SensingSystem = 'robot';
% para_.SensingSystem = 'simulation';
% para_.SensingSystem = 'robot-sampling';
para_.SensingSystem = 'simulated-sampling';


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%       ENVIRONMENT MAP TO PUBLISH
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
para_.EnvMapToPublish = 'raycasting';
% para_.EnvMapToPublish = 'coverage';
% para_.EnvMapToPublish = 'navigation';
% para_.EnvMapToPublish = 'conf';


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%       RESULT FILE NAMING
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% para_.RobotPlanFileName = 'robot_plan.dat';
% para_.RobotPlanFileName = 'planned_confs.dat';
para_.RobotGlobalPlanFileName = 'planned_confs_global.dat';
para_.RobotLocalPlanFileName = 'planned_confs_local.dat';
para_.RobotHcsFileName = 'hotspots.dat';




%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%       SENSING COVERAGE PARAMETER
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% para_.EffectiveToFeasibleCoverageIndRatioOfConf = 0.15;
para_.ExclusiveToTotalCoverageThreshold = 0.15;

para_.minSensingCoverage1t = 0.95 %0.95; %0.85;
para_.minSensingCoverage2t = 1.00; %0.85;



%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%       MAXIMUM ERQs (SIMULATION RESULTS)
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
para_.MapQuality2Conf = 0.6265;
para_.MapQuality3Conf = 0.7847;
para_.MapQuality4Conf = 0.8502;
para_.MapQuality5Conf = 0.8956;

para_.MapQualityRelative2Conf = 0.6995;
para_.MapQualityRelative3Conf = 0.8762;
para_.MapQualityRelative4Conf = 0.9493;
para_.MapQualityRelative5Conf = 1.0000;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%       DESIRED ERQs
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% para_.DesiredMapQuality = 0.50; %0.55; %0.65;
para_.desiredERQ1t = 0.60; %0.50; %0.55; %0.65;
%para_.desiredERQ2t = 0.55; % JFR experiments
%para_.desiredERQ2t = 0.75; % for example reconstruction
% para_.desiredERQ2t = 0.50; % for evaluation of with and without hotspt fusion
para_.desiredERQ2t = 0.65; % for evaluation of with and without hotspt fusion

para_.DesiredMapQualityRelative = 0.75;


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%       CROSS ANGLE ERGs (SIMULATION RESULTS)
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
para_.xGainMeu2Conf_DEG = mean([42.5,75])*[1];
para_.xGainMeu3Conf_DEG = mean([52.5,55])*[1,2];
para_.xGainMeu4Conf_DEG = mean([37.5,40])*[1,2,3];
para_.xGainMeu5Conf_DEG = mean([37.5,40])*[1,2,3,4];

para_.xGainSig2Conf_DEG = 20;
para_.xGainSig3Conf_DEG = 15;
para_.xGainSig4Conf_DEG = 10;
para_.xGainSig5Conf_DEG = 10;



%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%       MULTI-CRITERIA PARAMETER VALUES
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
para_.alpha = 0.95;
para_.beta  = 0.70;
para_.gamma = 0.20;


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%       MAXIMUM CANDIDATE CONFS FOR GDM PLANNING
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

para_.maxCandidateConf_Hotspot_1t = 30;%200; 
para_.deltaCandConf_Hotspot_1t = 5;

para_.maxCandidateConf_Fusion_1t = 15;%200; 
para_.deltaCandConf_Fusion_1t = 5;

para_.maxCandidateConf_Hotspot_2t = 100;%200; 
para_.deltaCandConf_Hotspot_2t = 15;

para_.maxCandidateConf_Fusion_2t = 100 %80; %200 %80; %50;%200; % for test
para_.deltaCandConf_Fusion_2t = 15;

% --------------------------------------------------------


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%       HOTSPOTS ESTIMATION
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Circular radius (m) to look for free cell, if hc is an occupied cell.
para_.radiusFreeCell_m = 10; %5; %2; 

%- discarding mean points - weight threshold
para_.HcsCutoffWeightFactor = 0.5;
%- discarding mean points - concentration threshold
para_.HcsCutoffPPMM = 25000; % sample-env-paper


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%       MAXIMUM CANDIDATE CONFS FOR GDM PLANNING
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% para_.maxConfLocalSolution1t = 2;%3;%4;
% para_.maxConfLocalSolution2t = 3;%4;

para_.maxConfToSelect_Hotspot_1t = 3; %2; %2;%3;%4;
para_.maxConfToSelect_Hotspot_2t = 5; %3;%4;


%-- minimum desired map quality improvement for the additional conf
para_.minERQImprovement_1t = 0.10;
para_.minERQImprovement_2t = 0.10;


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%       GDM PLANNING - FUSION PARAMETERS
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

para_.maxFusionIterations_1t = 07;  %10; %5
para_.maxFusionIterations_2t = 07; %20; %12; % sample-env-paper


% para_.toleranceERQ_Fusion_percentage = 05;
para_.toleranceERQFusion1t = 0.25; %0.15;
para_.toleranceERQFusion2t = 0.25; %0.15; % to be fixed


para_.desiredERQFusion1t = 0.40; %0.50;
para_.desiredERQFusion2t = 0.45; %0.55;





%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%       PUBLISH PARAMETERS
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
para_.PublishTextWidth = 'single-column';
% para_.PublishTextWidth = 'two-columns';

% para_.PublishUpperPPM = 10000;
%para_.PublishUpperPPM = 20000; %30000;  % sample-env-paper (coarse map)
 para_.PublishUpperPPM = 6000; %30000;  % sample-env-paper (coarse map)

%para_.PublishUpperPPM = inf; % sample-env-paper (ground truth, final recon)


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%       EVALUATION PARAMETERS
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
para_.concentrationThresholdEvaluation = 500; %250;%100; %5;

para_.SCsCutoffWeightFactor = 1; %0.8;  %0.5;
para_.SCsCutoffPPMM = 2500; %1000;
para_.SCsClusterDist_cells = 10;
% para_.SCsCombiningDist_cells = 05; %10
para_.SCsCombiningDist_cells = 14; % for sample-env-paper %05 for paper evaluation 
para_.ScCutOffClusters = 3;

para_.ScPositiveHits_cells = 06; %06; %05;


para_

end

function fTwoStepExplorationICRA16( para_ )
%
%



alpha = 0.75;
beta  = 0.50;
gamma = 0.00;


%//////////////////////////////////////
%   DIRECTORIES
%//////////////////////////////////////
para_.MissionStrategy = '2t-armex-icra16';
[ dir_ ] = fDirTwoStepXXX( para_,alpha,beta,gamma );



% -----------------------------------------------------------------------------------
% 
% -----------------------------------------------------------------------------------
FoV           = para_.FieldOfViewDEG;
FoVragd       = para_.FieldOfViewDEG;%/2;
SensorRange.m = para_.SensingRangeM;

% -----------------------------------------------------------------------------------
%                                   YAML file data
% -----------------------------------------------------------------------------------
% ParamPreprocess.MapYAML
% fileYAML = ReadYamlRaw(ParamPreprocess.MapYAML);
% fileYAML = ReadYamlRaw('map.yaml');

dataYAML = ReadYamlRaw([dir_.env,para_.EnvironmentYamlFile]);
dataYAML_origin = cell2mat(dataYAML.origin);

robot_origin_x = abs(dataYAML_origin(1)/dataYAML.resolution);
robot_origin_y = abs(dataYAML_origin(2)/dataYAML.resolution);
% map_pix_size   = dataYAML.resolution;

% robot_origin_in_map_x = 1;
% robot_origin_in_map_y = 1;
% dataYAML.resolution   = 1;
% -----------------------------------------------------------------------------------

fprintf(1,['\n\n',...
           '*********************************************************************\n\n',...
           '                  PLANNED+ADAPTIVE EXPLORATION \n',...
           '             Two-step exploration strategy (ICRA16) \n\n',...
           '*********************************************************************\n\n'])


% ------------------ Gas Detection -----------------------
solve_detection.preprocess                      = 0;
solve_detection.spp                             = 0;
solve_detection.tsp                             = 0;
solve_detection.robot_plan                      = 0;
solve_detection.ground_truth                    = 0;
solve_detection.measurements                    = 0;
solve_detection.recons                          = 0; % build gdm from gas detection.

publish_detection.spp                           = 0;
publish_detection.recon                         = 0;
publish_detection.coverage                      = 0;
publish_detection.gt                            = 0;
publish_detection.env                           = 0;
publish_detection.detection_events              = 0;
 
% ------------------ Gas Tomography ----------------------
solve_tomography.preprocess_maps                = 1; % preprocess: maps
solve_tomography.spp_hotspots                   = 0; % preprocess: coverage and cross angles of hotspots
solve_tomography.fusion                         = 0;
solve_tomography.tsp_spp                        = 0; % run tsp on the final solution

solve_tomography.measurements_spp_adaptive      = 0;
solve_tomography.measurements_spp_fixed         = 0;

solve_tomography.recons_spp_adaptive            = 0;
solve_tomography.recons_spp_fixed               = 0;

solve_tomography.evaluation_spp_adaptive        = 0;
solve_tomography.evaluation_spp_fixed           = 0;
solve_tomography.evaluation_spp_adaptive2       = 0;
solve_tomography.evaluation_spp_fixed2          = 0;

publish_tomography.maps                         = 0; % publish generated maps
publish_tomography.spp_hotspots                 = 0; % publish spp solutions for hotspots
publish_tomography.spp_hotspots_all             = 0; % publish spp overall solution for hotspots

publish_tomography.spp_fusion                   = 0; % publsih spp solution for fusion
publish_tomography.spp_solution                 = 0; % publsih final spp solution

publish_tomography.spp_plan_adaptive            = 0;

publish_tomography.coverage_spp_adaptive        = 0;
publish_tomography.coverage_spp_fixed           = 0;

publish_tomography.recon_spp_adaptive           = 0;
publish_tomography.recon_spp_fixed              = 0;

publish_tomography.exploration                  = 0; % publsih final spp solution






%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%                   GAS DETECTION
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

fSolveDetectionPlan( solve_detection,para_,dir_ )
fPublishDetectionPlan( publish_detection,para_,dir_ )

para_
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%                   GAS TOMOGRAPHY
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
fSolveTomographyPlanICRA16( alpha,beta,gamma,solve_tomography,para_,dir_ )
fPublish_GDMPlanPlanICRA16( publish_tomography,para_,dir_ )
            





end


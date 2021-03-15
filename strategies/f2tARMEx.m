function f2tARMEx( para_ )
%
%=====================================================================
%                       MAIN OF 2t-ARMEx
%=====================================================================
%

fprintf(1,['\n\n',...
'*********************************************************************\n\n',...
'                            2t-ARMEx \n',... 
'                   TWO-TOUR MISSION STRATEGY \n\n',...
'*********************************************************************\n\n'])


%===================================================
%      Gas Detection
%==================================================

% Select what you want me to do from the following options.

%--- For planning
%{ 
    solve_gd.*
    
        preprocess: To perform preprocessing for the selection of sensing
        configurations, i.e., playing with the geometric map resolutions,
        defiining candidate sensing configurations, computing visibilty
        matric etc. 

        spp: To solve sensor placement problem for gas detection, i.e.,
        selecting a minimal set of sensing confiugrations that provides
        full sensing coverage.

        tsp: To sort out the selected configurations using TSP.

        robot_plan: to translate the selected positions into the robot's
        understandble vectors (for ROS implementations). 

        ground_truth: Ground truth of gas distributions.

        measurements:Collecting measurements according to the plan.

        recons: To create (coarse) gas distribution from the collected
        measurements. 

        coverage: To compute coverage (after plan execution).

        traveling_path_dist: To compute the path distance (after plan
        execution).

%}

%--- For visualization
%{
    publish_gd.*
    
        env: To publish geometric map of the environment.

        gt: To publish ground truth gas distributions.

        spp: To publish the measurement plan.

        recon: To Publish the reconstructed gas distrubution map.

        coverage: To publish the provided sensing coverage.

        detection_events: To mark sensing detection events.

        rwl1_iterations: To publish detailed results of re-weighted l1
        minimization. 
%}

%----------- The Plan ----------------------------------
solve_gd.preprocess                      = 0; 
solve_gd.spp                             = 0;
solve_gd.tsp                             = 0;
solve_gd.robot_plan                      = 0;
solve_gd.ground_truth                    = 0;
solve_gd.measurements                    = 0;
solve_gd.recons                          = 0;
solve_gd.coverage                        = 0;
solve_gd.traveling_path_dist             = 0; 
%----------- Publish Results -----------------------
publish_gd.env                           = 0;
publish_gd.gt                            = 0;
publish_gd.spp                           = 0;
publish_gd.recon                         = 0;
publish_gd.coverage                      = 0;
publish_gd.detection_events              = 0;
publish_gd.rwl1_iterations               = 0;

%===================================================
%      Gas Distribution Mapping
%===================================================

% Select what you want me to do from the following options.

%--- For planning
%{ 
    solve_gdm.*
    
        preprocess_maps: To perform preprocessing and hotspots estimation
        on the corase GDM from gas detection.

        optimal_hotspot_geometries: To find optimal sensing geometries for
        each hotspot. 

        tsp: To sort out the selected configurations using TSP.

        measurements: Collecting measurements according to the plan. 

        recons: To create (accurate) gas distribution from the collected
        measurements. 

        coverage: To compute sensing coverage (after plan execution).

        traveling_path_dist: To compute the path distance (after plan
        execution).

        evaluation/evaluation2/evaluation3:To evaluate the reconstructed
        gas distribution against true gas distributions.

%}

%--- For visualization
%{
    publish_gdm.*
    
        ground_truth: To publish ground truth gas distributions.

        maps: To publish hotspot estimations on the corase GDM.

        spp_hotspots: To individually publish optimal sensing geometries
        for each hotspot.

        spp_hotspots_all: To publish all the selected sensing geometries.

        spp_fusion: To publish the fused geometries.

        spp_solution: To publish the final plan.

        recon: To Publish the reconstructed gas distrubution map.

        spp_plan_adaptive: adaptive sensor planning for gas distribution
        mapping.

        coverage: To publish the provided sensing coverage.

        exploration: To publsih final spp solution.
%}


%----------- Plan ----------------------------------
solve_gdm.preprocess_maps                = 0; 
solve_gdm.optimal_hotspot_geometries     = 0;  
solve_gdm.fusion_optimal_geometries      = 1;
solve_gdm.tsp_spp                        = 0; 
solve_gdm.measurements                   = 0;
solve_gdm.recons                         = 0;
solve_gdm.coverage                       = 0;
solve_gdm.traveling_path_dist            = 0;
%----------- Evaluation ----------------------------
solve_gdm.evaluation                     = 0;
solve_gdm.evaluation2                    = 0;
solve_gdm.evaluation3                    = 0;
%----------- Publish Results -----------------------
publish_gdm.ground_truth                 = 0;
publish_gdm.maps                         = 0; 
publish_gdm.spp_hotspots                 = 0; 
publish_gdm.spp_hotspots_all             = 0; 
publish_gdm.spp_fusion                   = 0; 
publish_gdm.spp_solution                 = 0; 
publish_gdm.recon                        = 0;
publish_gdm.spp_plan_adaptive            = 0;
publish_gdm.coverage                     = 0;
publish_gdm.exploration                  = 0; 

%===================================================
%      Common
%===================================================
% Select if you want me to create video for 2t-ARMEx.
publish_video                            = 0;
%--------------------------------------------------


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%       DIRECTORIES
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
para_.MissionStrategy = '2t-armex';
para_.FilePrefix = sprintf('%s-2t',para_.FilePrefix2);
[ dir_ ] = fDir2tARMEx( para_ );


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%       DIARY
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DateTimeNow = datestr(now,'yyyy-mm-dd-HH:MM:SS');
diaryfilename = sprintf('diary_%s',DateTimeNow);
diary([dir_.diary,diaryfilename]);


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%       GAS DETECTION
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
fPlan_2tARMEx_GD   (   solve_gd,para_,dir_ )
fPublish_2tARMEx_GD( publish_gd,para_,dir_ )

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%       GAS DISTRIBUTION MAPPING
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
fPlan_2tARMEx_GDM   (   solve_gdm,para_,dir_ )
fPublish_2tARMEx_GDM( publish_gdm,para_,dir_ )


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%       PUBLISH - VIDEO
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if publish_video == 1        
    fprintf(1,['\n\n',...
    '*********************************************\n',...
    '      PUBLISH 2t-ARMEx EXPLORATION VIDEO \n',...
    '*********************************************\n\n'])
    % publish
    % ---------------------
    fVideo_ARMEx( para_,dir_ )
end


diary save
diary off

end


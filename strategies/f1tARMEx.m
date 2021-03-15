function f1tARMEx( para_ )
%
% ONE-TOUR MISSION STRATEGY (1t-ARMEx)
%
% @ Asif Arain
%
%---------------------------------------------------------------------


fprintf(1,['\n\n',...
'*********************************************************************\n\n',...
'                            1t-ARMEx \n',... 
'                   ONE-TOUR MISSION STRATEGY \n\n',...
'*********************************************************************\n\n'])

%-------------------------------------- 
% INTIAL EXPLORATION PLAN
%-------------------------------------- 
solve_detection_initial_preprocess  = 0;
solve_detection_initial_spp         = 0;
solve_detection_initial_tsp         = 0;
solve_detection_initial_ees         = 0; %exploitation_enabled_sort
%-------------------------------------- 
% ADAPTIVE EXPLORATION AND EXPLOITATION
%-------------------------------------- 
solve_adaptive_plans                = 1;
solve_reconstruction                = 1;
solve_evaluation                    = 0;
solve_coverage                      = 0;
solve_traveling_path_dist           = 0;
%--------------------------------------
% PUBLISH RESULTS
%-------------------------------------- 
publish_executed_plan               = 1;
publish_reconstruction              = 1;
publish_sensing_coverage            = 0;
publish_video                       = 0;

publish_step_wise_plan_for_presentation = 0;
% ------------------------------------- 



%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%       DIRECTORIES
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
para_.MissionStrategy = '1t-armex';
para_.FilePrefix = sprintf('%s-1t',para_.FilePrefix2);
[ dir_ ] = fDir1tARMEx( para_ );


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
%       MAP AND RELATED STUFF
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[ map_env,...
  map_conf,...
  origin_env,...
  origin_conf,...
  cellsize_env,...
  cellsize_conf,...
  cellsize_org ] = fMapStuff( para_,dir_ );



% ----------------------
%   SENSING PARAMETERS
% ----------------------
FoV           = para_.FieldOfViewDEG;
%FoVragd       = para_.FieldOfViewDEG;%/2;
SensorRange.m = para_.SensingRangeM;
SensorRange.cell = SensorRange.m/cellsize_env;



%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%       INITIAL GAS DETECTION PLAN
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++       
       
%*************************************************
%   PREPROCESS
%*************************************************
if solve_detection_initial_preprocess == 1
    
    fprintf(1,['\n\n',...
    '*********************************************\n',...
    '   PREPROCESS - INITIAL GAS DETECTION PLAN \n',...
    '*********************************************\n\n'])
    visualize = 0;
    
    % preprocessing parametes:
    % --------------------------
    para_.GraphCnt                 = 8; % graph connectivity (4/8).
    para_.numFoVs                  = 8; % nummber of elementray sensing action per cell.
    para_.Ct                       = 0; % travelling cost matrix (yes/no).    
    para_.TravCost                 = 10; % travelling cost (sec) for one meter travel.    
    para_.A                        = 0; % generate connectivity matrix A (yes/no).
    para_.AWM                      = 0; % artificial world map (yes/no).
    para_.do_nTCell                = 0;
    para_.nTLength                 = 3; % non traversable area around an obstacle (meters).
    %para_.SensorDomainComputation = 'ComputeNew';
    para_.SensorDomainComputation  = 'UseEarlier';
    
    % preprocess
    [ SensorRange,...
      OBSenv,...
      OBSconf,...
      V,...
      oV,...
      FoVfaces,...
      confKept,...
      confRemoved,...
      cellKept,...
      cellRemoved,...
      T,...
      E,...
      eE,...
      A,...
      oA,...
      Ct,...
      o,...
      P,...
      tPreprocess ] = fPreprocessDetection( FoV,...
                                            SensorRange,...
                                            para_,...
                                            dir_,...
                                            visualize );
    
    save([dir_.sppdetection,'preprocess_detection.mat'],...
        'SensorRange',...
        'OBSenv',...
        'OBSconf',...
        'V',...
        'oV',...
        'FoVfaces',...
        'confKept',...
        'confRemoved',...
        'cellKept',...
        'cellRemoved',...
        'T',...
        'E',...
        'eE',...
        'A',...
        'oA',...
        'Ct',...
        'o',...
        'P',...
        'tPreprocess');
      
end


%*************************************************
%	SPP
%*************************************************
if solve_detection_initial_spp == 1
    
    fprintf(1,['\n\n',...
    '*********************************************\n',...
    '        SPP - INITIAL GAS DETECTION PLAN \n',...
    '*********************************************\n\n'])
    
    visualize = 0;
    load([dir_.sppdetection,'preprocess_detection.mat'],...
        'V',...
        'oV',...
        'FoVfaces',...
        'confKept');
    
    %----------------------------------------------------------------
    fprintf('\n******** SOLUTION FOR THE FULL COVERAGE ********\n');
    %----------------------------------------------------------------
    
    % spp parameters
    % --------------------------------------------------------------
    % para_.Algorithm = 'Comb';  % Options: 'Comb','convSPP'
    % para_.AlgoRev   = 2;       % Revision.
    para_.Algorithm = 'convSPP';  % Options: 'Comb','convSPP'
    para_.AlgoRev   = 5;       % Revision.
    para_.Prog      = 'GRB';   % Options: 'CVX','GRB','MOT'
    para_.numHistIR = 5;       % 
    para_.maxIterat = 150;     % Maximum number of iterations.
    para_.numCu     = 80;      % Number of uncertain confs. to terminate RWL1.
    
    % spp
    [ C,...
      C0,...
      C1,...
      Cu,...
      combC,...
      lBound,...
      uBound,...
      logs_rwl1,...
      tSPP_rwl1,...
      tSPP_Comb,...
      tSPP ] = fSPPDetection( V,para_ );
    
    %------------------------------------------------------------------
    fprintf('\n******** SOLUTION FOR THE REDUCED COVERAGE ********\n');
    %------------------------------------------------------------------
    [ C,...
      reducedCoverage ] = fReducedCoverageConfs( C,...
                                                 V,...
                                                 oV,...
                                                 FoV,...
                                                 FoVfaces,...
                                                 confKept,...
                                                 SensorRange,...
                                                 map_env,...
                                                 para_,...
                                                 dir_,...
                                                 visualize );
    
    save([dir_.sppdetection,'detection_spp_initial.mat'],...
        'C',...
        'C0',...
        'C1',...
        'Cu',...
        'combC',...
        'lBound',...
        'uBound',...
        'reducedCoverage',...
        'logs_rwl1',...
        'tSPP_rwl1',...
        'tSPP_Comb',...
        'tSPP');
end


%*************************************************
%	TSP
%*************************************************
if solve_detection_initial_tsp == 1

    %profile on;
    fprintf(1,['\n\n',...
    '*********************************************\n',...
    '        TSP - INITIAL EXPLORATION PLAN \n',...
    '*********************************************\n\n'])
    load([dir_.sppdetection,'preprocess_detection.mat'],...        
            'oV',...
            'o',...
            'E',...
            'T',...
            'OBSenv',...
            'FoVfaces',...
            'confKept');
    load([dir_.sppdetection,'detection_spp_initial.mat'],'C');

    Conf = zeros(size(oV,2),1);
    Conf(confKept) = C;

    % -- conf vector for tsp.
    ConfTSP = Conf';
    f = FoVfaces.num;
    robotStartPosition_ind = sub2ind(size(map_env),origin_conf(1),...
                                                   origin_conf(2));
    robotStartConf_ind = ((robotStartPosition_ind-1)*f)+1;
    ConfTSP(robotStartConf_ind) = 1;
    % robotStartPosition_ind = sub2ind(size(map_env),round(confs_executed(end,1)-0.5),...
    %                                                round(confs_executed(end,2)-0.5));
    % robotStartConf_ind = ((robotStartPosition_ind-1)*f)+1;
    % ConfTSP(robotStartConf_ind) = 1;

    % -- removing conf on occupied cells
    % --------------------------------
    % obs-conf nums
    confOBS = zeros(numel(OBSenv.ind)*f,1);
    for i = 1:numel(OBSenv.ind)
        confOBS(((i-1)*f)+1:((i-1)*f)+f) = ...
            (((OBSenv.ind(i)-1)*f)+1):(((OBSenv.ind(i)-1)*f)+f);
    end
    ConfTSP(:,confOBS) = [];
    ConfTSP = ConfTSP';

    % -- TSP PARAMETERS
    % --------------
    para_.TravCost   = 10; %ParamPreprocess.TravCost;
    para_.ConfOBS    = 1;
    para_.EoTinAstar = 1; %0; 1 for simplified path A*, 0 for complex path A*

    [ tRoute,...
      tCost,...
      tTSP ] = fTSP2( o,...
                      FoVfaces,...
                      ConfTSP,...
                      E,...
                      T,...
                      map_env,...
                      OBSenv,...
                      para_ );

    % -- traveling route with sequence starting from initial position            
    % -------------------------------------
    tRouteOrderedWRTStart = tRoute(1:end-1,:);
    tRouteOrderedWRTStart = circshift(tRouteOrderedWRTStart,numel(tRouteOrderedWRTStart)-...
        find(tRouteOrderedWRTStart==robotStartConf_ind)+1);

    % -- conf sequence number in the order of conf vector index
    % (start position is not a conf) 
    % -------------------------------------
    [~,confSequenceNum,~] = intersect(tRouteOrderedWRTStart,find(Conf));

    % -- since start position is not a conf,
    % -------------------------------------
    confSequenceNum = confSequenceNum-1;

    % -- save results
    % -------------------------------------
    file_name = sprintf('detection_tsp_initial.mat');
    save([dir_.sppdetection,file_name],...
        'tRoute',...
        'tCost',...
        'tTSP',...
        'confSequenceNum',...
        'tRouteOrderedWRTStart');
    %---------------
    %profile_name = 'profile_tsp22';
    %mkdir(profile_name)
    %profile off;
    %profsave(profile('info'),profile_name);
end


%*************************************************
%   EXPLOITATION ENABLED SORTED CONFIGURATIONS
%*************************************************
if solve_detection_initial_ees == 1

    %profile on;
    
    fprintf('\n\n ****************************************************\n');
    fprintf('   EXPLORATION ENABLED SORT - INITIAL GAS DETECTION \n');
    fprintf(' ****************************************************\n');
    
    load([dir_.sppdetection,'preprocess_detection.mat'],...        
            'oV',...
            'o',...
            'E',...
            'T',...
            'OBSenv',...
            'FoVfaces',...
            'confKept');

    load([dir_.sppdetection,'detection_spp_initial.mat'],'C');
    
    confs_vec = zeros(size(oV,2),1);
    confs_vec(confKept) = C;
    Confs_num = find(confs_vec);
    
    f = FoVfaces.num;
    confCell_ind = fix(Confs_num/f)+(~~mod(Confs_num,f));
    [confCell_r,confCell_c] = ind2sub(size(map_env),confCell_ind);
    
    sensing_confs = [confCell_r,confCell_c];
    initial_position = origin_conf;
    
    [ confSequenceNum,...
      sensing_confs_sorted ] = fExploitationSortedConfs( sensing_confs,...
                                                         initial_position,...
                                                         map_env );
                         
    % -- save results
    % -------------------------------------
    file_name = sprintf('detection_ees_initial.mat');
    save([dir_.sppdetection,file_name],...        
        'confSequenceNum',...
        'sensing_confs_sorted');

end



%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%       ADAPTIVE MEASUREMENT PLANNING
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if solve_adaptive_plans == 1
    
    fprintf(1,['\n\n',...
    '****************************************************\n',...
    '              PLANNING 1t-ARMEx \n',...
    '****************************************************\n\n']);
    
    visualize = 0;   

    load([dir_.sppdetection,'preprocess_detection.mat'],'SensorRange');

    % parameters:
    % --------------------------
    % nummber of elementray sensing action per cell.
    para_.numFoVs = 360; 
    para_.SensorDomainComputation = 'UseEarlier'; % 'ComputeNew'
    % circular area for the candidate conf around redundant confs
    para_.candConfRadiusForRedConfs_m = 1.5*SensorRange.cell; 
    % min radius where sensing conf are not allowed.
    para_.innerRadius_m = 1; 
    % number of allowed cand conf.
    para_.NumOfCandConf = 100;%200; 
    % min distance to check redundant confs.
    para_.DistanceForRedundantConfs_m = 3*para_.SensingRangeM; %2;
    % replanning threshould (cells)
    % para_.ReplanningThreshould = 2;
    para_.ReplanningThreshould = SensorRange.cell/2;
    para_.ResEstimatedHcOffset = SensorRange.cell/1.5; %SensorRange.cell/2;
    para_.ReEstimatedHcpath2euclRatio = 1.5;    
    % -- cocentration threshould to perform gas tomography
    para_.ConcentrationThreshould4TomographyPPM = 500;
    
    %para_.TravelingType1SE = 'TSP';
    %para_.TravelingType1SE = 'TSP2';
    para_.TravelingType1SE = 'ExploitationEnabledSort';
    %6 .96
    retrieve1SE_pts = [];
    %retrieve1SE_pts = load([dir_.str,'BackupPtsSequence1SE.dat']);
    % profile on;    

    fPlan_1tARMEx( FoV,...
                   cellsize_env,...
                   retrieve1SE_pts,...
                   dir_,...
                   para_,...
                   visualize );
    % profile_name = sprintf('1SE_tmp');
    % mkdir(profile_name)            
    % profile off;
    % profsave(profile('info'),profile_name);
    file_name = sprintf('%s-online-plot',para_.ExperimentTitle);
    export_fig([dir_.str,file_name], '-pdf','-r500');
    
    fprintf(1,['\n\n',...
    '*************** end of adaptive strategy ****************\n']);
end



%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%       FINAL ACCURATE GDM
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if solve_reconstruction == 1 
    
    fprintf(1,['\n\n',...
    '****************************************************\n',...
    '                  ACCURATE GDM \n',...
    '****************************************************\n\n']);
       
    %dir_.logs = dir_.MeasurementLogs; 
    dir_.Recon = dir_.recon;
    visualize = 0;
    
    [ M,...
      mean_map,...
      cell_coords_x,...
      cell_coords_y ] = fReconstructionFixed( para_,...
                                              dir_,...
                                              visualize );
                                               
    save([dir_.recon,'accurate_gdm.mat'],...
        'M','mean_map','cell_coords_x','cell_coords_y');
    
    
end



%***************************************************************
% 
%                 COVERAGE CALCULATION
% 
%***************************************************************
if solve_coverage == 1 
    
    fprintf(1,['\n\n',...
    '*********************************************\n',...
    '          SPP SENSING COVERAGE \n',...
    '*********************************************\n\n'])

    para_.ConfType = 'Planned'; %'Executed'
    %para_.ConfType = 'Executed';
    
    sensingCoverage = fCalculateCoverage( map_env,SensorRange,FoV,para_,dir_ );
    save([dir_.str,'sensing_coverage.mat'],'sensingCoverage');
    
end



%***************************************************************
% 
%                 FINAL TRAVELING PATH/DISTANCE
% 
%***************************************************************
if solve_traveling_path_dist == 1 
    
    fprintf(1,['\n\n',...
    '*********************************************\n',...
    '         SPP TRAVELING PATH/DISTANCE \n',...
    '*********************************************\n\n'])
       
    para_.ConfType = 'Planned'; %'Executed'
    
    [ travelingPath,travelingDist ] = fTravelingPathDist( map_env,para_,dir_ );
    save([dir_.str,'traveling_path_dist.mat'],...
        'travelingPath','travelingDist');
    
end



%***************************************************************
% 
%                        EVALUATION
% 
%***************************************************************
if solve_evaluation == 1
    
    fprintf(1,['\n\n',...
           '*********************************************\n',...
           '                EVALUATION \n',...
           '*********************************************\n\n'])
       
    visualize = 1;
    %recon_filename = 'reconstruction_OneStepExploration.mat';
    recon_filename = 'accurate_gdm.mat';
    
    %dir_.recon = dir_.Reconstruction;
    
    
    
    
    dir_
    
    dir_.ROSLogsThis = dir_.ROSLogs
    
    dir_.recon = dir_.recon %dir_.gdmplan_recon;
    
    
    [ sc_true,...
       sc_estimated,...
       true_positives,...
       false_positives,...
       num_of_true_positives,...
       num_of_false_positives,...
       true_positives_dist,...
       false_positives_dist,...
       dist_estimatedSources2trueSources,...
       precision,...
       recall,...
       f_measure ] = fEvaluation4from3old( recon_filename,...
                                   map_env,...
                                   cellsize_env,...
                                   para_,...
                                   dir_,...
                                   visualize )
%                                pause
    
%     [ sc_true,...
%            sc_estimated,...
%            true_positives,...
%            false_positives,...
%            num_of_true_positives,...
%            num_of_false_positives,...
%            true_positives_dist,...
%            false_positives_dist,...
%            dist_estimatedSources2trueSources,...
%            precision,...
%            recall,...
%            f_measure ] = fEvaluation3old( recon_filename,...
%                                        map_env,...
%                                        cellsize_env,...
%                                        para_,...
%                                        dir_,...
%                                        visualize )
%     pause
%     
%     
%     
%     
%     
%     
%     
%     [ true_positives,...
%       false_positives,...
%       num_of_true_positives,...
%       num_of_false_positives,...
%       precision,...
%       recall,...
%       f_measure ] = fEvaluation3( recon_filename,...
%                                   map_env,...
%                                   cellsize_env,...
%                                    para_,...
%                                    dir_ );
                               
    save([dir_.evaluation,'evaluation.mat'],...
        'true_positives',...
        'false_positives',...      
        'num_of_true_positives',...
        'num_of_false_positives',...
        'precision',...
        'recall',...
        'f_measure');
    
    % [ sc_true,...
    %   sc_estimated,...
    %   true_positives,...
    %   false_positives,...
    %   num_of_true_positives,...
    %   num_of_false_positives,...
    %   true_positives_dist,...
    %   false_positives_dist,...
    %   dist_estimatedSources2trueSources,...
    %   precision,...
    %   recall,...
    %   f_measure ] = fEvaluation3( recon_filename,...
    %                               map_env,...
    %                               cellsize_env,...
    %                               para_,...
    %                               dir_,...
    %                               visualize );
    % 
    % save([dir_.Evaluation,'evaluation.mat'],...
    %     'sc_true',...
    %     'sc_estimated',...
    %     'true_positives',...
    %     'false_positives',...      
    %     'num_of_true_positives',...
    %     'num_of_false_positives',...
    %     'true_positives_dist',...
    %     'false_positives_dist',...
    %     'dist_estimatedSources2trueSources',...
    %     'precision',...
    %     'recall',...
    %     'f_measure');
end




%***************************************************************
% 
%                PUBLISH - EXECUTED PLAN (STEP-WISE FOR PRESENTATION)
%
%***************************************************************
if publish_step_wise_plan_for_presentation == 1
        
    fprintf(1,['\n\n',...
           '*********************************************\n',...
           '        PUBLISH 1T-RMS EXECUTED PLAN \n',...
           '*********************************************\n\n'])
    
    load([dir_.sppdetection,'preprocess_detection.mat'],'SensorRange','FoVfaces');
    
    %***********************************
    % ONE
    %***********************************
    hcs_sub = load([dir_.str,'hc_list_final.txt']);
    hcs_ind = sub2ind(size(map_env),hcs_sub(:,1),hcs_sub(:,2));    
    executed_confs = load([dir_.str,'executed_confs_final.txt']);    
    executed_confs_ind = sub2ind(size(map_env),...
                                 round(executed_confs(:,1)-0.5),...
                                 round(executed_confs(:,2)-0.5) );
    executed_confs_orn = executed_confs(:,3);
    
    
    conf_ind = [executed_confs_ind(1,:)]
    conf_ord = [executed_confs_orn(1,:)]
    
    publish_gdm_for_executed_conf_num = 1
    
    % publish parameters
    %--------------------------------------------------------------
    para_.PublishFoV                    = 0;  % draw field of view.
    para_.PublishSymbolicFoV            = 1;
    para_.PublishCandConf               = 0;
    para_.PublishGDM                    = 1;
    para_.PublishHotspots               = 1;
    para_.PublishSubH                   = 0;
    para_.PublishStartPosition          = 0;
    para_.PublishConfPosition           = 1;
    para_.PublishArrow                  = 1;
    para_.PublishConfNum                = 1;
    para_.PublishFusedConfArrow         = 0;
    para_.PublishRedundantConfArrow     = 0;
    para_.PublishWeakConf               = 0;
    para_.PublishCondConf               = 0;
    para_.PublishHotspotsInOutTSPAngles = 0;
    para_.PublishReferenceDistPoint     = 0;
    para_.PublishGroundTruth            = 0;
    para_.PublishSensingCoverage        = 0;
    para_.PublishRedundantConfArrow     = 0;
    para_.PublishFusedConfArrow         = 0;
    para_.PublishAdaptivePlannedConf    = 0;
    para_.PublishExploration            = 0;
    para_.ReconstructionFile            = 'accurate_gdm.mat';
    %para_.PublishTextWidth              = 'two-columns';
    para_.PublishFigureTitle            = '1t-ARMEX step 1';

    % data to publish
    % -------------------------------------------------------------
    dataPublish.FoV                     = FoV;
    dataPublish.FoVfaces                = FoVfaces;    
    dataPublish.SensorRange             = SensorRange;
    dataPublish.robotStartCellInMap     = origin_conf;
    dataPublish.ConfNum                 = conf_ind;
    dataPublish.ConfOrn                 = conf_ord;
    dataPublish.confSequenceNum         = 1:numel(executed_confs_ind);
    dataPublish.hotspots                = hcs_ind(1);
    dataPublish.map_env                 = map_env;
    dataPublish.publish_gdm_for_executed_conf_num = publish_gdm_for_executed_conf_num
    
    % publish
    %--------------------------------------------------------------
    fPublish_GDMPlan( dataPublish,para_,dir_ )
        
    % save plots
    %--------------------------------------------------------------
    file_name = sprintf('%s-1t-ARMEx-step1',para_.ExperimentTitle);
    %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    %print('-painters','-dpdf','-r500',[dir_.TomographyResults,file_name]); % pdf
    %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png
    export_fig([dir_.str,file_name], '-pdf','-r500')
    
    %pause
    
    %************************************************************
    % TWO
    %************************************************************
    hcs_sub = load([dir_.str,'hc_list_final.txt']);
    hcs_ind = sub2ind(size(map_env),hcs_sub(:,1),hcs_sub(:,2));    
    executed_confs = load([dir_.str,'executed_confs_final.txt']);    
    executed_confs_ind = sub2ind(size(map_env),...
                                 round(executed_confs(:,1)-0.5),...
                                 round(executed_confs(:,2)-0.5) )
    executed_confs_orn = executed_confs(:,3)
    
    %--
       
    reselectedConf = load([dir_.sppgdm,'optimal_geometry_hc01_spp.dat']);
    
    load([dir_.sppgdm,'optimal_geometry_hc01_xvd_LastSolution.mat']);
    load([dir_.sppgdm,'optimal_geometry_hc01_vd.mat']);
    load([dir_.sppgdm,'optimal_geometry_hc01_xvd.mat']);
       
    reselectedConf_num = find(reselectedConf)
    reselectedConf_ind = confKept(reselectedConf_num)
    reselectedConf_orn = ConfOrns(reselectedConf_num) %oConfOrientations(confKept(reselectedConf_num));
    
    
    conf_ind = [executed_confs_ind(1:2,:);reselectedConf_ind']
    conf_ord = [executed_confs_orn(1:2,:);reselectedConf_orn']
    
    publish_gdm_for_executed_conf_num = 2
    
    % publish parameters
    %--------------------------------------------------------------
    para_.PublishFoV                    = 0;  % draw field of view.
    para_.PublishSymbolicFoV            = 1;
    para_.PublishCandConf               = 0;
    para_.PublishGDM                    = 1;
    para_.PublishHotspots               = 1;
    para_.PublishSubH                   = 0;
    para_.PublishStartPosition          = 0;
    para_.PublishConfPosition           = 1;
    para_.PublishArrow                  = 1;
    para_.PublishConfNum                = 1;
    para_.PublishFusedConfArrow         = 0;
    para_.PublishRedundantConfArrow     = 0;
    para_.PublishWeakConf               = 0;
    para_.PublishCondConf               = 0;
    para_.PublishHotspotsInOutTSPAngles = 0;
    para_.PublishReferenceDistPoint     = 0;
    para_.PublishGroundTruth            = 0;
    para_.PublishSensingCoverage        = 0;
    para_.PublishRedundantConfArrow     = 0;
    para_.PublishFusedConfArrow         = 0;
    para_.PublishAdaptivePlannedConf    = 0;
    para_.PublishExploration            = 0;
    para_.ReconstructionFile            = 'accurate_gdm.mat';
    %para_.PublishTextWidth              = 'two-columns';
    para_.PublishFigureTitle            = 'One step exploration plan for gas tomography';

    % data to publish
    % -------------------------------------------------------------
    dataPublish.FoV                     = FoV;
    dataPublish.FoVfaces                = FoVfaces;    
    dataPublish.SensorRange             = SensorRange;
    dataPublish.robotStartCellInMap     = origin_conf;
    dataPublish.ConfNum                 = conf_ind;
    dataPublish.ConfOrn                 = conf_ord;
    dataPublish.confSequenceNum         = 1:numel(executed_confs_ind);
    dataPublish.hotspots                = hcs_ind;
    dataPublish.map_env                 = map_env;
    dataPublish.publish_gdm_for_executed_conf_num = publish_gdm_for_executed_conf_num
    
    % publish
    %--------------------------------------------------------------
    fPublish_GDMPlan( dataPublish,para_,dir_ )
        
    % save plots
    %--------------------------------------------------------------
    file_name = sprintf('%s-1t-ARMEx-step2',para_.ExperimentTitle);
    %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    %print('-painters','-dpdf','-r500',[dir_.TomographyResults,file_name]); % pdf
    %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png
    export_fig([dir_.str,file_name], '-pdf','-r500')
    
    
    %************************************************************
    % THREE
    %************************************************************
    hcs_sub = load([dir_.str,'hc_list_final.txt']);
    hcs_ind = sub2ind(size(map_env),hcs_sub(:,1),hcs_sub(:,2));    
    executed_confs = load([dir_.str,'executed_confs_final.txt']);    
    executed_confs_ind = sub2ind(size(map_env),...
                                 round(executed_confs(:,1)-0.5),...
                                 round(executed_confs(:,2)-0.5) )
    executed_confs_orn = executed_confs(:,3)
    
    %--
       
    reselectedConf = load([dir_.sppgdm,'optimal_geometry_hc02_spp.dat']);
    
    load([dir_.sppgdm,'optimal_geometry_hc02_xvd_LastSolution.mat']);
    load([dir_.sppgdm,'optimal_geometry_hc02_vd.mat']);
    load([dir_.sppgdm,'optimal_geometry_hc02_xvd.mat']);
    
    reselectedConf_num = find(reselectedConf)
    reselectedConf_ind = confKept(reselectedConf_num)
    reselectedConf_orn = ConfOrns(reselectedConf_num) %oConfOrientations(confKept(reselectedConf_num));
        
    conf_ind = [executed_confs_ind(1:3,:);reselectedConf_ind']
    conf_ord = [executed_confs_orn(1:3,:);reselectedConf_orn']
    
    publish_gdm_for_executed_conf_num = 3
    
    % publish parameters
    %--------------------------------------------------------------
    para_.PublishFoV                    = 0;  % draw field of view.
    para_.PublishSymbolicFoV            = 1;
    para_.PublishCandConf               = 0;
    para_.PublishGDM                    = 1;
    para_.PublishHotspots               = 1;
    para_.PublishSubH                   = 0;
    para_.PublishStartPosition          = 0;
    para_.PublishConfPosition           = 1;
    para_.PublishArrow                  = 1;
    para_.PublishConfNum                = 1;
    para_.PublishFusedConfArrow         = 0;
    para_.PublishRedundantConfArrow     = 0;
    para_.PublishWeakConf               = 0;
    para_.PublishCondConf               = 0;
    para_.PublishHotspotsInOutTSPAngles = 0;
    para_.PublishReferenceDistPoint     = 0;
    para_.PublishGroundTruth            = 0;
    para_.PublishSensingCoverage        = 0;
    para_.PublishRedundantConfArrow     = 0;
    para_.PublishFusedConfArrow         = 0;
    para_.PublishAdaptivePlannedConf    = 0;
    para_.PublishExploration            = 0;
    para_.ReconstructionFile            = 'accurate_gdm.mat';
    %para_.PublishTextWidth              = 'two-columns';
    para_.PublishFigureTitle            = 'One step exploration plan for gas tomography';

    % data to publish
    % -------------------------------------------------------------
    dataPublish.FoV                     = FoV;
    dataPublish.FoVfaces                = FoVfaces;    
    dataPublish.SensorRange             = SensorRange;
    dataPublish.robotStartCellInMap     = origin_conf;
    dataPublish.ConfNum                 = conf_ind;
    dataPublish.ConfOrn                 = conf_ord;
    dataPublish.confSequenceNum         = 1:numel(executed_confs_ind);
    dataPublish.hotspots                = hcs_ind;
    dataPublish.map_env                 = map_env;
    dataPublish.publish_gdm_for_executed_conf_num = publish_gdm_for_executed_conf_num
    
    % publish
    %--------------------------------------------------------------
    fPublish_GDMPlan( dataPublish,para_,dir_ )
        
    % save plots
    %--------------------------------------------------------------
    file_name = sprintf('%s-1t-ARMEx-step3',para_.ExperimentTitle);
    %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    %print('-painters','-dpdf','-r500',[dir_.TomographyResults,file_name]); % pdf
    %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png
    export_fig([dir_.str,file_name], '-pdf','-r500')
    %pause
    
end


%***************************************************************
% 
%                PUBLISH - EXECUTED PLAN
%
%***************************************************************
if publish_executed_plan == 1
        
    fprintf(1,['\n\n',...
           '*********************************************\n',...
           '        PUBLISH 1T-RMS EXECUTED PLAN \n',...
           '*********************************************\n\n'])
    
    % load preprocessing of first center of mass for initialization.
    %load([dir_.Solutions,'detection_initial_preprocess.mat'],...
    %    'SensorRange','FoVfaces');
    load([dir_.sppdetection,'preprocess_detection.mat'],'SensorRange','FoVfaces');
    
    hcs_sub = load([dir_.str,'hc_list_final.txt']);
    hcs_ind = sub2ind(size(map_env),hcs_sub(:,1),hcs_sub(:,2));    
    executed_confs = load([dir_.str,'executed_confs_final.txt']);    
    executed_confs_ind = sub2ind(size(map_env),...
                                 round(executed_confs(:,1)-0.5),...
                                 round(executed_confs(:,2)-0.5) );
    executed_confs_orn = executed_confs(:,3);
    
    % publish parameters
    %--------------------------------------------------------------
    para_.PublishFoV                    = 0;  % draw field of view.
    para_.PublishSymbolicFoV            = 1;
    para_.PublishCandConf               = 0;
    para_.PublishGDM                    = 1;
    para_.PublishHotspots               = 1;
    para_.PublishSubH                   = 0;
    para_.PublishStartPosition          = 1;
    para_.PublishConfPosition           = 1;
    para_.PublishArrow                  = 1;
    para_.PublishConfNum                = 1;
    para_.PublishFusedConfArrow         = 0;
    para_.PublishRedundantConfArrow     = 0;
    para_.PublishWeakConf               = 0;
    para_.PublishCondConf               = 0;
    para_.PublishHotspotsInOutTSPAngles = 0;
    para_.PublishReferenceDistPoint     = 0;
    para_.PublishGroundTruth            = 0;
    para_.PublishSensingCoverage        = 0;
    para_.PublishRedundantConfArrow     = 0;
    para_.PublishFusedConfArrow         = 0;
    para_.PublishAdaptivePlannedConf    = 0;
    para_.PublishExploration            = 0;
    para_.ReconstructionFile            = 'accurate_gdm.mat';
    %para_.PublishTextWidth              = 'two-columns';
    para_.PublishFigureTitle            = 'One step exploration plan for gas tomography';

    % data to publish
    % -------------------------------------------------------------
    dataPublish.FoV                     = FoV;
    dataPublish.FoVfaces                = FoVfaces;    
    dataPublish.SensorRange             = SensorRange;
    dataPublish.robotStartCellInMap     = origin_conf;
    dataPublish.ConfNum                 = executed_confs_ind;
    dataPublish.ConfOrn                 = executed_confs_orn;
    dataPublish.confSequenceNum         = 1:numel(executed_confs_ind);
    %dataPublish.AdaptivePlannedConf_ind = plannedConfs_io(:,1);
    %dataPublish.AdaptivePlannedConf_orn = plannedConfs_io(:,2);
    %dataPublish.AdaptivePlannedConf_seq = 1:numel(plannedConfs_io(:,1));
    %dataPublish.AdaptiveHotspots        = hotspots_executed;
    dataPublish.hotspots                = hcs_ind;
    dataPublish.map_env                 = map_env;
    
    % directories
    %--------------------------------------------------------------
    %dir_.logs  = dir_.TomographyLogs;
    %dir_.recon = dir_.Reconstruction;
        
    % publish
    %--------------------------------------------------------------
    fPublish_GDMPlan( dataPublish,para_,dir_ )
        
    % save plots
    %--------------------------------------------------------------
    file_name = sprintf('%s-one-step-exploration-plan',para_.ExperimentTitle);
    %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    %print('-painters','-dpdf','-r500',[dir_.TomographyResults,file_name]); % pdf
    %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png
    export_fig([dir_.str,file_name], '-pdf','-r500')
end




%***************************************************************
% 
%                  PUBLISH - VIDEO
%
%***************************************************************
if publish_video == 1
        
    fprintf(1,['\n\n',...
           '*********************************************\n',...
           '      PUBLISH 1t-ARMEx EXPLORATION VIDEO \n',...
           '*********************************************\n\n'])
    
    % publish
    % ---------------------
    fVideo_ARMEx( para_,dir_ )
end


%***************************************************************
% 
%                PUBLISH - FINAL RECONSTRUCTION
% 
%***************************************************************
if publish_reconstruction == 1
    fprintf(1,['\n\n',...
           '*********************************************\n',...
           '     PUBLISH 1T-RMS FINAL RECONSTRUCTION \n',...
           '*********************************************\n\n'])
    % load preprocess maps
    load([dir_.sppdetection,'preprocess_detection.mat'],'SensorRange','FoVfaces');
    
    % publish parameters
    % -------------------------------------------------------------
    para_.PublishFoV                    = 0;  % draw field of view.
    para_.PublishSymbolicFoV            = 0;
    para_.PublishCandConf               = 0;
    para_.PublishGDM                    = 1;
    para_.PublishHotspots               = 0;
    para_.PublishSubH                   = 0;
    para_.PublishStartPosition          = 0;
    para_.PublishConfPosition           = 0;
    para_.PublishArrow                  = 0;
    para_.PublishConfNum                = 0;
    para_.PublishWeakConf               = 0;
    para_.PublishCondConf               = 0;
    para_.PublishHotspotsInOutTSPAngles = 0;
    para_.PublishReferenceDistPoint     = 0;
    para_.PublishGroundTruth            = 0;
    para_.PublishSensingCoverage        = 0;
    para_.PublishRedundantConfArrow     = 0;
    para_.PublishFusedConfArrow         = 0;
    para_.PublishAdaptivePlannedConf    = 0;
    para_.PublishExploration            = 0;
    para_.ReconstructionFile            = 'accurate_gdm.mat';
    %para_.PublishTextWidth              = 'two-columns';
    para_.PublishFigureTitle            = 'Tomographic reconstruction';
        
    % directories
    % -------------------------------------------------------------    
    %dir_.recon = dir_.Reconstruction;

    % data to publish
    % -------------------------------------------------------------
    dataPublish.FoVfaces            = FoVfaces;
    dataPublish.map_env             = map_env;
    
    % publish
    % ---------------------------------------------------
    fPublish_GDMPlan( dataPublish,para_,dir_ )
    % h=gcf;
    % set(h,'PaperPositionMode','auto');         
    % set(h,'PaperOrientation','landscape');
    % set(h,'Position',[50 50 1200 800]);
    % print(gcf, '-dpdf', 'test1.pdf')
    
    % save plots
    % -------------------------------------------------------------    
    file_name = sprintf('%s-one-step-exploration-reconstruction',para_.ExperimentTitle);
    %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    %print('-painters','-dpdf','-r500',[dir_.TomographyResults,file_name]); % pdf
    %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png    
    export_fig([dir_.str,file_name], '-pdf','-r500')
end


%***************************************************************
% 
%                   PUBLISH - FINAL SENSING COVERAGE
% 
%***************************************************************
if publish_sensing_coverage == 1
    fprintf(1,['\n\n',...
           '*********************************************\n',...
           '    PUBLISH 1T-RMS FINAL SENSING COVERAGE \n',...
           '*********************************************\n\n'])
    % load preprocess maps    
    load([dir_.Solutions,'detection_initial_preprocess.mat'],...
        'SensorRange','FoVfaces');
    
    % publish parameters
    % -------------------------------------------------------------
    para_.PublishFoV                    = 0;  % draw field of view.
    para_.PublishSymbolicFoV            = 0;
    para_.PublishCandConf               = 0;
    para_.PublishGDM                    = 0;
    para_.PublishHotspots               = 0;
    para_.PublishSubH                   = 0;
    para_.PublishStartPosition          = 0;
    para_.PublishConfPosition           = 0;
    para_.PublishArrow                  = 0;
    para_.PublishConfNum                = 0;
    para_.PublishWeakConf               = 0;
    para_.PublishCondConf               = 0;
    para_.PublishHotspotsInOutTSPAngles = 0;
    para_.PublishReferenceDistPoint     = 0;
    para_.PublishGroundTruth            = 0;
    para_.PublishSensingCoverage        = 1;
    para_.PublishRedundantConfArrow     = 0;
    para_.PublishFusedConfArrow         = 0;
    para_.PublishAdaptivePlannedConf    = 0;
    para_.PublishExploration            = 0;
    para_.ReconstructionFile            = 'reconstruction_OneStepExploration.mat';
    %para_.PublishTextWidth              = 'two-columns';
    para_.PublishFigureTitle            = 'Sensing Coverage';
    
    % directories
    % -------------------------------------------------------------
    dir_.recon = dir_.Reconstruction;

    % data to publish
    % -------------------------------------------------------------
    dataPublish.FoVfaces            = FoVfaces;
    dataPublish.map_env             = map_env;
    
    % publish
    % ---------------------------------------------------
    fPublish_GDMPlan( dataPublish,para_,dir_ )
    % h=gcf;
    % set(h,'PaperPositionMode','auto');         
    % set(h,'PaperOrientation','landscape');
    % set(h,'Position',[50 50 1200 800]);
    % print(gcf, '-dpdf', 'test1.pdf')
    
    % save plots
    % -------------------------------------------------------------    
    file_name = sprintf('%s-one-step-exploration-sensing-coverage',para_.ExperimentTitle);
    %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    %print('-painters','-dpdf','-r500',[dir_.TomographyResults,file_name]); % pdf
    %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png
    export_fig([dir_.Solutions,file_name], '-pdf','-r500')
end

%---- parameters
save([dir_.str,'expriment_paramets.mat'],'para_');



diary save
diary off

end


function fHumanExplorationJFR( para_ )
%=====================================================================
%   HUMAN EXPERT EXPLORATION
%=====================================================================
% 




solve_human_exploration   = 0;
solve_reconstruction      = 0;
solve_coverage            = 0;
solve_traveling_path_dist = 0;
solve_evaluation          = 1;

publish_exploration       = 0;
publish_reconstruction    = 0;
publish_sensing_coverage  = 0;




%=============================
%---- directories
%=============================
[ dir_ ] = fDirHumanJFR( para_ );
dir_.Recon = dir_.Solutions;
dir_.recon = dir_.Solutions;
dir_.logs = dir_.MeasurementLogs;


%============================
%   DIARY
%============================
DateTimeNow = datestr(now,'yyyy-mm-dd-HH:MM:SS');
diaryfilename = sprintf('diary_%s',DateTimeNow);
diary([dir_.Solutions,diaryfilename]);

%-- plan type
para_.PlanType = 'human-exploration';
dir_.ROSLogsThis = sprintf('%s%s/',dir_.ROSLogs,para_.PlanType);

para_.SamplingType = 'adaptive';
para_.PublishUpperPPM = 10000;


%************************************************
%               MAP:
%************************************************
file__map_env  = sprintf('%s_map_coverage.dat',para_.environment);
map_env        = load([dir_.env,file__map_env]);

file__cell_env = sprintf('%s_cellsize_coverage.dat',para_.environment);
cellsize_env   = load([dir_.env,file__cell_env]);



%--------------------------
%   SENSING PARAMETERS
%--------------------------
SensorRange.m = para_.SensingRangeM;
SensorRange.cell = SensorRange.m/cellsize_env;
FoV = para_.FieldOfViewDEG;


%*************************************************
% 
%              EXPLORATION
% 
%*************************************************
if solve_human_exploration == 1
    
    initialConfNum = 11; %23; %13;%13;

    
    fHumanExpertExploration( initialConfNum,dir_,para_ )
end




%*************************************************
% 
%              RECONSTRUCTION
% 
%*************************************************
if solve_reconstruction == 1 
    
    fprintf(1,['\n\n',...
           '*********************************************\n',...
           '               Reconstrction \n',...
           '*********************************************\n\n'])
    
    
    %dir_.logs = dir_.MeasurementLogs; 
    %dir_.Recon = dir_.Reconstruction;
    visualize = 0;
    
    [ M,...
      mean_map,...
      cell_coords_x,...
      cell_coords_y ] = fReconstructionFixed( para_,...
                                              dir_,...
                                              visualize );
                                               
    save([dir_.Solutions,'reconstruction_HumanExploration.mat'],...
        'M','mean_map','cell_coords_x','cell_coords_y');
    
    
end


%*************************************************
% 
%              COVERAGE
% 
%*************************************************
if solve_coverage == 1 
    
    fprintf(1,['\n\n',...
           '*********************************************\n',...
           '               SENSING COVERAGE \n',...
           '*********************************************\n\n'])
       
    % para_.ConfType = 'Planned'; %'Executed'
    para_.ConfType = 'Executed';
    
    sensingCoverage = fCalculateCoverage( map_env,SensorRange,FoV,para_,dir_ );
    save([dir_.Solutions,'sensing_coverage.mat'],'sensingCoverage');
    
end



%*************************************************
% 
%              TRAVELING PATH/DISTANCE
% 
%*************************************************
if solve_traveling_path_dist == 1 
    
    fprintf(1,['\n\n',...
           '*********************************************\n',...
           '         TRAVELING PATH/DISTANCE \n',...
           '*********************************************\n\n'])
       
    para_.ConfType = 'Executed';
    
    [ travelingPath,travelingDist ] = fTravelingPathDist( map_env,para_,dir_ );
    save([dir_.Solutions,'traveling_path_dist.mat'],...
        'travelingPath','travelingDist');
    
end





%*************************************************
% 
%              EVALUATION
% 
%*************************************************
if solve_evaluation == 1
    
    fprintf(1,['\n\n',...
           '*********************************************\n',...
           '                EVALUATION \n',...
           '*********************************************\n\n'])
    
    visualize = 1;
    recon_filename = 'reconstruction_HumanExploration.mat';
    
    dir_.recon = dir_.Solutions;
    
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
      f_measure ] = fEvaluation3( recon_filename,...
                                  map_env,...
                                  cellsize_env,...
                                  para_,...
                                  dir_,...
                                  visualize )
                                   
    save([dir_.Evaluation,'evaluation.mat'],...
        'sc_true',...
        'sc_estimated',...
        'true_positives',...
        'false_positives',...      
        'num_of_true_positives',...
        'num_of_false_positives',...
        'true_positives_dist',...
        'false_positives_dist',...
        'dist_estimatedSources2trueSources',...
        'precision',...
        'recall',...
        'f_measure');
end










%*************************************************
% 
%              PUBLISH - EXPLORATION
% 
%*************************************************
if publish_exploration == 1
    
        
    fprintf(1,'\n\n--> Publish Human Exploration \n')
    
    
end


%*************************************************
% 
%              PUBLISH - RECONSTRUCTION
% 
%*************************************************
if publish_reconstruction == 1
    
    fprintf(1,'\n\n--> Publish Tomographic Reconstruction - Human \n')

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
    para_.ReconstructionFile            = 'reconstruction_HumanExploration.mat';
    para_.PublishFigureTitle            = 'Tomographic reconstruction - Human Expert';
    
    
    % directories
    % -------------------------------------------------------------
    %dir_.logs  = dir_.TomographyLogs;
    %dir_.recon = dir_.Reconstruction;

    % data to publish
    % -------------------------------------------------------------
    %dataPublish.FoVfaces            = FoVfaces;
    dataPublish.map_env             = map_env;
    
    % publish
    % ---------------------------------------------------
    fPublish_GDMPlan( dataPublish,para_,dir_ )
    
    % save plots
    % -------------------------------------------------------------
    file_name = sprintf('%s-Human-exploration-reconstruction',para_.ExperimentTitle);
    %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    %print('-painters','-dpdf','-r500',[dir_.TomographyResults,file_name]); % pdf
    %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png
    export_fig([dir_.Solutions,file_name], '-pdf','-r500')
end


%*************************************************
% 
%              PUBLISH - SENSING COVERAGE
% 
%*************************************************
if publish_sensing_coverage == 1
    
        
    fprintf(1,'\n\n--> Publish Sensing Coverage - Human \n')
    
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
    para_.ReconstructionFile            = 'reconstruction_HumanExploration.mat';
    para_.PublishFigureTitle            = 'Sensing Coverage';
    
    % data to publish
    % -------------------------------------------------------------
    dataPublish.map_env             = map_env;
    
    % publish
    % ---------------------------------------------------
    fPublish_GDMPlan( dataPublish,para_,dir_ )

   
    % save plots
    % -------------------------------------------------------------    
    file_name = sprintf('%s-one-step-exploration-sensing-coverage',para_.ExperimentTitle);
    %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    %print('-painters','-dpdf','-r500',[dir_.TomographyResults,file_name]); % pdf
    %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png
    export_fig([dir_.Solutions,file_name], '-pdf','-r500')

end


diary save
diary off

end


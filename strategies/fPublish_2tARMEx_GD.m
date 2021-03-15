function fPublish_2tARMEx_GD( publish_gd,para_,dir_ )


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

% dataYAML = ReadYamlRaw([dir_.env,para_.EnvironmentYamlFile]);
% dataYAML_origin = cell2mat(dataYAML.origin);

% robot_origin_x = abs(dataYAML_origin(1)/dataYAML.resolution);
% robot_origin_y = abs(dataYAML_origin(2)/dataYAML.resolution);
% map_pix_size   = dataYAML.resolution;

% robot_origin_in_map_x = 1;
% robot_origin_in_map_y = 1;
% dataYAML.resolution          = 1;
% -----------------------------------------------------------------------------------



fprintf(1,['\n\n',...
           '*********************************************************************\n\n',...
           '                    PUBLISH GAS DETECTION \n\n',...
           '*********************************************************************\n\n'])

       
%************************************************
%               MAP AND REALTED:
%************************************************
[ map_env,...
  map_conf,...
  origin_env,...
  origin_conf,...
  cellsize_env,...
  cellsize_conf,...
  cellsize_org ] = fMapStuff( para_,dir_ );


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                              PUBLISH - SPP DETECTION
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if publish_gd.spp
    
    fprintf(1,'\n\n--> Publish SPP Detection \n')
    
    dir_.logs = dir_.gdplan_logs;    
    load([dir_.gdplan_spp,'preprocess_detection.mat'],...
            'oV',...
            'SensorRange',...
            'FoVfaces',...
            'confKept');    

        dir_.gdplan_spp
    load([dir_.gdplan_spp,'spp_detection_solution.mat'],'C')

    Conf = zeros(size(oV,2),1);
    Conf(confKept) = C;
    
    load([dir_.gdplan_spp,'tsp_detection_solution.mat']);

    % publish parameters
    % -------------------------------------------------------------
    para_.PublishFoV                    = 1;  % draw field of view.
    para_.PublishSymbolicFoV            = 1;
    para_.PublishCandConf               = 0;
    para_.PublishGDM                    = 0;
    para_.PublishHotspots               = 0;
    para_.PublishSubH                   = 0;
    para_.PublishStartPosition          = 0;
    para_.PublishConfPosition           = 1;
    para_.PublishArrow                  = 1;
    para_.PublishConfNum                = 0; 
    para_.PublishWeakConf               = 0;
    para_.PublishCondConf               = 0;
    para_.PublishHotspotsInOutTSPAngles = 0;
    para_.PublishReferenceDistPoint     = 0;
    para_.PublishGroundTruth            = 0;
    para_.PublishSensingCoverage        = 0;
    para_.PublishDetectionEvents        = 0;
    %para_.PublishTextWidth              = 'two-columns';
    para_.PublishFigureTitle            = 'An exploration plan for gas detection';

    % data to publish
    % -------------------------------------------------------------
    dataPublish.FoV                 = FoVragd;
    dataPublish.FoVfaces            = FoVfaces;
    dataPublish.Conf                = Conf;
    dataPublish.SensorRange         = SensorRange;
    dataPublish.confSequenceNum     = confSequenceNum;
    dataPublish.robotStartCellInMap = origin_conf;
    dataPublish.map_env             = map_env;
    
    % publish
    % -------------------------------------------------------------
    fPublish_GD( dataPublish,para_,dir_ )
    
    % save plots
    % -------------------------------------------------------------
    %pause
    file_name = sprintf('%s-detection-plan',para_.ExperimentTitle);
    %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    %print('-painters','-dpdf','-r500',[dir_.DetectionResults,file_name]); % pdf
    %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png    
    export_fig([dir_.gdplan,file_name], '-pdf','-r500')
    export_fig([dir_.gdplan,file_name], '-png','-r500')



end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                       PUBLISH - RECONSTRUCTION (DETECTION)
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if publish_gd.recon
    
    fprintf(1,'\n\n--> Publish Reconstruction Detection \n')
    
    load([dir_.gdplan_spp,'preprocess_detection.mat'],'FoVfaces')
    dir_.logs = dir_.gdplan_logs;
    
    
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
    para_.PublishDetectionEvents        = 0;
    para_.ReconstructionFile            = 'coarse_gdm.mat';
    %para_.PublishTextWidth              = 'two-columns';
    para_.PublishFigureTitle            = 'Coarse map';
    
    % data to publish
    % -------------------------------------------------------------
    dataPublish.FoVfaces            = FoVfaces;
    dataPublish.map_env             = map_env;


    % publish
    % ---------------------------------------------------
    fPublish_GD( dataPublish,para_,dir_ )

    % h=gcf;
    % set(h,'PaperPositionMode','auto');         
    % set(h,'PaperOrientation','landscape');
    % set(h,'Position',[50 50 1200 800]);
    % print(gcf, '-dpdf', 'test1.pdf')
    
    % save plots
    % -------------------------------------------------------------
    %pause
    file_name = sprintf('%s-detection-reconstruction',para_.ExperimentTitle);
    %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    %print('-painters','-dpdf','-r500',[dir_.DetectionResults,file_name]); % pdf
    %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png    
    export_fig([dir_.gdplan,file_name], '-pdf','-r500')
    export_fig([dir_.gdplan,file_name], '-png','-r500')

end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                       PUBLISH - SENSING COVERAGE (DETECTION)
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if publish_gd.coverage
    
    fprintf(1,'\n\n--> Publish Gas Detection Coverage \n')
    
    load([dir_.gdplan_spp,'preprocess_detection.mat'],'FoVfaces');
    dir_.logs = dir_.gdplan_logs;
    
    meas_files = 'measurement_conf*.dat*';
    
    
    
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
    para_.PublishDetectionEvents        = 0;
    para_.PublishHighConcentrationTH    = 7000;
    %para_.PublishTextWidth              = 'two-columns';
    para_.MeasurementFiles              = 'measurement_conf*.dat*';
    para_.PublishFigureTitle            = 'Sensing coverage for gas detection';
    
    % data to publish
    % -------------------------------------------------------------
    dataPublish.FoVfaces            = FoVfaces;
    dataPublish.meas_files          = meas_files;
    dataPublish.map_env             = map_env;

    % publish
    % -------------------------------------------------------------
    fPublish_GD( dataPublish,para_,dir_ )
    
    % save plots
    % -------------------------------------------------------------
    %pause
    file_name = sprintf('%s-detection-coverage',para_.ExperimentTitle);
    %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    %print('-painters','-dpdf','-r500',[dir_.DetectionResults,file_name]); % pdf
    %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png    
    export_fig([dir_.gdplan,file_name], '-pdf','-r500')
    export_fig([dir_.gdplan,file_name], '-png','-r500')
    
    
    
end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                       PUBLISH - GAS DETECTION EVENTS
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if publish_gd.detection_events
    
    
    fprintf(1,'\n\n--> Publish Gas Detection Events \n')
    
    dir_.logs = dir_.gdplan_logs;    
    load([dir_.gdplan_spp,'preprocess_detection.mat'],...
            'oV',...
            'SensorRange',...
            'FoVfaces',...
            'confKept',...
            'map_env',...
            'robotStartCellInMap');    
        
    load([dir_.gdplan_spp,'spp_detection_solution.mat'],'C')

    Conf = zeros(size(oV,2),1);
    Conf(confKept) = C;
    
    load([dir_.gdplan_spp,'tsp_detection_solution.mat']);
        
    dir_.logs = dir_.gdplan_logs;
    
    meas_files = 'measurement_conf*.dat*';
    
    
    
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
    para_.PublishGroundTruth            = 1;
    para_.PublishSensingCoverage        = 0;
    para_.PublishDetectionEvents        = 1;
    para_.PublishHighConcentrationTH    = 7000;
    %para_.PublishTextWidth              = 'two-columns';
    para_.MeasurementFiles              = 'measurement_conf*.dat*';
    para_.PublishFigureTitle            = 'Gas detection events';
    
    % data to publish
    % -------------------------------------------------------------
    dataPublish.FoVfaces            = FoVfaces;
    dataPublish.meas_files          = meas_files;
    dataPublish.map_env             = map_env;
    
    dataPublish.FoV                 = FoVragd;
    dataPublish.Conf                = Conf;
    dataPublish.SensorRange         = SensorRange;
    dataPublish.confSequenceNum     = confSequenceNum;

    
    % publish
    % -------------------------------------------------------------
    fPublish_GD( dataPublish,para_,dir_ )
    
    % save plots
    % -------------------------------------------------------------
    %pause
    file_name = sprintf('%s-detection-events',para_.ExperimentTitle);
    %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    %print('-painters','-dpdf','-r500',[dir_.DetectionResults,file_name]); % pdf
    %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png    
    export_fig([dir_.gdplan,file_name], '-pdf','-r500')
    export_fig([dir_.gdplan,file_name], '-png','-r500')    
end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                        PUBLISH - GROUND TRUTH ON ENVIRONMENT MAP
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if publish_gd.gt
    
    fprintf(1,'\n\n--> Ground Truth \n')
    
    load([dir_.gdplan_spp,'preprocess_detection.mat'],'FoVfaces');
    dir_.logs = dir_.gt;

    
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
    para_.PublishGroundTruth            = 1;
    para_.PublishSensingCoverage        = 0;
    para_.PublishDetectionEvents        = 0;
    %para_.PublishTextWidth              = 'two-columns';
    para_.PublishFigureTitle            = 'Ground truth (environment and gas distributions)';
    
    % data to publish
    % -------------------------------------------------------------
    dataPublish.FoVfaces                = FoVfaces;
    dataPublish.map_env                 = map_env;
    
    % publish
    % -------------------------------------------------------------
    fPublish_GD( dataPublish,para_,dir_ )
    
    
    % save plots
    % -------------------------------------------------------------
    %set(h,'Units','Inches');
    %pos = get(h,'Position');
    %set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(2), pos(4)])
    %pause
    %file_name = sprintf('gt_%s',para_.GroundTruthFile);
    file_name = sprintf('%s-ground-truth',para_.ExperimentTitle);
    %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    %print('-painters','-dpdf','-r500',[dir_.GroundTruthLogs,file_name]); % pdf
    %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png    
    export_fig([dir_.gt,file_name], '-pdf','-r500')
    export_fig([dir_.gt,file_name], '-png','-r500')
end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                        PUBLISH - ENVIRONMENT MAP
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if publish_gd.env
    
    fprintf(1,'\n\n--> Environment Map \n')
    
    %load([dir_.DetectionResults,'preprocess_detection.mat'],'FoVfaces','map_env');
    %dir_.logs = dir_.GroundTruthLogs;

    
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
    para_.PublishSensingCoverage        = 0;
    para_.PublishDetectionEvents        = 0;
    %para_.PublishTextWidth              = 'two-columns';
    para_.PublishFigureTitle            = 'Environment map';
    
    % data to publish
    % -------------------------------------------------------------
    %dataPublish.FoVfaces            = FoVfaces;
    dataPublish.map_env             = map_env;
    
    % publish
    % -------------------------------------------------------------
    fPublish_GD( dataPublish,para_,dir_ )
    
    
    % save plots
    % -------------------------------------------------------------
    %set(h,'Units','Inches');
    %pos = get(h,'Position');
    %set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(2), pos(4)])
    %pause
    %file_name = sprintf('gt_%s',para_.GroundTruthFile);
    file_name = sprintf('%s-environment',para_.ExperimentTitle);
    %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    %print('-painters','-dpdf','-r500',[dir_.GroundTruthLogs,file_name]); % pdf
    %print('-painters','-dpng','-r400',[dir_.GroundTruthLogs,file_name]); % png    
    export_fig([dir_.gt,file_name], '-pdf','-r500')
    export_fig([dir_.gt,file_name], '-png','-r500')
    

end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                              PUBLISH - RWL1 ITERATIONS
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if publish_gd.rwl1_iterations
    
    fprintf(1,'\n\n--> Publish RWL1 iterations \n')
    
    dir_.logs = dir_.gdplan_logs;    
    load([dir_.gdplan_spp,'preprocess_detection.mat'],...
            'oV',...
            'SensorRange',...
            'FoVfaces',...
            'confKept');    
        
    %load([dir_.DetectionResults,'spp_detection_solution.mat'],'C')
    %Conf = zeros(size(oV,2),1);
    %Conf(confKept) = C;
    
    %load([dir_.DetectionResults,'tsp_detection_solution.mat']);

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
    para_.PublishSensingCoverage        = 0;
    para_.PublishDetectionEvents        = 0;
    %para_.PublishTextWidth              = 'two-columns';
    para_.PublishFigureTitle            = 'RWL1 iterations';

    % data to publish
    % -------------------------------------------------------------
    %dataPublish.FoV                 = FoVragd;
    dataPublish.FoVfaces            = FoVfaces;
    %dataPublish.Conf                = Conf;
    %dataPublish.SensorRange         = SensorRange;
    %dataPublish.confSequenceNum     = confSequenceNum;
    %dataPublish.robotStartCellInMap = origin_conf;
    dataPublish.map_env             = map_env;
    
    % publish
    % -------------------------------------------------------------
    fPublish_RWL1Iterations( dataPublish,para_,dir_ )
    
    % save plots
    % -------------------------------------------------------------
    %pause
    % file_name = sprintf('%s-rwl1-iteration',para_.ExperimentTitle);
    % %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    % %print('-painters','-dpdf','-r500',[dir_.DetectionResults,file_name]); % pdf
    % %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png    
    % export_fig([dir_.DetectionResults,file_name], '-pdf','-r500')
    % export_fig([dir_.DetectionResults,file_name], '-png','-r500')
    
    
    


end



end

function fPublish_ICRA16_GDM( publish_tomography,para_,dir_ ) 



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
% dataYAML.resolution          = 1;
% -----------------------------------------------------------------------------------

fprintf(1,['\n\n',...
           '*********************************************************************\n\n',...
           '                     PUBLISH GAS TOMOGRAPHY (ICRA16) \n\n',...
           '*********************************************************************\n\n'])

%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%              PUBLISH ENVIRONMENT, GAS DISTRIBUTION MAP AND HOTSPOTS
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if publish_tomography.maps
    
    fprintf(1,'\n\n******************** Publish Maps *************************\n')
    
    dir_.logs = dir_.DetectionLogs;
    load([dir_.TomographyResults,'preprocess_maps.mat']);

    % publish parameters
    % ---------------------------------------------------
    para_.PublishFoV                    = 0;  % draw field of view.
    para_.PublishSymbolicFoV            = 0;
    para_.PublishCandConf               = 0;
    para_.PublishGDM                    = 1;
    para_.PublishHotspots               = 1;
    para_.PublishSubH                   = 0;
    para_.PublishStartPosition          = 1;
    para_.PublishConfPosition           = 0;
    para_.PublishArrow                  = 0;
    para_.PublishConfNum                = 0;
    para_.PublishWeakConf               = 0;
    para_.PublishCondConf               = 0;
    para_.PublishHotspotsInOutTSPAngles = 1;
    para_.PublishReferenceDistPoint     = 0;
    para_.PublishGroundTruth            = 0;
    para_.PublishSensingCoverage        = 0;
    para_.PublishRedundantConfArrow     = 0;
    para_.PublishFusedConfArrow         = 0;
    para_.ReconstructionFile            = 'reconstruction_detection.mat';
    para_.ReconstructionLocation        = dir_.DetectionResults;
    para_.PublishTextWidth              = 'single-column';
    
    % directories
    % -------------------------------------------------------------
    dir_.logs  = dir_.TomographyLogs;
    dir_.recon = dir_.DetectionResults;
    
    
    % data to publish
    % ---------------------------------------------------
    dataPublish.FoVfaces            = FoVfaces;
    dataPublish.hotspots            = hotspotSEQ;
    dataPublish.InOutTSPAnglesMean  = InOutTSPAnglesMean;
    dataPublish.robotStartCellInMap = robotStartCellInMap;


    fPublish_GDMPlanICRA16( dataPublish,para_,dir_ )
    
    pause(2)
    file_name = sprintf('maps');
    print('-painters','-dpdf','-r500',[dir_.TomographyResults,file_name]); % pdf

end         






%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                               LOCAL SOLUTIONS
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if publish_tomography.spp_hotspots
    
    fprintf(1,'\n\n****************** Publish Local Solutions ********************\n')
    
    
    load([dir_.TomographyResults,'preprocess_maps.mat'],...
        'hotspotSEQ',...
        'map_recon',...
        'InOutTSPAnglesMean',...
        'robotStartCellInMap');

    % publish subproblems
    % -----------------------------------------------------------------------
    for i = 1:numel(hotspotSEQ)

        % load preprocessing 
        load([dir_.TomographyLocalSol,'preprocess_XVD_hotspot',num2str(i),'.mat']);

        % load SPP solution
        selectedConf = load([dir_.TomographyLocalSol,...
            'selectedConf_hotspot',num2str(i),'.dat']);
        selectedConf = round(selectedConf);

        % conf vector with all cells in the map.
        Conf = zeros(size(oV,2),1);
        Conf(confKept) = selectedConf;
        
        % publish parameters
        % ---------------------------------------------------
        para_.PublishFoV                    = 0;  % draw field of view.
        para_.PublishSymbolicFoV            = 0;
        para_.PublishCandConf               = 1;
        para_.PublishGDM                    = 1;
        para_.PublishHotspots               = 1;
        para_.PublishSubH                   = 0;
        para_.PublishStartPosition          = 0;
        para_.PublishConfPosition           = 1;
        para_.PublishArrow                  = 1;
        para_.PublishConfNum                = 0;
        para_.PublishWeakConf               = 0;
        para_.PublishCondConf               = 0;
        para_.PublishHotspotsInOutTSPAngles = 0;
        para_.PublishReferenceDistPoint     = 1;
        para_.PublishGroundTruth            = 0;
        para_.PublishSensingCoverage        = 0;        
        para_.PublishRedundantConfArrow     = 0;
        para_.PublishFusedConfArrow         = 0;
        para_.ReconstructionFile            = 'reconstruction_detection.mat';
        para_.ReconstructionLocation        = dir_.DetectionResults;
        para_.PublishTextWidth              = 'single-column';
        
        % directories
        % -------------------------------------------------------------
        dir_.logs  = dir_.TomographyLogs;
        dir_.recon = dir_.DetectionResults;

         
        % data to publish
        % ---------------------------------------------------
        dataPublish = [];
        dataPublish.map_candConf        = map_candConf;
        dataPublish.map_env             = map_env;
        % dataPublish.map_coverage        = map_coverage;
        dataPublish.map_gd              = map_recon;
        dataPublish.FoV                 = FoV;
        dataPublish.FoVfaces            = FoVfaces;
        dataPublish.Conf                = Conf;
        dataPublish.SensorRange         = SensorRange;
        dataPublish.hotspots            = hotspotSEQ(i);
        % dataPublish.subHotspots         = subHotspots;
        % dataPublish.confSequenceNum     = confSequenceNum;
        % dataPublish.weakConfGlobal_ind  = RedundantConf_ind;
        % dataPublish.condensedConf_ind   = condensedConf_ind;
        % dataPublish.InOutTSPAnglesMean  = InOutTSPAnglesMean;
        dataPublish.refDistPoint        = refDistPoint;
        % dataPublish.robotStartCellInMap = robotStartCellInMap;

        % publish
        % ---------------------------------------------------
        fPublish_GDMPlanICRA16( dataPublish,para_,dir_ )
        
        
        pause(2)
        file_name = sprintf('spp_hotspot%02d',i);
        %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
        print('-painters','-dpdf','-r500',[dir_.TomographyLocalSol,file_name]); % pdf
        %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png
        
        pause(2)

    end

end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                               PUBLISH SPP HOTSPOTS (ALL)
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if publish_tomography.spp_hotspots_all

    fprintf(1,'\n\n************** Publish Local Solutions (all) ******************\n')
    
    dir_.logs = dir_.DetectionLogs;

    load([dir_.TomographyResults,'preprocess_maps.mat'],...
        'map_env',...
        'hotspotSEQ',...
        'map_recon',...
        'InOutTSPAnglesMean',...
        'robotStartCellInMap');
 
    % load preprocessing of first center of mass for initialization.
    load([dir_.TomographyLocalSol,'preprocess_XVD_hotspot',num2str(1),'.mat'],...
        'oV','map_candConf','map_coverage');

    % initialize conf selected and candidate conf map.
    Conf = zeros(size(oV,2),1);
    MAP_candConf = zeros(size(map_candConf));

    ConfAllHotspots = zeros(size(oV,2),numel(hotspotSEQ));
    refDistPoints = zeros(numel(hotspotSEQ),2);
    
    % for each center of mass
    for num = 1:numel(hotspotSEQ)

        % load preprocessing
        load([dir_.TomographyLocalSol,'preprocess_XVD_hotspot',num2str(num),'.mat']);

        MAP_candConf = MAP_candConf+map_candConf;

        % load SPP solution for the center of mass
        selectedConf = load([dir_.TomographyLocalSol,...
            'selectedConf_hotspot',num2str(num),'.dat']);

        selectedConf = round(selectedConf);

        % conf vector with all cells in the map.
        conf = zeros(size(oV,2),1);
        conf(confKept) = selectedConf;

        %find(conf)

        % update overall conf vector.
        Conf = Conf+conf;

        % update conf vector for all hotspots
        ConfAllHotspots(:,num) = conf;
        
        % reference point for distance
        refDistPoints(num,:) = refDistPoint;

    end
    
    
    % publish parameters
    % ---------------------------------------------------
    para_.PublishFoV                    = 0;  % draw field of view.
    para_.PublishSymbolicFoV            = 0;
    para_.PublishCandConf               = 1;
    para_.PublishGDM                    = 1;
    para_.PublishHotspots               = 1;
    para_.PublishSubH                   = 0;
    para_.PublishStartPosition          = 0;
    para_.PublishConfPosition           = 1;
    para_.PublishArrow                  = 1;
    para_.PublishConfNum                = 0;
    para_.PublishWeakConf               = 0;
    para_.PublishCondConf               = 0;
    para_.PublishHotspotsInOutTSPAngles = 0;
    para_.PublishReferenceDistPoint     = 1;
    para_.PublishGroundTruth            = 0;
    para_.PublishSensingCoverage        = 0;    
    para_.PublishRedundantConfArrow     = 0;
    para_.PublishFusedConfArrow         = 0;
    para_.ReconstructionFile            = 'reconstruction_detection.mat';
    para_.ReconstructionLocation        = dir_.DetectionResults;
    para_.PublishTextWidth              = 'single-column';
    
    % directories
    % -------------------------------------------------------------
    dir_.logs  = dir_.TomographyLogs;
    dir_.recon = dir_.DetectionResults;
    
    
    % data to publish
    % ---------------------------------------------------
    dataPublish.map_candConf        = MAP_candConf;
    dataPublish.FoV                 = FoV;
    dataPublish.FoVfaces            = FoVfaces;
    dataPublish.Conf                = Conf;
    dataPublish.SensorRange         = SensorRange;
    dataPublish.hotspots            = hotspotSEQ;
    dataPublish.refDistPoint        = refDistPoints;
    
        
    % publish
    % ---------------------------------------------------
    fPublish_GDMPlanICRA16( dataPublish,para_,dir_ )
    
    pause(2)
    file_name = sprintf('spp_hotspots_all');
    %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    print('-painters','-dpdf','-r500',[dir_.TomographyLocalSol,file_name]); % pdf
    %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png

    
    

end

%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                        PUBLISH LOCAL AND FUSED SOLUTIONS
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if publish_tomography.spp_fusion
    
    fprintf(1,'\n\n****************** Publish Fusion ********************\n')
    
    
    dir_.logs = dir_.DetectionLogs;
    
    % load preprocess maps    
    load([dir_.TomographyResults,'preprocess_maps.mat'],...
        'hotspotSEQ',...
        'map_env',...
        'map_coverage',...
        'map_recon',...
        'InOutTSPAnglesMean',...
        'robotStartCellInMap'); 
    % load preprocessing of first center of mass for initialization.
    load([dir_.TomographyLocalSol,'preprocess_XVD_hotspot',num2str(1),'.mat'],...
        'SensorRange',...
        'FoVfaces');
    
   
    % --- local solutions
    % load preprocessing of first center of mass for initialization.
    load([dir_.TomographyLocalSol,'preprocess_XVD_hotspot',num2str(1),'.mat'],...
        'oV','map_candConf','map_coverage');

    % initialize conf selected and candidate conf map.
    Conf = zeros(size(oV,2),1);
    
    refDistPoints = zeros(numel(hotspotSEQ),2);

    ConfAllHotspots = zeros(size(oV,2),numel(hotspotSEQ));
    % for each center of mass
    for num = 1:numel(hotspotSEQ)

        % load preprocessing
        load([dir_.TomographyLocalSol,'preprocess_XVD_hotspot',num2str(num),'.mat']);

        % load SPP solution for the center of mass
        selectedConf = load([dir_.TomographyLocalSol,...
            'selectedConf_hotspot',num2str(num),'.dat']);
        selectedConf = round(selectedConf);

        % conf vector with all cells in the map.
        conf = zeros(size(oV,2),1);
        conf(confKept) = selectedConf;
        
        % update overall conf vector.
        Conf = Conf+conf;
        
        % update conf vector for all hotspots
        ConfAllHotspots(:,num) = conf;
        
        % reference point for distance
        refDistPoints(num,:) = refDistPoint;

    end
    
    % -- local solutions are redundant conf
    RedundantConf_ind = find(Conf); 
    
    load([dir_.TomographyFusion,'integrated_solution_hotspots_fusions.mat'],'Conf');
    
    
    % publish parameters
    % ---------------------------------------------------
    para_.PublishFoV                    = 0;  % draw field of view.
    para_.PublishSymbolicFoV            = 1;
    para_.PublishCandConf               = 0;
    para_.PublishGDM                    = 1;
    para_.PublishHotspots               = 1;
    para_.PublishSubH                   = 0;
    para_.PublishStartPosition          = 0;
    para_.PublishConfPosition           = 1;
    para_.PublishArrow                  = 1;
    para_.PublishConfNum                = 0;
    para_.PublishFusedConfArrow         = 0;
    para_.PublishRedundantConfArrow     = 0;
    para_.PublishWeakConf               = 1;
    para_.PublishCondConf               = 0;
    para_.PublishHotspotsInOutTSPAngles = 0;
    para_.PublishReferenceDistPoint     = 0;
    para_.PublishGroundTruth            = 0;
    para_.PublishSensingCoverage        = 0;
    para_.PublishRedundantConfArrow     = 0;
    para_.PublishFusedConfArrow         = 0;
    para_.ReconstructionFile            = 'reconstruction_detection.mat';
    para_.ReconstructionLocation        = dir_.DetectionResults;
    para_.PublishTextWidth              = 'single-column';
    
    % directories
    % -------------------------------------------------------------
    dir_.logs  = dir_.TomographyLogs;
    dir_.recon = dir_.DetectionResults;
    
    % data to publish
    % ---------------------------------------------------
    dataPublish.FoV                 = FoV;
    dataPublish.FoVfaces            = FoVfaces;
    dataPublish.Conf                = Conf;
    dataPublish.SensorRange         = SensorRange;
    dataPublish.hotspots            = hotspotSEQ;
    dataPublish.refDistPoint        = refDistPoints;    
    dataPublish.weakConfGlobal_ind  = RedundantConf_ind;
    
    % publish
    % ---------------------------------------------------
    pause(2)
    fPublish_GDMPlanICRA16( dataPublish,para_,dir_ )


    % save plots
    % -------------------------------------------------------------
    pause(1)
    file_name = sprintf('spp_fusion_tomography');
    print('-painters','-dpdf','-r500',[dir_.TomographyFusion,file_name]); % pdf
    %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png
    



end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                        PUBLISH SPP TOMOGRAPHY SOLUTION
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if publish_tomography.spp_solution
    
    
    fprintf(1,'\n\n****************** Publish SPP Tomography ********************\n')
    
    dir_.logs = dir_.DetectionLogs;
    
    % load preprocess maps    
    load([dir_.TomographyResults,'preprocess_maps.mat'],...
        'hotspotSEQ',...
        'map_env',...
        'map_coverage',...
        'map_recon',...
        'InOutTSPAnglesMean',...
        'robotStartCellInMap'); 
    % load preprocessing of first center of mass for initialization.
    load([dir_.TomographyLocalSol,'preprocess_XVD_hotspot',num2str(1),'.mat'],...
        'SensorRange',...
        'FoVfaces');
    
    load([dir_.TomographyFusion,'integrated_solution_hotspots_fusions.mat'],'Conf');

    load([dir_.TomographyResults,'tsp_tomography_spp.mat'],'confSequenceNum');

    %{
    tCost/10
    transStep = 0;
    rotatStep = 0;
    f = 4;
    for l = 1:numel(TRoute)-1
        iCELL = fix(TRoute(l)/f)+(~~mod(TRoute(l),f));
        jCELL = fix(TRoute(l+1)/f)+(~~mod(TRoute(l+1),f));
        if iCELL == jCELL
            rotatStep = rotatStep+1;
        else
            transStep = transStep+1;
        end
    end
    transStep 
    rotatStep
    %}
    
    % publish parameters
    % ---------------------------------------------------
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
    para_.ReconstructionFile            = 'reconstruction_detection.mat';
    para_.ReconstructionLocation        = dir_.DetectionResults;
    para_.PublishTextWidth              = 'single-column';
    
    % directories
    % -------------------------------------------------------------
    dir_.logs  = dir_.TomographyLogs;
    dir_.recon = dir_.DetectionResults;
    
    % data to publish
    % ---------------------------------------------------
    dataPublish.FoV                 = FoV;
    dataPublish.FoVfaces            = FoVfaces;
    dataPublish.Conf                = Conf;
    dataPublish.SensorRange         = SensorRange;
    dataPublish.hotspots            = hotspotSEQ;
    %dataPublish.refDistPoint        = refDistPoint;
    dataPublish.confSequenceNum     = confSequenceNum;
    dataPublish.robotStartCellInMap = robotStartCellInMap;
    

    % publish
    % ---------------------------------------------------
    pause(2)
    fPublish_GDMPlanICRA16( dataPublish,para_,dir_ )


    % save plots
    % -------------------------------------------------------------
    pause(1)
    file_name = sprintf('robot_plan_tomography');
    print('-painters','-dpdf','-r500',[dir_.TomographyResults,file_name]); % pdf
    %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png
    



end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%              PUBLISH - COVERAGE MEASUREMENTS TOMOGRAPHY (FIXED)
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if publish_tomography.coverage_spp_fixed
    
    fprintf(1,'\n\n********* Publish Tomographic Coverage (Fixed) **************\n')
    
    % load preprocessing of first center of mass for initialization.
    load([dir_.TomographyLocalSol,'preprocess_XVD_hotspot',num2str(1),'.mat'],...
        'SensorRange','FoVfaces');
             
    dir_.logs = dir_.TomographyLogs;

    % publish parameters
    % -------------------------------------------------------------
    para_.PublishFoV                    = 0; % draw field of view.
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
    para_.PublishSensingCoverage        = 1;
    para_.PublishRedundantConfArrow     = 0;
    para_.PublishFusedConfArrow         = 0;
    para_.ReconstructionFile            = 'reconstruction_detection.mat';    
    para_.MeasurementFiles              = 'measurement_conf*_fixed.dat*';
    para_.PublishTextWidth              = 'single-column';

    % data to publish
    % -------------------------------------------------------------
    dataPublish.FoVfaces            = FoVfaces;
    
    % directories
    % -------------------------------------------------------------
    dir_.logs  = dir_.TomographyLogs;
    dir_.recon = dir_.DetectionResults;
    
    % publish
    % ---------------------------------------------------
    fPublish_GDMPlanICRA16( dataPublish,para_,dir_ )
    
    
    % save plots
    % -------------------------------------------------------------    
    file_name = sprintf('coverage_tomography_fixed');
    %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    print('-painters','-dpdf','-r500',[dir_.TomographyResults,file_name]); % pdf
    %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png

end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%              PUBLISH - FIXED RECONSTRUCTION (TOMOGRAPHY)
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if publish_tomography.recon_spp_fixed
    
    fprintf(1,'\n\n********* Publish Reconstruction (Fixed) **************\n')
    
    % load preprocess maps    
    % load([dir_.SolutionLoc,'preprocess_maps.mat'],'map_env'); 

    % load preprocessing of first center of mass for initialization.
    load([dir_.TomographyLocalSol,'preprocess_XVD_hotspot',num2str(1),'.mat'],...
        'SensorRange','FoVfaces');
             
    dir_.logs = dir_.TomographyLogs;
    %map_gd = mean_map;

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
    para_.ReconstructionFile            = 'reconstruction_tomography_fixed.mat';
    para_.PublishTextWidth              = 'single-column';
    
    % directories
    % -------------------------------------------------------------
    dir_.logs  = dir_.TomographyLogs;
    dir_.recon = dir_.TomographyResults;

    % data to publish
    % -------------------------------------------------------------
    dataPublish.FoVfaces            = FoVfaces;

    % publish
    % ---------------------------------------------------
    fPublish_GDMPlanICRA16( dataPublish,para_,dir_ )

    % h=gcf;
    % set(h,'PaperPositionMode','auto');         
    % set(h,'PaperOrientation','landscape');
    % set(h,'Position',[50 50 1200 800]);
    % print(gcf, '-dpdf', 'test1.pdf')
    
    % save plots
    % -------------------------------------------------------------
    
    file_name = sprintf('reconstruction_tomography_fixed');
    %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    print('-painters','-dpdf','-r500',[dir_.TomographyResults,file_name]); % pdf
    %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png

end



end

function fPublish_2tARMEx_GDM( publish_gdm,para_,dir_ ) 




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
'               PUBLISH GAS DISTRIBUTION MAPPING \n\n',...
'*********************************************************************\n\n'])

       
%**********************************************************************
% MAP:
%**********************************************************************
file__map_env   = sprintf('%s_map_coverage.dat',para_.environment);
file__map_conf  = sprintf('%s_map_conf.dat',para_.environment);

file__org_env   = sprintf('%s_origin_coverage.dat',para_.environment);
file__org_conf  = sprintf('%s_origin_conf.dat',para_.environment);

file__cell_env  = sprintf('%s_cellsize_coverage.dat',para_.environment);
file__cell_conf = sprintf('%s_cellsize_conf.dat',para_.environment);
file__cell_org  = sprintf('%s_cellsize_original.dat',para_.environment);

file__cut_org   = sprintf('%s_mapcut_original.dat',para_.environment);

map_env         = load([dir_.env,file__map_env]);
map_conf        = load([dir_.env,file__map_conf]);
origin_env      = load([dir_.env,file__org_env]);
origin_conf     = load([dir_.env,file__org_conf]);
cellsize_env    = load([dir_.env,file__cell_env]);
cellsize_conf   = load([dir_.env,file__cell_conf]);
cellsize_org    = load([dir_.env,file__cell_org]);
% mapCut          = load([dir_.env,file__cut_org]);
       
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                               PUBLISH GROUND TRUTH
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if publish_gdm.ground_truth
            
    fprintf(1,'\n\n--> Publish Maps \n')
    
    load([dir_.TomographyMaps,'preprocess_maps.mat']);

    % publish parameters
    % ---------------------------------------------------
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
    para_.PublishRedundantConfArrow     = 0;
    para_.PublishFusedConfArrow         = 0;
    para_.PublishAdaptivePlannedConf    = 0;
    para_.PublishExploration            = 0;
    %para_.ReconstructionFile            = 'reconstruction_detection.mat';
    %para_.PublishTextWidth              = 'two-columns';
    para_.PublishFigureTitle            = 'Ground truth';
    
    % directories
    % -------------------------------------------------------------    
    dir_.logs  = dir_.DetectionLogs;
    dir_.recon = dir_.DetectionResults;
    
    % data to publish
    % ---------------------------------------------------
    dataPublish.FoVfaces            = FoVfaces;
    dataPublish.hotspots            = hotspotSEQ;
    dataPublish.InOutTSPAnglesMean  = InOutTSPAnglesMean;
    dataPublish.map_env             = map_env;
    
    % publish
    % ---------------------------------------------------
    fPublish_GDMPlan( dataPublish,para_,dir_ )
    
    pause(2)
        
    file_name = sprintf('%s-ground-truth',para_.ExperimentTitle);
    %print('-painters','-dpdf','-r500',[dir_.TomographyResults,file_name]); % pdf    
    export_fig([dir_.TomographyResults,file_name], '-pdf','-r500')
    export_fig([dir_.TomographyResults,file_name], '-png','-r500')
    
end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                               PUBLISH ENVIRONMENTS
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if publish_gdm.maps
            
    fprintf(1,'\n\n--> Publish Maps \n')
    
    load([dir_.TomographyMaps,'preprocess_maps.mat']);

    % publish parameters
    % ---------------------------------------------------
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
    para_.ReconstructionFile            = 'reconstruction_detection.mat';
    %para_.PublishTextWidth              = 'two-columns';
    para_.PublishFigureTitle            = 'A coarse map with hotspots';
    
    % directories
    % -------------------------------------------------------------    
    dir_.logs  = dir_.DetectionLogs;
    dir_.recon = dir_.DetectionResults;
    
    % data to publish
    % ---------------------------------------------------
    dataPublish.FoVfaces            = FoVfaces;
    dataPublish.hotspots            = hotspotSEQ;
    dataPublish.InOutTSPAnglesMean  = InOutTSPAnglesMean;
    dataPublish.map_env             = map_env;
    
    % publish
    % ---------------------------------------------------
    fPublish_GDMPlan( dataPublish,para_,dir_ )
    
    set(gcf,'color','w');

    pause(2)
        
    file_name = sprintf('%s-maps',para_.ExperimentTitle);
    %print('-painters','-dpdf','-r500',[dir_.TomographyResults,file_name]); % pdf    
    export_fig([dir_.TomographyResults,file_name], '-pdf','-r500')
    export_fig([dir_.TomographyResults,file_name], '-png','-r500')
    
end




%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                               LOCAL SOLUTIONS
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if publish_gdm.spp_hotspots
    
    fprintf(1,'\n\n--> Publish Local Solutions \n')
    
    load([dir_.TomographyMaps,'preprocess_maps.mat'],...
        'hotspotSEQ',...
        'InOutTSPAnglesMean');

    % publish subproblems
    % -----------------------------------------------------------------------
    for hc_num = 1:numel(hotspotSEQ)

        % load preprocessing 
        filename = sprintf('local_hc%02d_xvd.mat',hc_num);
        load([dir_.TomographyLocalSol,filename]);

        % load SPP solution
        filename = sprintf('local_hc%02d_spp.dat',hc_num);
        selectedConf = load([dir_.TomographyLocalSol,filename]);
        selectedConf = round(selectedConf);

        % conf vector with all cells in the map.
        Conf = zeros(size(oV,2),1);
        Conf(confKept) = selectedConf;
        
        % -- conf numbers
        selectedConf_ind = confKept(selectedConf==1);
        
        % -- orientations
        selectedConf_orn = ConfOrns(selectedConf==1);
        
        
        figtitle = sprintf('Local solution for gas tomography (hc #%02d)',hc_num);
        
        % publish parameters
        % ---------------------------------------------------
        para_.PublishFoV                    = 0;  % draw field of view.
        para_.PublishSymbolicFoV            = 1;
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
        para_.PublishAdaptivePlannedConf    = 0;
        para_.PublishExploration            = 0;
        para_.ReconstructionFile            = 'reconstruction_detection.mat';
        %para_.PublishTextWidth              = 'two-columns';
        para_.PublishFigureTitle            = figtitle;
        
        % directories
        % -------------------------------------------------------------        
        %dir_.logs  = dir_.TomographyLogs;
        dir_.logs  = dir_.DetectionLogs;
        dir_.recon = dir_.DetectionResults;

        % data to publish
        % ---------------------------------------------------
        dataPublish.map_candConf        = map_candConf;
        dataPublish.FoV                 = FoV;
        dataPublish.FoVfaces            = FoVfaces;
        dataPublish.Conf                = Conf;
        dataPublish.SensorRange         = SensorRange;
        dataPublish.hotspots            = hotspotSEQ(hc_num);
        dataPublish.refDistPoint        = refDistPoint;
        dataPublish.ConfNum             = selectedConf_ind;
        dataPublish.ConfOrn             = selectedConf_orn;
        dataPublish.map_env             = map_env;

        % publish
        % ---------------------------------------------------
        fPublish_GDMPlan( dataPublish,para_,dir_ )
        
        pause(2)
        file_name = sprintf('%s-local-solution%02d',para_.ExperimentTitle,hc_num);
        %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
        %print('-painters','-dpdf','-r500',[dir_.TomographyResults,file_name]); % pdf
        %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png
        export_fig([dir_.TomographyResults,file_name], '-pdf','-r500')
        export_fig([dir_.TomographyResults,file_name], '-png','-r500')
        
        %{
        matlab2tikz([file_name,'.tex'],...
            'height','\figureheight',...
            'width','\figurewidth',...    
            'extraAxisOptions','ylabel near ticks');
        %}
        
        pause(2)

    end
    

end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                       PUBLISH SPP HOTSPOTS (ALL)
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if publish_gdm.spp_hotspots_all
        
    fprintf(1,'\n\n--> Publish Local Solutions (all) \n')

    load([dir_.TomographyMaps,'preprocess_maps.mat'],...
        'map_env',...
        'hotspotSEQ',...
        'map_recon',...
        'InOutTSPAnglesMean');
 
    % load preprocessing of first center of mass for initialization.
    %load([dir_.TomographyLocalSol,'preprocess_XVD_hotspot',num2str(1),'.mat'],...
    %    'oV','map_candConf','map_coverage');
    load([dir_.TomographyLocalSol,sprintf('local_hc%02d_xvd.mat',1)],...
        'oV',...
        'map_candConf',...
        'map_coverage');
    
    % initialize conf selected and candidate conf map.
    Conf = zeros(size(oV,2),1);
    MAP_candConf = zeros(size(map_candConf));
    
    selectedConf_ind = [];
    selectedConf_orn = [];

    ConfAllHotspots = zeros(size(oV,2),numel(hotspotSEQ));
    
    refDistPoints = zeros(numel(hotspotSEQ),2);
    
    % for each center of mass
    for num = 1:numel(hotspotSEQ)

        % load preprocessing
        %load([dir_.TomographyLocalSol,'preprocess_XVD_hotspot',num2str(num),'.mat']);
        load([dir_.TomographyLocalSol,sprintf('local_hc%02d_xvd.mat',num)])

        MAP_candConf = MAP_candConf+map_candConf;

        % load SPP solution for the center of mass
        %selectedConf = load([dir_.TomographyLocalSol,...
        %    'selectedConf_hotspot',num2str(num),'.dat']);
        selectedConf = load([dir_.TomographyLocalSol,...
            sprintf('local_hc%02d_spp.dat',num)]);
        selectedConf = round(selectedConf);

        % conf vector with all cells in the map.
        conf = zeros(size(oV,2),1);
        conf(confKept) = selectedConf;
        
        
        % -- conf numbers
        selectedConf_ind_this = confKept(selectedConf==1);
        selectedConf_ind = [selectedConf_ind,selectedConf_ind_this];
        
        % -- orientations
        selectedConf_Orn_this = ConfOrns(selectedConf==1);
        selectedConf_orn = [selectedConf_orn,selectedConf_Orn_this];

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
    para_.PublishCandConf               = 0;
    para_.PublishGDM                    = 1;
    para_.PublishHotspots               = 1;
    para_.PublishSubH                   = 0;
    para_.PublishStartPosition          = 0;
    para_.PublishConfPosition           = 0;
    para_.PublishArrow                  = 0;
    para_.PublishConfNum                = 0;
    para_.PublishWeakConf               = 1;
    para_.PublishCondConf               = 0;
    para_.PublishHotspotsInOutTSPAngles = 0;
    para_.PublishReferenceDistPoint     = 0;
    para_.PublishGroundTruth            = 0;
    para_.PublishSensingCoverage        = 0;
    para_.PublishRedundantConfArrow     = 0;
    para_.PublishFusedConfArrow         = 0;
    para_.PublishAdaptivePlannedConf    = 0;
    para_.PublishExploration            = 0;
    para_.ReconstructionFile            = 'reconstruction_detection.mat';
    %para_.PublishTextWidth              = 'two-columns';
    para_.PublishFigureTitle            = 'Local solutions for gas tomography';
    
    % directories
    % -------------------------------------------------------------
    %dir_.logs  = dir_.TomographyLogs;
    dir_.logs  = dir_.DetectionLogs;
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
    dataPublish.ConfNum             = selectedConf_ind;
    dataPublish.ConfOrn             = selectedConf_orn;
    dataPublish.map_env             = map_env;
    
    
    %--- for Achim presenstation
    dataPublish.RedundantConfNum    = selectedConf_ind;
    dataPublish.RedundantConfOrn    = selectedConf_orn;
    
    %dataPublish.confSequenceNum     = confSequenceNum;
    %dataPublish.weakConfGlobal_ind  = [];
    dataPublish.robotStartCellInMap = origin_conf;
    %dataPublish.map_env             = map_env;
    %------------------
    
    % publish
    % ---------------------------------------------------
    fPublish_GDMPlan( dataPublish,para_,dir_ )

    pause(2)
    file_name = sprintf('%s-local-solutions',para_.ExperimentTitle);
    %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    %print('-painters','-dpdf','-r500',[dir_.TomographyResults,file_name]); % pdf
    %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png    
    export_fig([dir_.TomographyResults,file_name], '-pdf','-r500')
    export_fig([dir_.TomographyResults,file_name], '-png','-r500')
    
end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                        PUBLISH LOCAL AND FUSED SOLUTIONS
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if publish_gdm.spp_fusion
    
    fprintf(1,'\n\n--> Publish Fusion of Local Solutions \n')
    % load preprocess maps    
    load([dir_.TomographyMaps,'preprocess_maps.mat'],...
        'hotspotSEQ',...
        'map_env',...
        'map_coverage',...
        'map_recon',...
        'InOutTSPAnglesMean'); 
    % load preprocessing of first center of mass for initialization.    
    load([dir_.TomographyLocalSol,sprintf('local_hc%02d_xvd.mat',1)],...
        'SensorRange',...
        'FoVfaces','oV','map_candConf','map_coverage');
    % load([dir_.SolutionLoc,'integrated_solution_hotspots_fusions_fusions.mat'])
    load([dir_.TomographyFusion,'spp_fusion.mat'],...
        'tspConf',...
        'selectedConfFusion_ioh');
    
    % initialize conf selected and candidate conf map.
    Conf = zeros(size(oV,2),1);
        
    selectedConf_ind = [];
    selectedConf_orn = [];

    ConfAllHotspots = zeros(size(oV,2),numel(hotspotSEQ));
    % for each center of mass
    for num = 1:numel(hotspotSEQ)

        % load preprocessing
        %load([dir_.TomographyLocalSol,'preprocess_XVD_hotspot',num2str(num),'.mat']);
        load([dir_.TomographyLocalSol,sprintf('local_hc%02d_xvd.mat',num)])

        % load SPP solution for the center of mass
        selectedConf = load([dir_.TomographyLocalSol,...
            sprintf('local_hc%02d_spp.dat',num)]);
        selectedConf = round(selectedConf);

        % conf vector with all cells in the map.
        conf = zeros(size(oV,2),1);
        conf(confKept) = selectedConf;
                
        % -- conf numbers
        selectedConf_ind_this = confKept(selectedConf==1);        
        selectedConf_ind = [selectedConf_ind,selectedConf_ind_this];
        
        % -- orientations
        selectedConf_Orn_this = ConfOrns(selectedConf==1);        
        selectedConf_orn = [selectedConf_orn,selectedConf_Orn_this];
        
        % update overall conf vector.
        Conf = Conf+conf;
        
        % update conf vector for all hotspots
        ConfAllHotspots(:,num) = conf;
    end
    
    
    
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
    para_.PublishAdaptivePlannedConf    = 0;
    para_.PublishExploration            = 0;
    para_.ReconstructionFile            = 'reconstruction_detection.mat';
    %para_.PublishTextWidth              = 'two-columns';
    para_.PublishFigureTitle            = 'Fusion of the local solutions for gas tomography';
    
    % directories
    % -------------------------------------------------------------
    %dir_.logs  = dir_.TomographyLogs;
    dir_.logs  = dir_.DetectionLogs;
    dir_.recon = dir_.DetectionResults;

    % data to publish
    % ---------------------------------------------------
    dataPublish.map_candConf        = [];
    dataPublish.FoV                 = FoV;
    dataPublish.FoVfaces            = FoVfaces;
    dataPublish.Conf                = [];
    dataPublish.SensorRange         = SensorRange;
    dataPublish.hotspots            = hotspotSEQ;
    %dataPublish.refDistPoint        = refDistPoint;
    dataPublish.ConfNum             = selectedConfFusion_ioh(:,1);
    dataPublish.ConfOrn             = selectedConfFusion_ioh(:,2);
    dataPublish.FusedConfNum        = [];
    dataPublish.FusedConfOrn        = [];
    dataPublish.RedundantConfNum    = selectedConf_ind;
    dataPublish.RedundantConfOrn    = selectedConf_orn;    
    dataPublish.robotStartCellInMap = origin_conf;
    dataPublish.map_env             = map_env;
    
    % publish
    % ---------------------------------------------------
    pause(2)
    fPublish_GDMPlan( dataPublish,para_,dir_ )


    % save plots
    % -------------------------------------------------------------
    pause(1)
    file_name = sprintf('%s-fusion',para_.ExperimentTitle);
    %print('-painters','-dpdf','-r500',[dir_.TomographyResults,file_name]); % pdf
    %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png    
    export_fig([dir_.TomographyResults,file_name], '-pdf','-r500')
    export_fig([dir_.TomographyResults,file_name], '-png','-r500')
    
end



%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                        PUBLISH SPP TOMOGRAPHY SOLUTION
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if publish_gdm.spp_solution
    
    fprintf(1,'\n\n--> Publish Tomographic Plan (initial) \n')
          
    % load preprocess maps    
    load([dir_.gdmplan_spp,'preprocess_maps.mat'],...
        'hotspotSEQ',...
        'map_env',...
        'map_coverage',...
        'map_recon',...
        'InOutTSPAnglesMean'); 
    % load preprocessing of first center of mass for initialization.
    %load([dir_.TomographyLocalSol,'preprocess_XVD_hotspot',num2str(1),'.mat'],...
    load([dir_.TomographyLocalSol,sprintf('local_hc%02d_xvd.mat',1)],...
        'SensorRange',...
        'FoVfaces');
    % load([dir_.SolutionLoc,'integrated_solution_hotspots_fusions_fusions.mat'])
    load([dir_.TomographyFusion,'spp_fusion.mat'],...
        'tspConf',...
        'selectedConfFusion_ioh');

    load([dir_.TomographyResults,'tsp_tomography_spp.mat']);
    
    
    %{
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
    
    
    % --- local solutions
    % initialize conf selected and candidate conf map.
    selectedLocalConf_ind = [];
    selectedLocalConf_orn = [];

    % for each center of mass
    for num = 1:numel(hotspotSEQ)

        % load preprocessing
        %load([dir_.TomographyLocalSol,'preprocess_XVD_hotspot',num2str(num),'.mat']);
        load([dir_.TomographyLocalSol,sprintf('local_hc%02d_xvd.mat',num)])

        % load SPP solution for the center of mass
        selectedLocalConf = load([dir_.TomographyLocalSol,...
            sprintf('local_hc%02d_spp.dat',num)]);
        selectedLocalConf = round(selectedLocalConf);
        
        % -- conf numbers
        selectedLocalConf_ind_this = confKept(selectedLocalConf==1);        
        selectedLocalConf_ind = [selectedLocalConf_ind,selectedLocalConf_ind_this];
        
        % -- orientations
        selectedLocalConf_Orn_this = ConfOrns(selectedLocalConf==1);        
        selectedLocalConf_orn = [selectedLocalConf_orn,selectedLocalConf_Orn_this];
    end
    
    
    
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
    para_.PublishWeakConf               = 0; %1
    para_.PublishCondConf               = 0;
    para_.PublishHotspotsInOutTSPAngles = 0;
    para_.PublishReferenceDistPoint     = 0;
    para_.PublishGroundTruth            = 0;
    para_.PublishSensingCoverage        = 0;
    para_.PublishRedundantConfArrow     = 0;
    para_.PublishFusedConfArrow         = 0;
    para_.PublishAdaptivePlannedConf    = 0;
    para_.PublishExploration            = 0;
    para_.ReconstructionFile            = 'reconstruction_detection.mat';
    %para_.PublishTextWidth              = 'two-columns';
    para_.PublishFigureTitle            = 'An initial exploration plan for gas tomography';
    
    % directories
    % -------------------------------------------------------------
    %dir_.logs  = dir_.TomographyLogs;
    dir_.logs  = dir_.DetectionLogs;
    dir_.recon = dir_.DetectionResults;

    % data to publish
    % ---------------------------------------------------
    dataPublish.FoV                 = FoV;
    dataPublish.FoVfaces            = FoVfaces;
    dataPublish.Conf                = [];
    dataPublish.SensorRange         = SensorRange;
    dataPublish.hotspots            = hotspotSEQ;
    %dataPublish.refDistPoint        = refDistPoint;
    dataPublish.ConfNum             = selectedConfFusion_ioh(:,1);
    dataPublish.ConfOrn             = selectedConfFusion_ioh(:,2);
    dataPublish.FusedConfNum        = [];
    dataPublish.FusedConfOrn        = [];
    dataPublish.RedundantConfNum    = selectedLocalConf_ind;
    dataPublish.RedundantConfOrn    = selectedLocalConf_orn;
    
    dataPublish.confSequenceNum     = confSequenceNum;
    dataPublish.weakConfGlobal_ind  = [];
    dataPublish.robotStartCellInMap = origin_conf;
    dataPublish.map_env             = map_env;
    

    % publish
    % ---------------------------------------------------
    pause(2)
    fPublish_GDMPlan( dataPublish,para_,dir_ )


    % save plots
    % -------------------------------------------------------------
    pause(1)
    file_name = sprintf('%s-tomographic-plan-initial',para_.ExperimentTitle);
    %print('-painters','-dpdf','-r500',[dir_.TomographyResults,file_name]); % pdf
    %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png    
    export_fig([dir_.TomographyResults,file_name], '-pdf','-r500')
    export_fig([dir_.TomographyResults,file_name], '-png','-r500')
    



end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                        PUBLISH TOMOGRAPHIC EXPLORATION
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if publish_gdm.exploration
    
    
    fprintf(1,'\n\n--> Publish Tomographic Exploration \n')
    
    % load preprocess maps    
    load([dir_.TomographyMaps,'preprocess_maps.mat'],...
        'hotspotSEQ',...
        'map_env',...
        'map_coverage',...
        'map_recon',...
        'InOutTSPAnglesMean'); 
    % load preprocessing of first center of mass for initialization.
    %load([dir_.TomographyLocalSol,'preprocess_XVD_hotspot',num2str(1),'.mat'],...
    load([dir_.TomographyLocalSol,sprintf('local_hc%02d_xvd.mat',1)],...
        'SensorRange',...
        'FoVfaces');
    % load([dir_.SolutionLoc,'integrated_solution_hotspots_fusions_fusions.mat'])
    load([dir_.TomographyFusion,'spp_fusion.mat'],...
        'tspConf',...
        'selectedConfFusion_ioh');

    load([dir_.TomographyResults,'tsp_tomography_spp.mat']);
    
    
    % publish parameters
    % ---------------------------------------------------
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
    para_.PublishExploration            = 1;
    para_.ReconstructionFile            = 'reconstruction_detection.mat';
    %para_.PublishTextWidth              = 'two-columns';    
    para_.MeasurementFiles              = 'measurement_conf*_adaptive.dat*';
    para_.PublishFigureTitle            = 'An exploration for gas tomography';
    
    % directories
    % -------------------------------------------------------------
    %dir_.logs  = dir_.TomographyLogs;
    dir_.logs  = dir_.DetectionLogs;
    dir_.recon = dir_.DetectionResults;

    % data to publish
    % ---------------------------------------------------
    dataPublish.FoV                 = FoV;
    dataPublish.FoVfaces            = FoVfaces;
    dataPublish.Conf                = [];
    dataPublish.SensorRange         = SensorRange;
    dataPublish.hotspots            = hotspotSEQ;
    %dataPublish.refDistPoint        = refDistPoint;
    dataPublish.ConfNum             = selectedConfFusion_ioh(:,1);
    dataPublish.ConfOrn             = selectedConfFusion_ioh(:,2);
    dataPublish.FusedConfNum        = [];
    dataPublish.FusedConfOrn        = [];
    dataPublish.RedundantConfNum    = [];
    dataPublish.RedundantConfOrn    = [];
    
    dataPublish.confSequenceNum     = confSequenceNum;
    dataPublish.weakConfGlobal_ind  = [];
    dataPublish.robotStartCellInMap = origin_conf;
    
    dataPublish.tRoute              = tRoute;
    dataPublish.TRoute              = tRoute;%TRoute;
    
    dataPublish.map_env             = map_env;

    % publish
    % ---------------------------------------------------
    pause(2)
    fPublish_GDMPlan( dataPublish,para_,dir_ )
    
    
    % save plots
    % -------------------------------------------------------------
    pause(1)
    file_name = sprintf('%s-exploration_tomography',para_.ExperimentTitle);
    %print('-painters','-dpdf','-r500',[dir_.TomographyResults,file_name]); % pdf
    %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png    
    export_fig([dir_.TomographyResults,file_name], '-pdf','-r500')
    export_fig([dir_.TomographyResults,file_name], '-png','-r500')
    

end



%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%              PUBLISH - ADAPTIVE MEASUREMENTS
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if publish_gdm.spp_plan_adaptive
    
        
    fprintf(1,'\n\n--> Publish Tomographic Plan (adaptive) \n')
        
    
    %load([dir_.TomographyMaps,'preprocess_maps.mat'],...    
    %    'robotStartCellInMap'); 
    
    % load preprocessing of first center of mass for initialization.
    %load([dir_.TomographyLocalSol,'preprocess_XVD_hotspot',num2str(1),'.mat'],...
    load([dir_.TomographyLocalSol,sprintf('local_hc%02d_xvd.mat',1)],...
        'SensorRange','FoVfaces');
    
    load([dir_.TomographyResults,'adaptive_replanning.mat'],...
        'plannedConfs_io',...
        'executedConfs_io',...
        'hotspots_executed');
    
             
    

    % publish parameters
    % -------------------------------------------------------------
    para_.PublishFoV                    = 0;  % draw field of view.
    para_.PublishSymbolicFoV            = 1;
    para_.PublishCandConf               = 0;
    para_.PublishGDM                    = 1;
    para_.PublishHotspots               = 0;
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
    para_.PublishAdaptivePlannedConf    = 1;
    para_.PublishExploration            = 0;
    para_.ReconstructionFile            = 'reconstruction_detection.mat';
    %para_.PublishTextWidth              = 'two-columns';
    para_.PublishFigureTitle            = 'An exploration for gas tomography with adaptive replanning';

    % data to publish
    % -------------------------------------------------------------
    dataPublish.FoV                     = FoV;
    dataPublish.FoVfaces                = FoVfaces;    
    dataPublish.SensorRange             = SensorRange;
    dataPublish.robotStartCellInMap     = origin_conf;
    dataPublish.ConfNum                 = executedConfs_io(:,1);
    dataPublish.ConfOrn                 = executedConfs_io(:,2);
    dataPublish.confSequenceNum         = 1:numel(executedConfs_io(:,1));
    dataPublish.AdaptivePlannedConf_ind = plannedConfs_io(:,1);
    dataPublish.AdaptivePlannedConf_orn = plannedConfs_io(:,2);
    dataPublish.AdaptivePlannedConf_seq = 1:numel(plannedConfs_io(:,1));
    dataPublish.AdaptiveHotspots        = hotspots_executed;
    dataPublish.map_env             = map_env;
    
    % directories
    % -------------------------------------------------------------
    %dir_.logs  = dir_.TomographyLogs;
    dir_.logs = dir_.TomographyLogsAdp;
    dir_.recon = dir_.DetectionResults;
    
    
    % publish
    % ---------------------------------------------------
    fPublish_GDMPlan( dataPublish,para_,dir_ )
    
    
    % save plots
    % -------------------------------------------------------------
    
    file_name = sprintf('%s-tomographic-plan-adaptive',para_.ExperimentTitle);
    %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    %print('-painters','-dpdf','-r500',[dir_.TomographyResults,file_name]); % pdf
    %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png    
    export_fig([dir_.TomographyResults,file_name], '-pdf','-r500')
    export_fig([dir_.TomographyResults,file_name], '-png','-r500')

end




%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%              PUBLISH - COVERAGE MEASUREMENTS TOMOGRAPHY (ADAPTIVE)
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
%{
if publish_tomography.coverage_spp_adaptive
    
        
    fprintf(1,'\n\n--> Publish Tomographic Coverage (adaptive) \n')
    
    
    % load preprocessing of first center of mass for initialization.
    %load([dir_.TomographyLocalSol,'preprocess_XVD_hotspot',num2str(1),'.mat'],...
    load([dir_.TomographyLocalSol,sprintf('local_hc%02d_xvd.mat',1)],...
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
    para_.PublishGroundTruth            = 1;
    para_.PublishSensingCoverage        = 1;
    para_.PublishRedundantConfArrow     = 0;
    para_.PublishFusedConfArrow         = 0;
    para_.PublishAdaptivePlannedConf    = 0;
    para_.PublishExploration            = 0;
    para_.ReconstructionFile            = 'reconstruction_detection.mat';    
    para_.MeasurementFiles              = 'measurement_conf*.dat';
    %para_.PublishTextWidth              = 'two-columns';
    para_.PublishFigureTitle            = 'Sensing coverage for gas tomography with adaptive replanning';

    % data to publish
    % -------------------------------------------------------------
    dataPublish.FoVfaces            = [];
    dataPublish.map_env             = map_env;
    
    % directories
    % -------------------------------------------------------------
    dir_.logs  = dir_.TomographyLogsAdp;
    dir_.recon = dir_.DetectionResults;
    
    % publish
    % ---------------------------------------------------
    fPublish_GDMPlan( dataPublish,para_,dir_ )
    
    
    % save plots
    % -------------------------------------------------------------
    
    file_name = sprintf('%s-tomographic-coverage-adaptive',para_.ExperimentTitle);
    %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    %print('-painters','-dpdf','-r500',[dir_.TomographyResults,file_name]); % pdf
    %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png
    
    export_fig([dir_.TomographyResults,file_name], '-pdf','-r500')

end
%}


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%              PUBLISH - COVERAGE MEASUREMENTS TOMOGRAPHY (FIXED)
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% if publish_tomography.coverage_spp_fixed
if publish_gdm.coverage
    
        
    fprintf(1,'\n\n--> Publish Tomographic Coverage (initial) \n')
    
    % load preprocessing of first center of mass for initialization.
    %load([dir_.TomographyLocalSol,'preprocess_XVD_hotspot',num2str(1),'.mat'],...
    load([dir_.TomographyLocalSol,sprintf('local_hc%02d_xvd.mat',1)],...
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
    para_.PublishGroundTruth            = 1;
    para_.PublishSensingCoverage        = 1;
    para_.PublishRedundantConfArrow     = 0;
    para_.PublishFusedConfArrow         = 0;
    para_.PublishAdaptivePlannedConf    = 0;
    para_.PublishExploration            = 0;
    para_.ReconstructionFile            = 'reconstruction_detection.mat';
    para_.MeasurementFiles              = 'measurement_conf*.dat';
    %para_.PublishTextWidth              = 'two-columns';
    para_.PublishFigureTitle            = 'Sensing coverage for gas tomography without adaptive replanning';
    

    % data to publish
    % -------------------------------------------------------------
    dataPublish.FoVfaces            = [];
    dataPublish.map_env             = map_env;
    
    % directories
    % -------------------------------------------------------------
    dir_.logs  = dir_.TomographyLogsFix;
    dir_.recon = dir_.DetectionResults;
    
    % publish
    % ---------------------------------------------------
    fPublish_GDMPlan( dataPublish,para_,dir_ )
    
    
    % save plots
    % -------------------------------------------------------------    
    file_name = sprintf('%s-tomographic-coverage-initial',para_.ExperimentTitle);
    %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    %print('-painters','-dpdf','-r500',[dir_.TomographyResults,file_name]); % pdf
    %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png    
    export_fig([dir_.TomographyResults,file_name], '-pdf','-r500')
    export_fig([dir_.TomographyResults,file_name], '-png','-r500')

end



%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%              PUBLISH - ADAPTIVE RECONSTRUCTION (TOMOGRAPHY)
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
%{
if publish_tomography.recon_spp_adaptive
    
        
    fprintf(1,'\n\n--> Publish Tomographic Reconstruction (adaptive) \n')
    
    % load preprocess maps    
    % load([dir_.SolutionLoc,'preprocess_maps.mat'],'map_env'); 

    % load preprocessing of first center of mass for initialization.
    %load([dir_.TomographyLocalSol,'preprocess_XVD_hotspot',num2str(1),'.mat'],...
    load([dir_.TomographyLocalSol,sprintf('local_hc%02d_xvd.mat',1)],...
        'SensorRange','FoVfaces');
             
    
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
    para_.PublishAdaptivePlannedConf    = 0;
    para_.PublishExploration            = 0;
    para_.ReconstructionFile            = 'reconstruction_OneStepExploration.mat';
    %para_.PublishTextWidth              = 'two-columns';
    para_.PublishFigureTitle            = 'Tomographic reconstruction with adaptive replanning';
    
    % directories
    % -------------------------------------------------------------
    dir_.logs  = dir_.TomographyLogsAdp;
    dir_.recon = dir_.TomographyReconAdp;

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
    
    file_name = sprintf('%s-tomographic-reconstruction-adaptive',para_.ExperimentTitle);
    %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    %print('-painters','-dpdf','-r500',[dir_.TomographyResults,file_name]); % pdf
    %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png
    
    export_fig([dir_.TomographyResults,file_name], '-pdf','-r500')

end
%}

%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%              PUBLISH - FIXED RECONSTRUCTION (TOMOGRAPHY)
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
%if publish_tomography.recon_spp_fixed
if publish_gdm.recon
    
        
    fprintf(1,'\n\n--> Publish Tomographic Reconstruction (initial) \n')
    
    % load preprocess maps    
    % load([dir_.SolutionLoc,'preprocess_maps.mat'],'map_env'); 

    % load preprocessing of first center of mass for initialization.
    %load([dir_.TomographyLocalSol,'preprocess_XVD_hotspot',num2str(1),'.mat'],...
    load([dir_.TomographyLocalSol,sprintf('local_hc%02d_xvd.mat',1)],...
        'SensorRange','FoVfaces');
             
    %dir_.logs = dir_.TomographyLogs;
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
    para_.PublishAdaptivePlannedConf    = 0;
    para_.PublishExploration            = 0;
    para_.ReconstructionFile            = 'reconstruction.mat';
    %para_.PublishTextWidth              = 'two-columns';
    para_.PublishFigureTitle            = 'Tomographic reconstruction without adaptive replanning';
    
    % directories
    % -------------------------------------------------------------
    dir_.logs  = dir_.TomographyLogsFix;
    dir_.recon = dir_.TomographyReconFix;

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
    
    file_name = sprintf('%s-tomographic-reconstruction-initial',para_.ExperimentTitle);
    %print('-painters','-dpdf','-r1049',[dir_.results,file_name]); % pdf
    %print('-painters','-dpdf','-r500',[dir_.TomographyResults,file_name]); % pdf
    %print('-painters','-dpng','-r400' ,[dir_.results,file_name]); % png    
    export_fig([dir_.TomographyResults,file_name], '-pdf','-r500')
    export_fig([dir_.TomographyResults,file_name], '-png','-r500')

end

end

function [ tspConf,...
           refPointDist,...
           selectedConfFusion_ioh,...
           tFusion ] = fFusionHotspotGeometries2t( FoV,...
                                                   SensorRange,...
                                                   cellsize_env,...
                                                   dir_,...
                                                   para_,...
                                                   visualize )
%
% 


tFusioni = tic; % initialize pre computational time.

% {

%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
%      load preprocess maps
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
fprintf(1,'-- load preprocess \n')

load([dir_.gdmplan_spp,'preprocess_maps.mat'],...
        'FoVfaces',...
        'hotspotSEQ',...
        'map_env',...
        'map_recon',...
        'map_coverage',...
        'OBSenv',...
        'OBSconf',...
        'SensorRange');

hcsSEQ_ind = hotspotSEQ;

%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
%    initial setup for fusion
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
fprintf(1,'-- initial setup for fusion \n')

% { 
[ localConfs,...
  localERGall,...
  localERQall,...
  localNumOfSelectedConfs,...
  refPointDist ] = fInitialFusionSetup2t( hcsSEQ_ind,...
                                          map_env,...
                                          dir_,...
                                          para_,...
                                          visualize );
save([dir_.gdmplan_spp,'initial_fusion_setup.mat'],...
        'localConfs',...
        'localERGall',...
        'localNumOfSelectedConfs',...
        'refPointDist');
%}
% load([dir_.TomographyFusion,'initial_fusion_setup.mat']);

%=======================================================================
% Publish Loal Solutions
%=======================================================================

if visualize == 1
    
    % {
    
    selectedConf_io = unique([(cell2mat(localConfs.ind))',...
                              (cell2mat(localConfs.orn))'],'rows')
    selectedConf_ind = selectedConf_io(:,1);
    selectedConf_orn = selectedConf_io(:,2);
    

    fprintf(1,'\n::::: publish local solutions \n')

    % publish parameters
    % ---------------------------------------------------
    para_.PublishFoV                    = 0;  % draw field of view.
    para_.PublishSymbolicFoV            = 0;
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
    para_.PublishWeakConf               = 0;
    para_.PublishCondConf               = 0;
    para_.PublishHotspotsInOutTSPAngles = 0;
    para_.PublishReferenceDistPoint     = 0;
    para_.PublishGroundTruth            = 0;
    para_.PublishSensingCoverage        = 0;
    para_.ReconstructionFile            = 'coarse_gdm.mat';
    para_.PublishFigureTitle            = 'Local solutions';
    para_.PublishExploration            = 0;
    para_.PublishAdaptivePlannedConf    = 0;
    % directories
    % -------------------------------------------------------------
    %dir_.logs  = dir_.TomographyLogs;
    dir_.recon = dir_.gdplan_recon;

    % data to publish
    % ---------------------------------------------------
    dataPublish.map_candConf        = [];
    dataPublish.FoV                 = FoV;
    dataPublish.FoVfaces            = FoVfaces;
    dataPublish.Conf                = [];
    dataPublish.SensorRange         = SensorRange;
    dataPublish.hotspots            = hotspotSEQ;
    %dataPublish.refDistPoint        = refDistPoint;
    dataPublish.ConfNum             = selectedConf_ind(:);
    dataPublish.ConfOrn             = selectedConf_orn(:);
    dataPublish.FusedConfNum        = [];
    dataPublish.FusedConfOrn        = [];
    dataPublish.RedundantConfNum    = [];
    dataPublish.RedundantConfOrn    = [];
    dataPublish.map_env             = map_env;

    % publish
    % ---------------------------------------------------
    fPublish_GDMPlan( dataPublish,para_,dir_ )

    %pause(3)
    %file_name = sprintf('spp_local_solutions');
    %print('-painters','-dpdf','-r500',[dir_.TomographyFusion,file_name]); % pdf
    %pause(2)

    pause
    %}
end    

%***********************************************************************
%               First conf pair for fusion
%***********************************************************************
fprintf(1,'-- initialize iterative fusion \n')

% selectedFusedConf_ind = localConfs_ind;
% selectedFusedConf_orn = localConfs_orn;

currentFusedConfs = localConfs;
PairNumberToSelect = 1;
[ ConfPair_ioh ] = fIterativeFusionSetup2t( map_env,...
                                            PairNumberToSelect,...
                                            currentFusedConfs,...
                                            cellsize_env,...
                                            para_ );
                                            
%***********************************************************************
%               Pair-wise fusion
%***********************************************************************
iteration_num = 1;
while ...
        size(ConfPair_ioh,1)>1 ...
        && ...
        iteration_num <= para_.maxFusionIterations_2t
    
    fprintf(1,['\n***********************************\n'...
               ':: iteration %02d \n'...
               '***********************************\n'],iteration_num)
    
           
    fprintf(1,'-- preprocess pairwise fusion \n')
    
    %///////////////////////////////////
    %       PREPROCESSING
    %///////////////////////////////////
    visualize = 0;
    % {
    [ wV,...
      D,...
      oV,...
      Conf_DistOrnGainV,...
      fusionHots1,...
      fusionHots2,...
      confKept,...
      map_candConf,...
      fixedConfs_num,...
      fixedConfsAll_num,...
      fixedConfs1_ind,...
      fixedConfs2_ind,...
      NumOfFixeConf,...
      redundantConfs_ind,...
      redundantConfs_orn,...
      oConfOrns,...
      hcs_fusion_ind,...
      hcs_fusion_num,...
      FusionNotPossibleDueToInsufficientCandConf,...
      conf_crossAngles_G ] = fPairwiseFusion2t( ConfPair_ioh,...
                                                currentFusedConfs,...
                                                localNumOfSelectedConfs,...
                                                refPointDist,...
                                                hotspotSEQ,...
                                                map_env,...
                                                map_recon,...
                                                FoV,...
                                                map_coverage,...
                                                SensorRange,...
                                                OBSenv,...
                                                OBSconf,...
                                                cellsize_env,...
                                                visualize,...
                                                para_ );
                                                
    save([dir_.gdmplan_spp,'pairwise_fusion_preprocess.mat'],...
        'wV',...
        'D',...
        'oV',...
        'Conf_DistOrnGainV',...
        'fusionHots1',...
        'fusionHots2',...
        'confKept',...
        'map_candConf',...
        'fixedConfs_num',...
        'fixedConfsAll_num',...
        'NumOfFixeConf',...
        'oConfOrns',...
        'hcs_fusion_ind',...
        'hcs_fusion_num',...
        'conf_crossAngles_G');
    
    if FusionNotPossibleDueToInsufficientCandConf == 1
        break
    end
     
    %///////////////////////////////////
    %       FUSION OPTIMIZATION
    %///////////////////////////////////
    
    cond_NumOfFixeConf = 1;
    for i = 1:size(wV,2)
        cond_NumOfFixeConf = cond_NumOfFixeConf*size(wV{i},2)>NumOfFixeConf;
    end
    
    if cond_NumOfFixeConf %size(wV1,2)>NumOfFixeConf && size(wV2,2)>NumOfFixeConf

        %fprintf(1,'\n:: **** optimization for pairwise fusion **** \n\n')
        
        
        % -- objective values / ERQ
        G  = cell(1,numel(hcs_fusion_ind));
        Uc = cell(1,numel(hcs_fusion_ind));
        Ud = cell(1,numel(hcs_fusion_ind));
        U  = cell(1,numel(hcs_fusion_ind));
        
        for i = 1:numel(hcs_fusion_ind)
            nc    = localNumOfSelectedConfs(hcs_fusion_num(i));
            G{i}  = (para_.alpha*para_.beta)*(1/(nc*(nc-1)))*conf_crossAngles_G{i};
            Uc{i} = (para_.alpha*(1-para_.beta))/nc * ...
                    ((para_.gamma*sum(wV{i},1))+((1-para_.gamma)*(Conf_DistOrnGainV{i})') );
            Ud{i} = (1-para_.alpha)/nc*(1-(D{i}/max(D{i})));
            U{i}  = (Uc{i}'+Ud{i});
            
        end
        
        NumOfAllowedConfFusion = NumOfFixeConf+1;
        strategy = para_.MissionStrategy;
        EXP_TITLE = para_.ExperimentTitle;
        sensing_system = para_.SensingSystem;
        alpha = para_.alpha;
        beta = para_.beta;
        gamma = para_.gamma;
        FilePath = dir_.gdmplan_spp;
        save('DATA_SPP_FUSION_2t.mat',...
             'EXP_TITLE',...
             'strategy',...
             'sensing_system',...
             'NumOfAllowedConfFusion',...             
             'alpha',...
             'beta',...
             'gamma',...
             'G',...
             'U',...
             'FilePath');
        
        % -- solve optimization problem
        PythonTime_start = tic;
        
        system('python optimization/mOptimizationFusedGeometry2t.py');
               
        tPythonTime = 1e-4*(round(toc(PythonTime_start)*1e4));
        fprintf(1,'--- Computation time (from Matlab): %0.4f sec \n',tPythonTime)

        
        fprintf(1,'\n results \n')

        % load SPP solution
        %selectedConf = load([dir_.FusionProcess,...
        %    'spp_selectedconf_pairwise_fusion',num2str(PairNum),'.dat']);
        filename = 'selectedConf_PairwiseFusion.dat';
        selectedConf = load([dir_.gdmplan_spp,filename]);
        selectedConf = round(selectedConf);

        C = selectedConf;
        
        fusionERGs = cell(1,numel(hcs_fusion_ind));
        fusionERQs = cell(1,numel(hcs_fusion_ind));
        for i = 1:numel(hcs_fusion_ind)
            %fusionERGs{i} = C'*(G{i}*C+U{i});
            thisERG = C'*(G{i}*C+U{i});
            
            %-- current ERQ (map quality)
            %---------------------------------------------
            nc = localNumOfSelectedConfs(hcs_fusion_num(i));
            if     nc == 2
                thisERQ = para_.MapQuality2Conf*thisERG;
            elseif nc == 3
                thisERQ = para_.MapQuality3Conf*thisERG;
            elseif nc == 4
                thisERQ = para_.MapQuality4Conf*thisERG;
            elseif nc == 5
                thisERQ = para_.MapQuality5Conf*thisERG;
            end
            
            fprintf('---- ERQ fused #%d: %f\n',i,thisERQ);
            
            %-- update
            fusionERGs{i} = thisERG;
            fusionERQs{i} = thisERQ;
        end
        
        % -- ERQ/ERG difference
        localERGs = cell(1,numel(hcs_fusion_ind));
        localERQs = cell(1,numel(hcs_fusion_ind));
        for i = 1:numel(hcs_fusion_ind)
            %ERQ___local(hcs_fusion_num(i))
            localERGs{i} = localERGall(hcs_fusion_num(i));
            localERQs{i} = localERQall(hcs_fusion_num(i));
            
            fprintf('---- ERQ local #%d: %f\n',i,localERQs{i});
            
        end
        
                      
        ERQ_acceptable = cell(1,numel(hcs_fusion_ind));
        for i = 1:numel(hcs_fusion_ind)
            
            %ERQ_acceptable{i} = (1-para_.toleranceERQFusion2SE)*localERQs{i};
            ERQ_acceptable{i} = ...
                min( [(1-para_.toleranceERQFusion2t)*localERQs{i},...
                       para_.desiredERQFusion2t] );
            
            fprintf('---- ERQ acceptable #%d: %f\n',i,ERQ_acceptable{i});
        end
        
        %fprintf('\n-----\n');
        

        % -- this solution
        selectedConf_ind_this = confKept(selectedConf==1);
        selectedConf_orn_this = oConfOrns(confKept(selectedConf==1));
                
        fixedConfs_ind = unique([fixedConfs1_ind,fixedConfs2_ind]);
        fusedConf_ind  = selectedConf_ind_this(ismember(selectedConf_ind_this,fixedConfs_ind)==0);
        fusedConf_orn  = oConfOrns(fusedConf_ind);
        
        
        %visualize = 1;
        
        cond_ERQimprovement = 1;
        for i = 1:numel(hcs_fusion_ind)
                       
            cond_ERQimprovement = ...
                cond_ERQimprovement * ...
                (fusionERQs{i} >= ERQ_acceptable{i});
            
            fprintf('---- ERQ cond #%d: %d\n',...
                i,(fusionERQs{i} >= ERQ_acceptable{i}));
                        
        end
        
        
        %///////////////////////////////////////////////////////////////
        % ADOPT THE SOLUTION IF IT ADDS SIGNIGICANT IMPROVEMENT
        %///////////////////////////////////////////////////////////////

        %if all(ERQ_Fusion1 >= ERQ_Local1-ERQ_delta1) && all(ERQ_Fusion2 >= ERQ_Local2-ERQ_delta2)
        if cond_ERQimprovement

            %**************************************************************
            % -- update selected fused conf matrix
            %**************************************************************
            
            for i = 1:numel(hcs_fusion_ind)
                ind_num = ismember( currentFusedConfs.ind{hcs_fusion_num(i)},...
                                    redundantConfs_ind ) == 1;
                currentFusedConfs.ind{hcs_fusion_num(i)}(ind_num) = fusedConf_ind;
                currentFusedConfs.orn{hcs_fusion_num(i)}(ind_num) = fusedConf_orn;
            end
            
            
            fprintf(1,'      -- THE CONF PAIR IS FUSED -- \n');

            
            %**************************************************************
            % -- PUBLISH
            %**************************************************************
            if visualize == 1

            % {
            fprintf(1,'\n:: publish this fusion \n')

            % -- orientations
            %ConfOrns = oConfOrns(confKept);            
            %selectedConf_Orn_this = ConfOrns(selectedConf==1);

            fig_title = sprintf('Accepted fusion (Pair#%02d)',iteration_num);
            
            % publish parameters
            % ---------------------------------------------------
            para_.PublishFoV                    = 0;  % draw field of view.
            para_.PublishSymbolicFoV            = 0;
            para_.PublishCandConf               = 1;
            para_.PublishGDM                    = 0;
            para_.PublishHotspots               = 1;
            para_.PublishSubH                   = 0;
            para_.PublishStartPosition          = 0;
            para_.PublishConfPosition           = 1;
            para_.PublishArrow                  = 1;
            para_.PublishConfNum                = 0;
            para_.PublishFusedConfArrow         = 1;
            para_.PublishRedundantConfArrow     = 1;
            para_.PublishWeakConf               = 0;
            para_.PublishCondConf               = 0;
            para_.PublishHotspotsInOutTSPAngles = 0;
            para_.PublishReferenceDistPoint     = 0;
            para_.PublishGroundTruth            = 0;
            para_.PublishSensingCoverage        = 0;
            para_.PublishAdaptivePlannedConf    = 0;
            para_.ReconstructionFile            = 'coarse_gdm.mat';
            para_.PublishFigureTitle            = fig_title;
            para_.PublishExploration            = 0;
            para_.PublishAdaptivePlannedConf    = 0;
            
            % directories
            % -------------------------------------------------------------
            %dir_.logs  = dir_.TomographyLogs;
            dir_.recon = dir_.gdplan_recon;

            % data to publish
            % ---------------------------------------------------
            dataPublish.map_candConf        = map_candConf;
            dataPublish.FoV                 = FoV;
            dataPublish.FoVfaces            = FoVfaces;
            %dataPublish.Conf                = Conf;
            dataPublish.SensorRange         = SensorRange;
            dataPublish.hotspots            = hotspotSEQ([fusionHots1,fusionHots2]);
            %dataPublish.refDistPoint        = refDistPoint;
            dataPublish.ConfNum             = selectedConf_ind_this;
            dataPublish.ConfOrn             = selectedConf_orn_this;
            dataPublish.FusedConfNum        = fusedConf_ind;
            dataPublish.FusedConfOrn        = fusedConf_orn;
            dataPublish.RedundantConfNum    = redundantConfs_ind;
            dataPublish.RedundantConfOrn    = redundantConfs_orn;
            dataPublish.map_env             = map_env;
            
            % publish
            % ---------------------------------------------------
            fPublish_GDMPlan( dataPublish,para_,dir_ )

            %pause
            %pause(3)
            %file_name = sprintf('spp_fusion%02d',PairNum);
            %print('-painters','-dpdf','-r500',[dir_.TomographyFusion,file_name]);
            %}
            end

            % -- update pair number to select
            PairNumberToSelect = 1;
            
        else

            fprintf(2,'      -- THE CONF PAIR IS NOT FUSED -- \n')


            if visualize == 1

            % {
            fprintf(1,'\n:: publish this not accepted fusion \n')

            % -- orientations
            %ConfOrns = oConfOrns(confKept);
            %selectedConf_Orn_this = ConfOrns(selectedConf==1);

            fig_title = sprintf('Not accepted fusion (Pair#%02d)',iteration_num);
            
            % publish parameters
            % ---------------------------------------------------
            para_.PublishFoV                    = 0;  % draw field of view.
            para_.PublishSymbolicFoV            = 0;
            para_.PublishCandConf               = 1;
            para_.PublishGDM                    = 0;
            para_.PublishHotspots               = 1;
            para_.PublishSubH                   = 0;
            para_.PublishStartPosition          = 0;
            para_.PublishConfPosition           = 1;
            para_.PublishArrow                  = 1;
            para_.PublishConfNum                = 0;
            para_.PublishFusedConfArrow         = 1;
            para_.PublishRedundantConfArrow     = 1;
            para_.PublishWeakConf               = 0;
            para_.PublishCondConf               = 0;
            para_.PublishHotspotsInOutTSPAngles = 0;
            para_.PublishReferenceDistPoint     = 0;
            para_.PublishGroundTruth            = 0;
            para_.PublishSensingCoverage        = 0;
            para_.ReconstructionFile            = 'coarse_gdm.mat';
            para_.PublishFigureTitle            = fig_title;
            para_.PublishExploration            = 0;
            para_.PublishAdaptivePlannedConf    = 0;
            
            % directories
            % -------------------------------------------------------------
            %dir_.logs  = dir_.TomographyLogs;
            dir_.recon = dir_.gdplan_recon;

            % data to publish
            % ---------------------------------------------------
            dataPublish.map_candConf        = map_candConf;
            dataPublish.FoV                 = FoV;
            dataPublish.FoVfaces            = FoVfaces;
            %dataPublish.Conf                = Conf;
            dataPublish.SensorRange         = SensorRange;
            dataPublish.hotspots            = hotspotSEQ([fusionHots1,fusionHots2]);
            dataPublish.ConfNum             = selectedConf_ind_this;
            dataPublish.ConfOrn             = selectedConf_orn_this;
            dataPublish.FusedConfNum        = fusedConf_ind;
            dataPublish.FusedConfOrn        = fusedConf_orn;
            dataPublish.RedundantConfNum    = redundantConfs_ind;
            dataPublish.RedundantConfOrn    = redundantConfs_orn;
            dataPublish.map_env             = map_env;

            % publish
            % ---------------------------------------------------
            fPublish_GDMPlan( dataPublish,para_,dir_ )

            %pause
            %pause(3)
            %file_name = sprintf('spp_fusion%02d',PairNum);
            %print('-painters','-dpdf','-r500',[dir_.TomographyFusion,file_name]);
            %}
            end

            % -- update pair number to select
            PairNumberToSelect = PairNumberToSelect+1;
        end

    else
        % -- update pair number to select
        PairNumberToSelect = PairNumberToSelect+1;
    end
    
    %pause
    iteration_num = iteration_num+1;    
    
    fprintf(1,'\n---- conf pair for the next iteration \n')
    
    % -- new conf pair for fusion    
    [ ConfPair_ioh ] = fIterativeFusionSetup2t( map_env,...
                                                PairNumberToSelect,...
                                                currentFusedConfs,...
                                                cellsize_env,...
                                                para_ );
    
end

%-- final list/solution
%--------------------------------
% selectedFusedConf_ind = unique(cell2mat(currentFusedConfs.ind));
% selectedFusedConf_orn = oConfOrns(selectedFusedConf_ind);

% selectedFusedConf_ind = cell2mat(currentFusedConfs.ind);
% selectedFusedConf_orn = cell2mat(currentFusedConfs.orn);

selectedFusedConf_io = unique([(cell2mat(currentFusedConfs.ind))',...
                               (cell2mat(currentFusedConfs.orn))'],'rows');
selectedFusedConf_ind = selectedFusedConf_io(:,1);
selectedFusedConf_orn = selectedFusedConf_io(:,2);


% selectedFusedConf_orn
%**********************************
% save pair-wise fusion results
%**********************************
save([dir_.gdmplan_spp,'pairwise_fusion_results.mat'],...
        'currentFusedConfs',...
        'localConfs',...
        'localERGall',...
        'localNumOfSelectedConfs',...
        'hotspotSEQ',...
        'refPointDist',...
        'selectedFusedConf_ind',...
        'selectedFusedConf_orn',...
        'selectedFusedConf_io');
   
%==================================
% Publish Final Solution
%==================================
visualize = 1;
if visualize == 1

% {
fprintf(1,'\n:: publish final fused solution \n')

% publish parameters
% ---------------------------------------------------
para_.PublishFoV                    = 0;  % draw field of view.
para_.PublishSymbolicFoV            = 0;
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
para_.PublishWeakConf               = 0;
para_.PublishCondConf               = 0;
para_.PublishHotspotsInOutTSPAngles = 0;
para_.PublishReferenceDistPoint     = 0;
para_.PublishGroundTruth            = 0;
para_.PublishSensingCoverage        = 0;
para_.PublishAdaptivePlannedConf    = 0;
para_.ReconstructionFile            = 'coarse_gdm.mat';
para_.PublishFigureTitle            = 'Pairwise fusion result';
para_.PublishExploration            = 0;

% directories
% -------------------------------------------------------------
% dir_.logs  = dir_.TomographyLogs;
dir_.recon = dir_.gdplan_recon;

% data to publish
% ---------------------------------------------------
% dataPublish.map_candConf        = map_candConf;
dataPublish.FoV                 = FoV;
dataPublish.FoVfaces            = FoVfaces;
%dataPublish.Conf                = Conf;
dataPublish.SensorRange         = SensorRange;
dataPublish.hotspots            = hotspotSEQ;
%dataPublish.refDistPoint        = refDistPoint;
dataPublish.ConfNum             = selectedFusedConf_ind(:);
dataPublish.ConfOrn             = selectedFusedConf_orn(:);
dataPublish.FusedConfNum        = [];
dataPublish.FusedConfOrn        = [];
dataPublish.RedundantConfNum    = [];
dataPublish.RedundantConfOrn    = [];
dataPublish.map_env             = map_env;

% publish
% ---------------------------------------------------
fPublish_GDMPlan( dataPublish,para_,dir_ )

pause(3)

% file_name = sprintf('spp_final_solution');
file_name = sprintf('pairwise_fusion_solution');
print('-painters','-dpdf','-r500',[dir_.gdmplan_spp,file_name]); % pdf
end

%}



fprintf(1,'\n--> confs for TSP \n')


selectedConfFusion_ioh = zeros(numel(selectedFusedConf_ind),2+size(hotspotSEQ,1));

for i = 1:numel(selectedFusedConf_ind)
    selectedConfFusion_ioh(i,1:2) = selectedFusedConf_io(i,1:2);
    for j = 1:numel(hotspotSEQ)
        if any( ismember(currentFusedConfs.ind{j},selectedConfFusion_ioh(i,1)) &...
                ismember(currentFusedConfs.orn{j},selectedConfFusion_ioh(i,2)) )
            
            selectedConfFusion_ioh(i,2+j) = 1;
        end
    end
end

disp('----- confs vs hcs')
% sum(selectedConfFusion_ioh(:,3:end),1)
% sum(selectedConfFusion_ioh(:,3:end),2)

%{
selectedConfFusion_ioh = zeros(numel(selectedFusedConf_ind(:)),2+size(hotspotSEQ,1));
for i = 1:numel(hotspotSEQ)
    % -- update new selected conf for ind ont and hotspots
    selectedConfFusion_ioh(((i-1)*para_.NumberOfConf_SPPH)+(1:3),1) = selectedFusedConf_ind(i,:);
    selectedConfFusion_ioh(((i-1)*para_.NumberOfConf_SPPH)+(1:3),2) = selectedFusedConf_orn(i,:);
    selectedConfFusion_ioh(((i-1)*para_.NumberOfConf_SPPH)+(1:3),2+i) = 1;
    
end

uniqConf = unique(selectedConfFusion_ioh(:,1:2),'rows');
selectedConfFusion_ioh_unq = zeros(size(uniqConf,1),size(selectedConfFusion_ioh,2));

for i = 1:size(uniqConf)
    
    selectedConfFusion_ioh_unq(i,1:2) = uniqConf(i,:);
    
    dual_ind = (logical(selectedConfFusion_ioh(:,1)==uniqConf(i,1)) .*...
        logical(selectedConfFusion_ioh(:,2)==uniqConf(i,2)));
    
    selectedConfFusion_ioh_unq(i,3:end) = ...
        sum(selectedConfFusion_ioh(find(dual_ind),3:end),1);
    
end

selectedConfFusion_ioh = selectedConfFusion_ioh_unq;
%}

tspConf = ((selectedConfFusion_ioh(:,1)-1)*FoVfaces.num)+ ...
    (((selectedConfFusion_ioh(:,2)>=045 & selectedConfFusion_ioh(:,2)<=135))*1) +...
    (((selectedConfFusion_ioh(:,2)>=136 & selectedConfFusion_ioh(:,2)<=225))*2) +...
    (((selectedConfFusion_ioh(:,2)>=226 & selectedConfFusion_ioh(:,2)<=315))*3) +...
    (((selectedConfFusion_ioh(:,2)>=316 | selectedConfFusion_ioh(:,2)<=044))*4);


%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% Check for the duplicate confs
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
if numel(tspConf)~=numel(unique(tspConf))
    
    %find(hist(tspConf,unique(tspConf))>1)    
    %A = [7 18 27 42 65 49 54 65 78 82 87 98];
    [n, bin] = histc(tspConf, unique(tspConf));
    multiple = find(n > 1);
    index    = find(ismember(bin, multiple));
    
    %selectedConfFusion_ioh(:,2),045 & selectedConfFusion_ioh(:,2)<=135))*1)
        
    unqConf = unique(tspConf(index));
    
    
    for duplicate_num = 1:numel(unqConf)
        
        index_this = index(tspConf(index) == unqConf(duplicate_num));
        
        offset = ...
            [rad2deg(angleAbsDiff(deg2rad(selectedConfFusion_ioh(index_this,2)),deg2rad(090))),...
             rad2deg(angleAbsDiff(deg2rad(selectedConfFusion_ioh(index_this,2)),deg2rad(180))),...
             rad2deg(angleAbsDiff(deg2rad(selectedConfFusion_ioh(index_this,2)),deg2rad(270))),...
             rad2deg(angleAbsDiff(deg2rad(selectedConfFusion_ioh(index_this,2)),deg2rad(000)))];
         
         [val,ind] = sort(min(offset,[],2));
         
         for i = 1:numel(index_this)
             
             [val,num] = min(offset(ind(i),:));
             tspConf(index_this(ind(i))) = ...
                 ((selectedConfFusion_ioh(index_this(ind(i)),1)-1)*FoVfaces.num)+num;
             offset(:,num) = inf;
         end
    end
end


%**********************************
% save pair-wise fusion results
%**********************************
save([dir_.gdmplan_spp,'pairwise_fusion_results.mat'],...
        'currentFusedConfs',...
        'localConfs',...
        'localERGall',...
        'localNumOfSelectedConfs',...
        'hotspotSEQ',...
        'refPointDist',...
        'selectedFusedConf_ind',...
        'selectedFusedConf_orn',...
        'selectedFusedConf_io',...
        'selectedConfFusion_ioh',...
        'tspConf');
    
    

%**********************************************************************
% Total computation time.
%**********************************************************************
tFusion = 1e-4*(round(toc(tFusioni)*1e4));

fprintf(1,'\nComputation time: %0.4f sec \n',tFusion)


end


function fPublishSetup( dataPublish,para_,dir_ )
% 
% 


zoom_in_environment = 750;
cell_marker_size = 4.85;
color_scale = 'linear'; %
%beam_color = [255,204,153]/256; % [204,229,255]/256
marker_size_cand_conf = 2.5;
marker_size_hotspot = 5; %7;
%marker_size_weak_conf = 4;
%marker_size_cond_conf = 5;
marker_size_conf_position = 2; %5; % icra-submission size: 3
line_width_fov = 1; %2.5
symbolic_fov_length = 2.0; %3;
line_width_symbolic_fov = 0.5; %1.25;

arrow_length = 3;
line_width_arrow = 0.5; %0.25
arrow_head_size = 20;
marker_size_arrow = 1; %2

font_size_conf_num = 6;%16
font_size_conc_scale = 10;
%conf_num_dist = 0.5; % 0.3
%conf_num_dist_adaptive_planned = 0.75;




color_conf                  = [0.0 0.0 0.6];
color_conf_redundatn        = [135,206,235]/255; %[0.7,0.7,0.7];
color_conf_adaptive_planned = [0.4,0.4,0.4];

distance_conf_seq                  = 0.5; % 0.3
distance_conf_seq_adaptive_planned = 0.6;


% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             ENVIRONMENT MAP:
% 
% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% map_env = dlmread([dir_.env,para_.EnvironmentFigFile]);
% map_env = 1-map_env;
map_env = dataPublish.map_env;
[map_row,map_col] = size(map_env);

h = figure('name',para_.PublishFigureTitle); 
imshow(map_env','InitialMagnification',zoom_in_environment);
hold on;
set(gca,'YDir','normal');
colormap('gray'); 

%caxis([-5 2]) % for resulting gdm
%caxis([-12 4]) % for sensor placement solution
% caxis([-0.25 1]) % for sensor placement solution
% caxis([-2 2])
caxis([-1 1.25])
%caxis([-10 1.5]) % for gdm insets


% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             GAS CONCENTRATION MAP
% 
% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
    
    % -- reconstruction 
    load([dir_.recon,para_.ReconstructionFile],...
        'mean_map',...
        'cell_coords_x',...
        'cell_coords_y');
    map_recon = mean_map';
    
    
    % ------------- select color map for gdm --------------
    %gdm_colormap = flipud(colormap('hot')); % for main fig. 
    gdm_colormap = flipud(hot(512)); % for main fig. 
    %gdm_colormap = flipud(colormap('summer')); % for insets
    % -----------------------------------------------------
        
    % ------ color step size --------------------------------
    %max_concentration = 10000; % inset gdm
    %gt = dlmread([dir_.GroundTruthLogs,'map_con.dat']);
    %max_concentration = max(gt(:))
    max_concentration = max(map_recon(:)); % main fig.
    delta = max_concentration/(size(gdm_colormap,1)-1);    
    %delta = max(map_gd(:))/(size(gdm_colormap,1)-1);
    %delta = 500/(size(gdm_colormap,1)-1); % lets limit to 500ppm
    
    % ------------------- nonlinear color code ---------------------  
    if strcmp(color_scale,'nonlinear')
        %b = (exp(256/1))/255;
        %a = 256 /(log(b*256));

        x = 1:1:256;
        %epsilon = 1/(exp(1)-1); % standard
        %epsilon = 1/(exp(0.05)-1); % slower
        epsilon = 1/(exp(0.10)-1); % slower
        %epsilon = 1/(exp(0.15)-1); % slower
        %epsilon = 1/(exp(0.20)-1); % slower
        %epsilon = 1/(exp(0.90)-1); % slower
        %epsilon = 1/(exp(1.5)-1); % standard
        normlogsum_penalty = log(1+(abs(x)/epsilon));
        %normlogsum_penalty = 1.05.^(1+(abs(x)/epsilon));
        normlogsum_penalty = normlogsum_penalty/(normlogsum_penalty(x==max(x))); % normalized
        normlogsum_penalty = round(256*normlogsum_penalty);
        normlogsum_penalty(1)
        normlogsum_penalty(end)
        %figure; plot(x,normlogsum_penalty); % to test

        % publish nonlinear steps
        linear_percentage    = 10:10:100;
        linear_percentage_ind  = round((linear_percentage/100)*numel(x));
        nonlinear_percentage_ind = normlogsum_penalty(linear_percentage_ind);
        nonlinear_percentage = (nonlinear_percentage_ind/numel(x))*100;
        nonlinear_percentage = round(nonlinear_percentage);
        %disp('nonlinear_percentage = ',num2str(nonlinear_percentage))
        disp_stng = ['nonlinear color code scale is: ',num2str(nonlinear_percentage)];
        %disp(disp_stng)
        disp(disp_stng)
    end
    % ----------------------------------------------------------------
    
    
    for i = 1:numel(cell_coords_x)%1:numel(cell_coords_y)%map_row
        for j = 1:numel(cell_coords_y)%1:numel(cell_coords_x)%map_col
            % if positive concentration
            
            if map_recon(i,j)>=0
                
                if map_env(i,j)>0
                % ---------------------- useful trash -------------------
                %col = gdm_colormap(round(map_gd(i,j)/delta)+1,:);
                
                % --- if concentration is greater than 1000 ppm, its still 1000 ppm
                %col = gdm_colormap( round(min(map_gd(i,j),3000)/delta)+1, :);
                %col = gdm_colormap( round(min(map_gd(i,j),1.78e+04)/delta)+1, :);
                %col = gdm_colormap( round(min(map_gd(i,j),inf)/delta)+1, :);
                % -------------------------------------------------------
                
                % ----------------------------------------------------
                % linear color code
                % ----------------------------------------------------
                if strcmp(color_scale,'linear')
                % {
                linear_color_number = round( min(map_recon(i,j),max_concentration)/delta)+1;
                col = gdm_colormap( linear_color_number,:);                
                %}
                end
                % ----------------------------------------------------
                
                % ----------------------------------------------------
                % nonlinear color code
                % ----------------------------------------------------
                if strcmp(color_scale,'nonlinear')
                % {
                linear_color_number = round( min(map_recon(i,j),max_concentration)/delta)+1;
                %nonlinear_color_number = round( a*log(b*linear_color_number) );
                %nonlinear_color_number = linear_color_number;
                nonlinear_color_number = normlogsum_penalty(linear_color_number);
                %if nonlinear_color_number <1
                %    nonlinear_color_number = 1;
                %end
                col = gdm_colormap( nonlinear_color_number,:);                
                %}
                end
                % ----------------------------------------------------
                
                %----- to avoid black borders of earch cell -----
                % {
                if sum(col) == 3
                    col = col-(1e-10);
                end
                %}
                % ----- plot ------------
                plot(cell_coords_x(i)+0.0,cell_coords_y(j)+0.0,...
                    's',...
                    'MarkerEdgeColor',col,...
                    'MarkerFaceColor',col,...
                    'MarkerSize',cell_marker_size);
                end
            end
        end
    end
    
    % --- concentration scale.
    conc_text = sprintf('Concentration scale: %d -- %.2d ppm',0,max_concentration);
    %{
    text(map_row-2,2,conc_text,...
            'HorizontalAlignment','right',...
            'VerticalAlignment','bottom',...
            'FontWeight','demi',...
            'FontSize',font_size_conc_scale,...
            'Color',[0,0,0],...
            'BackgroundColor',[1,1,1],...
            'Interpreter','latex',...
            'EdgeColor',[0,0,0],...
            'Margin',1,...
            'FontWeight','bold');
    %}

%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             CANDIDATE CONFIGURATIONS MAP
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishCandConf
    
    map_candConf = dataPublish.map_candConf;
    
    candConf_ind = find(map_candConf);
    [candConf_r,candConf_c] = ind2sub(size(map_env),candConf_ind);
    plot(candConf_r,candConf_c,...
        'x',...
        'Color',[102,178,255]/255,...
        'MarkerSize',marker_size_cand_conf);
end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             HOTSPOTS:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
   
hotspots = dataPublish.hotspots;

[hot_r,hot_c] = ind2sub(size(map_env),hotspots);
plot(hot_r,hot_c,...
    'ok',...
    'MarkerSize',marker_size_hotspot,...
    'MarkerFaceColor',[1,0,0]);




%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             MARK SUBJECT HOTSPOTS:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishSubH
    
    subHotspots = dataPublish.subHotspots;    
    plot(subHotspots.sub(:,1),subHotspots.sub(:,2),...
        'og','MarkerSize',marker_size_hotspot);
end



%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                        MARK REDUNDANT CONFIGURATIONS:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
    
RedundantConf_ind = dataPublish.RedundantConfNum;
RedundantConf_orn = dataPublish.RedundantConfOrn;

redConfCell_num = RedundantConf_ind;
[redConfCell_r,redConfCell_c] = ind2sub(size(map_env),redConfCell_num);

% -- position
plot(redConfCell_r,redConfCell_c,...
    'o',...
    'Color',color_conf_redundatn,...
    'MarkerFaceColor',color_conf_redundatn,...
    'MarkerSize',marker_size_conf_position+1);

% -- arrow
for i = 1:numel(RedundantConf_ind)

    quiver(redConfCell_r(i),redConfCell_c(i),...
           arrow_length*cosd(RedundantConf_orn(i)),...
           arrow_length*sind(RedundantConf_orn(i)),...
           'LineWidth',line_width_arrow+1,...
           'color',color_conf_redundatn,...
           'MaxHeadSize',arrow_head_size,...
           'MarkerSize',marker_size_arrow);
end


FoV = dataPublish.FoV;

%arrow_length = 2;
%line_length  = 2;    

% -- symbolic field of view
for i = 1:numel(RedundantConf_ind)

    % -- drat sumbolic field of view
    start_angle  = RedundantConf_orn(i)-FoV/2; %FoVfaces.ang(nCONF(i),1);
    end___angle  = RedundantConf_orn(i)+FoV/2; %start_angle+FoV;
    r            = symbolic_fov_length;
    x1           = r*cosd(start_angle);
    y1           = r*sind(start_angle);
    x2           = r*cosd(end___angle);
    y2           = r*sind(end___angle);
    plot([redConfCell_r(i),x1+redConfCell_r(i)],...
         [redConfCell_c(i),y1+redConfCell_c(i)],...
          'color',color_conf_redundatn,...
          'LineWidth',line_width_symbolic_fov+1); 
    plot([redConfCell_r(i),x2+redConfCell_r(i)],...
         [redConfCell_c(i),y2+redConfCell_c(i)],...
          'color',color_conf_redundatn,...
          'LineWidth',line_width_symbolic_fov+1);
end



%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             MARK CONFIGURATIONS:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishConfPosition
    
    ConfNum = dataPublish.ConfNum;
    ConfCell_num = ConfNum;
    [ConfCell_r,ConfCell_c] = ind2sub(size(map_env),ConfCell_num);
    plot(ConfCell_r,ConfCell_c,...
        'o',...
        'Color',color_conf,...
        'MarkerFaceColor',color_conf,...
        'MarkerSize',marker_size_conf_position);
end




%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             DRAW SYMBOLIC FIELD OF VIEW:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishSymbolicFoV
    
    FoV = dataPublish.FoV;
    
    %arrow_length = 2;
    %line_length  = 2;    
    
    ConfNum = dataPublish.ConfNum;
    ConfCell_num = ConfNum;
    [ConfCell_r,ConfCell_c] = ind2sub(size(map_env),ConfCell_num);
    ConfOrn = dataPublish.ConfOrn;
    
    for i = 1:numel(ConfNum)
        
        % -- drat sumbolic field of view
        start_angle  = ConfOrn(i)-FoV/2; %FoVfaces.ang(nCONF(i),1);
        end___angle  = ConfOrn(i)+FoV/2; %start_angle+FoV;
        r            = symbolic_fov_length;
        x1           = r*cosd(start_angle);
        y1           = r*sind(start_angle);
        x2           = r*cosd(end___angle);
        y2           = r*sind(end___angle);
        plot([ConfCell_r(i),x1+ConfCell_r(i)],...
             [ConfCell_c(i),y1+ConfCell_c(i)],...
              'color',color_conf,'LineWidth',line_width_symbolic_fov); 
        plot([ConfCell_r(i),x2+ConfCell_r(i)],...
             [ConfCell_c(i),y2+ConfCell_c(i)],...
              'color',color_conf,'LineWidth',line_width_symbolic_fov); 
    end
end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             DRAW CONF ARROWS:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishArrow
    
    ConfNumbers      = dataPublish.ConfNum;
    ConfOrn = dataPublish.ConfOrn;
    
    [ConfCell_r,ConfCell_c] = ind2sub(size(map_env),ConfNumbers);
    
    for i = 1:length(ConfNumbers)
        
        conf_to_plot = [ConfNumbers(i),ConfCell_r(i),ConfCell_c(i),ConfOrn(i)];
        
        %FoVfaces.lgt
        %nCONF(i)
        %FoVfaces.lgt(nCONF(i))
        
        % Draw arrow
        quiver(ConfCell_r(i),ConfCell_c(i),...
               arrow_length*cosd(ConfOrn(i)),...
               arrow_length*sind(ConfOrn(i)),...
               'LineWidth',line_width_arrow,...
               'color',color_conf,...
               'MaxHeadSize',arrow_head_size,... %'Marker','o',...
               'MarkerSize',marker_size_arrow);
        %pause       
        
        %text(ConfCell_r(i),ConfCell_c(i),'\rightarrow','FontSize',15,'Rotation',ConfOrn(i));        
        %pause
        
        %arrow([ConfCell_r(i),ConfCell_c(i)],...
        %    [ConfCell_r(i)+arrow_length*cosd(ConfOrn(i)),ConfCell_c(i)+arrow_length*sind(ConfOrn(i))])
        %pause
        
    end
end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             REDUNDANT CONF ARROWS:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishRedundantConfArrow
    
    
    redundantConfNum = dataPublish.RedundantConfNum;
    redundantConfOnt = dataPublish.RedundantConfOrn;
    
    [redundantConfCell_r,redundantConfCell_c] = ind2sub(size(map_env),redundantConfNum);
    
    for i = 1:length(redundantConfNum)
        
        % Draw arrow
        quiver(redundantConfCell_r(i),redundantConfCell_c(i),...
               2*cosd(redundantConfOnt(i)),...
               2*sind(redundantConfOnt(i)),...
               'LineWidth',line_width_arrow,...
               'color','m',...
               'MaxHeadSize',arrow_head_size,...
               'Marker','o',...
               'MarkerSize',marker_size_arrow);
        %pause
    end
end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             FUSED CONF ARROWS:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishFusedConfArrow
    
    
    fusedNum = dataPublish.FusedConfNum;
    fusedOnt = dataPublish.FusedConfOrn;
    
    [fusedCell_r,fusedCell_c] = ind2sub(size(map_env),fusedNum);
    
    for i = 1:length(fusedNum)
        
        % Draw arrow
        quiver(fusedCell_r(i),fusedCell_c(i),...
               2.5*cosd(fusedOnt(i)),...
               2.5*sind(fusedOnt(i)),...
               'LineWidth',line_width_arrow,...
               'color','g',...
               'MaxHeadSize',arrow_head_size,...
               'Marker','o',...
               'MarkerSize',marker_size_arrow);
        %pause
    end
end

    

% hold off
% alpha(.1)

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(2), pos(4)])












end


% --------------------------- End of the document ----------------------------------------


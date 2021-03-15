function [ reselectedConf2Exe_ind,...
           reselectedConf2Exe_orn,...
           tFusion ] = fFusionHotspotGeometries1t( FoV,...
                                                   hc_current,...
                                                   hc_nums_local_this,...
                                                   hc_block_num,...
                                                   SensorRange,...
                                                   FoVfaces,...
                                                   map_env,...
                                                   map_gd_current,...
                                                   map_coverage,...
                                                   confs_executed_hc,...
                                                   refDistPoint,...
                                                   OBSenv,...
                                                   OBSconf,...
                                                   cellsize_env,...
                                                   dir_,...
                                                   para_,...
                                                   visualize )
%
% 


tFusioni = tic; % initialize pre computational time.


fprintf('--> hc block num: %02d\n',hc_block_num);

hcs_sub = hc_current;
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
%      load preprocess maps
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% fprintf(1,'\n:: load preprocess \n')
% 
% load([dir_.TomographyMaps,'preprocess_maps.mat'],...
%         'FoVfaces',...
%         'hotspotSEQ',...
%         'map_env',...
%         'map_recon',...
%         'map_coverage',...
%         'OBSenv',...
%         'OBSconf',...
%         'SensorRange');


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
%    initial setup for fusion
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
fprintf(1,'---- initial setup for fusion \n');


% {
[ localConfs,...
  localERGall,...
  localERQall,...
  localNumOfSelectedConfs] = fInitialFusionSetup1t( hcs_sub,...
                                                    hc_nums_local_this,...                                                     
                                                    map_env,...
                                                    dir_,...
                                                    para_,...
                                                    visualize );

filename = sprintf('initial_fusion_setup_block%02d.mat',hc_block_num);
save([dir_.sppgdm,filename],...
        'localConfs',...
        'localERGall',...
        'localERQall',...
        'localNumOfSelectedConfs');
%}
% load([dir_.FusionProcess,'initial_fusion_setup.mat']);


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
    para_.PublishFigureTitle            = 'Local solutions 1t-rms';
    para_.PublishExploration            = 0;
    para_.PublishAdaptivePlannedConf    = 0;

    % directories
    % -------------------------------------------------------------
    %dir_.logs  = dir_.TomographyLogs;
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


    %}
end    

%***********************************************************************
%               First conf pair for fusion
%***********************************************************************
fprintf(1,'---- initialize iterative fusion \n');

% selectedFusedConf_ind = selectedConf_ind;
% selectedFusedConf_orn = selectedConf_ont;
% PairNumberToSelect = 1;
% [ ConfFusionPair_ioh ] = fIterativeFusionSetup( map_env,...
%                                                 PairNumberToSelect,...
%                                                 selectedFusedConf_ind,...
%                                                 selectedFusedConf_orn,...
%                                                 dataYAML,...
%                                                 para_ );

currentFusedConfs = localConfs;
PairNumberToSelect = 1;                                            
[ ConfPair_ioh ] = fIterativeFusionSetup1t( map_env,...
                                            PairNumberToSelect,...
                                            confs_executed_hc,...
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
        iteration_num <= para_.maxFusionIterations_1t %20
    
    fprintf(1,['\n***********************************\n'...
               ':: iteration %02d \n'...
               '***********************************\n'],iteration_num)
           
           
    fprintf(1,'---- preprocess pairwise fusion \n')
    
    
    %///////////////////////////////////
    %       PREPROCESSING
    %///////////////////////////////////
    
    visualize = 0;
    
    [ wV,...
      D,...
      oV,...
      Conf_DistOrnGainV,...
      fusionHcs1,...
      fusionHcs2,...
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
      conf_crossAngles_G ] = fPairwiseFusion1t( ConfPair_ioh,...
                                                currentFusedConfs,...
                                                localNumOfSelectedConfs,...
                                                hcs_sub,...
                                                refDistPoint,...
                                                map_env,...
                                                map_gd_current,...
                                                FoV,...
                                                map_coverage,...
                                                SensorRange,...
                                                OBSenv,...
                                                OBSconf,...
                                                cellsize_env,...
                                                visualize,...
                                                para_ );
    
    filename = sprintf('pairwise_fusion_block%02d_preprocess.mat',hc_block_num);
    save([dir_.sppgdm,filename],...
        'wV',...
        'D',...
        'oV',...
        'Conf_DistOrnGainV',...
        'fusionHcs1',...
        'fusionHcs2',...
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

        %fprintf(1,'\n:: -- optimization for pairwise fusion -- \n\n')
        
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
        FilePath = dir_.sppgdm;
        save('DATA_SPP_FUSION_1t.mat',...
             'EXP_TITLE',...
             'strategy',...
             'sensing_system',...
             'NumOfAllowedConfFusion',...
             'hc_block_num',...
             'alpha',...
             'beta',...
             'gamma',...
             'G',...
             'U',...
             'FilePath');
         
         
        % -- solve optimization problem
        PythonTime_start = tic;
        
        system('python optimization/mOptimizationFusedGeometry1t.py');
        
        tPythonTime = 1e-4*(round(toc(PythonTime_start)*1e4));
        fprintf(1,'--- Computation time (from Matlab): %0.4f sec \n',tPythonTime)
        

        %fprintf(1,'\n:: results \n')
        fprintf(1,'\n results: \n')

        % load SPP solution
        filename = sprintf('pairwise_fusion_block%02d_spp.dat',hc_block_num);
        selectedConf = load([dir_.sppgdm,filename]);
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
            
            %ERQ_acceptable{i} = (1-para_.toleranceERQFusion1SE)*localERQs{i};
            ERQ_acceptable{i} = ...
                min( [(1-para_.toleranceERQFusion1t)*localERQs{i},...
                       para_.desiredERQFusion1t] );
            
            fprintf('---- ERQ acceptable #%d: %f\n',i,ERQ_acceptable{i});
        end
        
        
        selectedConf_ind_this = confKept(selectedConf==1);
        selectedConf_orn_this = oConfOrns(confKept(selectedConf==1));
        
        fixedConfs_ind = unique([fixedConfs1_ind,fixedConfs2_ind]);
        fusedConf_ind = selectedConf_ind_this(ismember(selectedConf_ind_this,fixedConfs_ind)==0);
        fusedConf_orn = oConfOrns(fusedConf_ind);
                                
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
        %if cond_ERGimprovement
        if cond_ERQimprovement
            
            %**************************************************************
            % -- update selected fused conf matrix
            %**************************************************************

            for i = 1:numel(hcs_fusion_ind)
                ind_num = ismember(currentFusedConfs.ind{i},redundantConfs_ind)==1;
                currentFusedConfs.ind{i}(ind_num) = fusedConf_ind;
                currentFusedConfs.orn{i}(ind_num) = fusedConf_orn;
            end
            
            
            %fprintf(1,'\n:: ***** CONFIGURATIONS PAIR IS FUSED SUCCESSFULLY ****** \n')
            %fprintf(1,'\n  -- THE CONF PAIR IS FUSED -- \n')
            fprintf(1,'\n      -- THE CONF PAIR IS FUSED -- \n');
            
            %=======================================================================
            % Publish
            %=======================================================================
            
            
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
            dir_.recon = dir_.Reconstruction;

            % data to publish
            % ---------------------------------------------------
            dataPublish.map_candConf        = map_candConf;
            dataPublish.FoV                 = FoV;
            dataPublish.FoVfaces            = FoVfaces;
            %dataPublish.Conf                = Conf;
            dataPublish.SensorRange         = SensorRange;
            dataPublish.hotspots            = hcs_fusion_ind;
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

            %fprintf(2,'\n:: ***** CONFIGURATIONS PAIR IS NOT FUSED ****** \n')
            %fprintf(2,'\n  -- THE CONF PAIR IS NOT FUSED -- \n')
            fprintf(2,'\n      -- THE CONF PAIR IS NOT FUSED -- \n')

            if visualize == 1

            % {
            fprintf(1,'\n:: publish this fusion \n')

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
            dir_.recon = dir_.Reconstruction;

            % data to publish
            % ---------------------------------------------------
            dataPublish.map_candConf        = map_candConf;
            dataPublish.FoV                 = FoV;
            dataPublish.FoVfaces            = FoVfaces;
            %dataPublish.Conf                = Conf;
            dataPublish.SensorRange         = SensorRange;
            dataPublish.hotspots            = hcs_fusion_ind;
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
    
    %fprintf(1,'\n:: conf pair for the next iteration \n')
    fprintf(1,'\n---- conf pair for the next iteration \n')
    
    % -- new conf pair for fusion
    [ ConfPair_ioh ] = fIterativeFusionSetup1t( map_env,...
                                                PairNumberToSelect,...
                                                confs_executed_hc,...
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
filename = sprintf('pairwise_fusion_block%02d_results.mat',hc_block_num);
save([dir_.sppgdm,filename],...
        'currentFusedConfs',...
        'localConfs',...
        'localERGall',...
        'localNumOfSelectedConfs',...
        'refDistPoint',...
        'selectedFusedConf_ind',...
        'selectedFusedConf_orn',...
        'selectedFusedConf_io');

    
    
%*************************************************************************
%  SORT SELECTED CONF WITH TRAVELING DISTANCE FROM THE LAST EXECUTED CONF
%*************************************************************************

[selectedFusedConf_r,selectedFusedConf_c] = ind2sub(size(map_env),selectedFusedConf_ind);

dist2selectedConf = fPathDist2(refDistPoint,[selectedFusedConf_r,selectedFusedConf_c],map_env);
[val,ind] = sort(dist2selectedConf);

reselectedConf2Exe_ind = selectedFusedConf_ind(ind);
reselectedConf2Exe_orn = selectedFusedConf_orn(ind);


%=======================================================================
% Publish Final Solution
%=======================================================================
% visualize = 1;
if visualize == 1

% {
fprintf(1,'\n:: publish final fused solution \n')

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
dir_.recon = dir_.Reconstruction;

% data to publish
% ---------------------------------------------------
dataPublish.map_candConf        = map_candConf;
dataPublish.FoV                 = FoV;
dataPublish.FoVfaces            = FoVfaces;
%dataPublish.Conf                = Conf;
dataPublish.SensorRange         = SensorRange;
dataPublish.hotspots            = hcs_fusion_ind;
%dataPublish.refDistPoint        = refDistPoint;
dataPublish.ConfNum             = reselectedConf2Exe_ind;
dataPublish.ConfOrn             = reselectedConf2Exe_orn;
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
print('-painters','-dpdf','-r500',[dir_.sppgdm,file_name]); % pdf
end



%******************************************************* 
%-- remove already executed confs for this hc
%******************************************************* 
% confs_executed_hc
executedConfs_ind = sub2ind(size(map_env),...
                            round(confs_executed_hc(:,1)-0.5),...
                            round(confs_executed_hc(:,2)-0.5));
executedConfs_orn = confs_executed_hc(:,3);

%-- adjust 0 and 360 dg issue
reselectedConf2Exe_orn(reselectedConf2Exe_orn==360) = 0;
executedConfs_orn(executedConfs_orn==360) = 0;


ind_to_remove = (ismember(reselectedConf2Exe_ind,executedConfs_ind) &...
    ismember(reselectedConf2Exe_orn,executedConfs_orn));


if size(find(ind_to_remove),2) < 1
    warning("The executed confs are not removed. Check the angles [0,360].");
end

reselectedConf2Exe_ind(ind_to_remove) = [];
reselectedConf2Exe_orn(ind_to_remove) = [];


disp(['--- Selected confs ind: ',num2str(reselectedConf2Exe_ind')]);




%**********************************************************************
% Total computation time.
%**********************************************************************
tFusion = 1e-4*(round(toc(tFusioni)*1e4));

% fprintf(1,'\nComputation time: %0.4f sec \n',tFusion)
fprintf(1,'\n---- Computation time: %0.4f sec \n',tFusion);


end

% --------------------------- End of the document ----------------------------------------


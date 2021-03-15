function [ redC,...
           reducedCoverage ] = fPlan_GD_ReducedCoverage( fullC,...
                                                          V,...
                                                          oV,...
                                                          FoV,...
                                                          FoVfaces,...
                                                          confKept,...
                                                          SensorRange,...
                                                          map_env,...
                                                          para_,...
                                                          dir_,...
                                                          visualize )

%
%
%


if visualize 
    %figure;
    para_.PublishFigureTitle = 'BEFORE reducing confs based on coverage ind.';
    fPublishThisConfs( fullC,oV,FoV,FoVfaces,confKept,SensorRange,map_env,para_,dir_ )
    pause
end



Vconf = V(:,fullC);

% global_fci = size(V,1);

Vconf = full(Vconf);

Confs_fci = sum(Vconf,1);

[fci_val,fci_ind] = sort(Confs_fci,'descend');

% currentMapCoverage = 1.00;

% C_fullCoverage = C;
% C_reducedCoverage = C;

lastMapCoverage = 1;

ConfRemoved = [];

% condKeepShrinking = 1;
find(fullC)
if numel(find(fullC))>1
    condKeepShrinking = 1
else
    condKeepShrinking = 0
end

num = 1;

while condKeepShrinking
    
    %num 
    %fci_ind(end)
    
    % Vconf except the conf subject to remove
    otherConfs_Vconf = Vconf;
    otherConfs_Vconf(:,fci_ind(end)) = [];
    
    thisConf_Vconf = Vconf(:,fci_ind(end));
    
    otherConfs_fcc = find(sum(otherConfs_Vconf,2));
    thisConf_fcc = find(sum(thisConf_Vconf,2));
    
    thisConf_ecc = thisConf_fcc;
    thisConf_ecc(ismember(thisConf_fcc,otherConfs_fcc)) = [];
    
    
    %thisConf_pci (potential coverage index)
    thisConf_fci = numel(thisConf_fcc);
    thisConf_eci = numel(thisConf_ecc);
    otherConfs_fci = numel(otherConfs_fcc);
    
    
    
    ratio_eci_fci = thisConf_eci/thisConf_fci;    
    currenMapCoverage = otherConfs_fci/size(V,1);
    
    
    %if ratio_eci_fci<=0.15 && currenMapCoverage>0.85
    if ratio_eci_fci <= para_.ExclusiveToTotalCoverageThreshold &&...
            currenMapCoverage > para_.minSensingCoverage1t
        
        condKeepShrinking = 1;
        
        Vconf(:,fci_ind(end)) = 0;
        
        ConfRemoved = [ConfRemoved,fci_ind(end)];
        
        fci_ind(end) = [];
        fci_val(end) = [];
        
        num = num+1;
        
        lastMapCoverage = currenMapCoverage;
        
    else
        condKeepShrinking = 0;
    end
    
    
    
    
end


selectedConfInd = find(fullC);

redC = fullC;
redC(selectedConfInd(ConfRemoved)) = 0;
% new_selectedConfInd = find(Cred);
% Cfull = Cred;
% numel(find(C))


if visualize 
    %figure;
    para_.PublishFigureTitle = 'AFTER reducing confs based on coverage ind.';
    fPublishThisConfs( redC,oV,FoV,FoVfaces,confKept,SensorRange,map_env,para_,dir_ )
    pause
end

reducedCoverage = lastMapCoverage;

fprintf('Num of confs. (full coverage):    %02d\n',numel(find(fullC)));
fprintf('Num of confs. (reduced coverage): %02d\n',numel(find(redC)));
fprintf('Coverage reduced to:              %0.2f%%\n',reducedCoverage*100);




   
end


function fPublishThisConfs( C,oV,FoV,FoVfaces,confKept,SensorRange,map_env,para_,dir_ )
%
%


Conf = zeros(size(oV,2),1);
Conf(confKept) = C;

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
% para_.PublishFigureTitle            = 'An exploration plan for gas detection';

% data to publish
% -------------------------------------------------------------
dataPublish.FoV                 = FoV;%FoVragd;
dataPublish.FoVfaces            = FoVfaces;
dataPublish.Conf                = Conf;
dataPublish.SensorRange         = SensorRange;
% dataPublish.confSequenceNum     = confSequenceNum;
% dataPublish.robotStartCellInMap = origin_conf;
dataPublish.map_env             = map_env;

% publish
% -------------------------------------------------------------
fPublish_GD( dataPublish,para_,dir_ )
% export_fig('test_detection_results_after', '-pdf','-r500')


end

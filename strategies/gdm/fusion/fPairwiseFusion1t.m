function [ wV,...
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
                                                     refPointDist,...
                                                     map_env,...
                                                     map_recon,...
                                                     FoV,...
                                                     map_coverage,...
                                                     SensorRange,...
                                                     OBSenv,...
                                                     OBSconf,...
                                                     cellsize_env,...
                                                     visualize,...
                                                     para_ )
%
% 

tCompute_i = tic; % initialize pre computational time.


%******************************************************************** 
% SENSING PARAMETERS AND CONSTANTS:
%********************************************************************
FoVfaces.num = para_.numFoVs;   % number of conf per cell.
n            = numel(map_env);  % number of cells in the map.
f            = para_.numFoVs;   % for simplification.
c            = n*f;             % total conf.


%********************************************************************
% confs for fusion
%********************************************************************

fprintf(1,'----- confs for fusion \n')

%/////////////////////////////////
% -- first conf
%/////////////////////////////////
fusionConf1 = ConfPair_ioh(1,1);
fusionHcs1 = find(ConfPair_ioh(1,3:end));

fixedConfs1_ind = cell2mat(currentFusedConfs.ind(fusionHcs1));% selectedFusedConf_ind(fusionHcs1,:);
fixedConfs1_orn = cell2mat(currentFusedConfs.orn(fusionHcs1));

ind1 = find(fixedConfs1_ind==fusionConf1);
fixedConfs1_ind(ind1) = [];
fixedConfs1_orn(ind1) = [];

% 
redundantConf1_ind = fusionConf1;
% redundantConf1_orn = selectedFusedConf_orn(fusionHots1,ind1);
redundantConf1_orn = ConfPair_ioh(1,2);

% refPointDist1 = refPointDist(fusionHcs1,:);
refPointDist1 = refPointDist;

%//////////////////////////////////
% -- second conf
%//////////////////////////////////
fusionConf2 = ConfPair_ioh(2,1);
fusionHcs2 = find(ConfPair_ioh(2,3:end));

% fixedConfs2_ind = selectedFusedConf_ind(fusionHcs2,:);
% fixedConfs2_orn = selectedFusedConf_orn(fusionHcs2,:);

fixedConfs2_ind = cell2mat(currentFusedConfs.ind(fusionHcs2));
fixedConfs2_orn = cell2mat(currentFusedConfs.orn(fusionHcs2));


ind2 = find(fixedConfs2_ind==fusionConf2);
fixedConfs2_ind(ind2) = [];
fixedConfs2_orn(ind2) = [];

redundantConf2_ind = fusionConf2;
% redundantConf2_orn = selectedFusedConf_orn(fusionHots2,ind2);
redundantConf2_orn = ConfPair_ioh(2,2);

% refPointDist2 = refPointDist(fusionHcs2,:);
refPointDist2 = refPointDist;


% ----------------
redundantConfs_ind = [redundantConf1_ind,redundantConf2_ind];
redundantConfs_orn = [redundantConf1_orn,redundantConf2_orn];


% -- candidate conf map based on sensing range
% hotspots = hotspotSEQ([fusionHcs1,fusionHcs2]);
% [hotspots_r,hotspots_c] = ind2sub(size(map_env),hotspots);

hcs_fusion_num = unique([fusionHcs1,fusionHcs2]);
hcs_fusion_sub = hcs_sub(hcs_fusion_num,:);

hcs_fusion_r = hcs_fusion_sub(:,1);
hcs_fusion_c = hcs_fusion_sub(:,2);

hcs_fusion_ind = sub2ind(size(map_env),hcs_fusion_r,hcs_fusion_c);



%---- NEW
fusionConfs_ind = ConfPair_ioh(:,1);

fixedConfs_ind = cell(1,numel(hcs_fusion_num));
fixedConfs_orn = cell(1,numel(hcs_fusion_num));



for i = 1:numel(hcs_fusion_ind)
    %fixedConfs_ind{i} = selectedFusedConf_ind(hcs_fusion_num(i),:);
    %fixedConfs_orn{i} = selectedFusedConf_orn(hcs_fusion_num(i),:);
    
    fixedConfs_ind{i} = currentFusedConfs.ind{hcs_fusion_num(i)};
    fixedConfs_orn{i} = currentFusedConfs.orn{hcs_fusion_num(i)};
    
    %fixedConfs_ind{i}
    %fixedConfs_orn{i}
    
    ind1 = find(fixedConfs_ind{i}==fusionConfs_ind(1));
    ind2 = find(fixedConfs_ind{i}==fusionConfs_ind(2));
    
    ind = unique([ind1,ind2]);
    
    fixedConfs_ind{i}(ind) = [];
    fixedConfs_orn{i}(ind) = [];
end

% refPointDist
refPointsD = refPointDist;%(hcs_fusion_num,:);%= refPointDist(hcs_fusion_num,:);


%-- Placement of previously selected confs
engagedConfs_ind = cell2mat(currentFusedConfs.ind);
%=======================================================================
% Candidate Conf map based on sensing range
%=======================================================================
fprintf(1,'----- cand conf based on sensing range \n')

[ map_candConf ] = fConfMapSensingRange( hcs_fusion_ind,...
                                         map_env,...
                                         SensorRange,...
                                         OBSconf,...
                                         refPointDist,...
                                         para_,...
                                         visualize);



%=======================================================================
% Visibility matrix
%=======================================================================

fprintf(1,'----- visibility matrix \n')

[ V,...
  oV,...
  wV,...
  Conf_DistOrnGainV,...
  confKept,...
  fixedConfs_num,...
  fixedConfXG_num,... 
  confRemoved,...
  cellKept,...   
  cellRemoved,...
  oConfOrns,...
  FusionNotPossibleDueToInsufficientCandConf,...
  FoVfaces ] = fV_FusionHotspotGeometries( map_candConf,...
                                           map_coverage,...
                                           map_env,...
                                           map_recon,...
                                           engagedConfs_ind,...
                                           n,...
                                           FoV,...
                                           FoVfaces,...
                                           SensorRange,...
                                           fixedConfs_ind,...
                                           fixedConfs_orn,...
                                           fusionConf1,...
                                           fusionConf2,...
                                           OBSenv,...
                                           OBSconf,...
                                           hcs_fusion_ind,...
                                           cellsize_env,...
                                           para_,...
                                           visualize );

  
fixedConfsAll_num = unique(cell2mat(fixedConfs_num));


if FusionNotPossibleDueToInsufficientCandConf == 1
    D                   = [];
    NumOfFixeConf       = [];
    hcs_fusion_ind      = [];
    hcs_fusion_num      = [];
    conf_crossAngles_G  = [];
    return
end


%=======================================================================
% -- CANDIDATE CONF MAP based on sensing range and visibility matrix
%=======================================================================
fprintf(1,'----- cand confs map \n')

% initialize final candidate conf map
map_candConf = zeros(size(map_env));
confCELL = confKept; % cell num for each conf num.

% update the map
map_candConf(confCELL) = 1;


% -- plot 
if visualize == 1
figure('name','cand conf map - FINAL');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])

map_candConf_ind = find(map_candConf(:)==1);
[map_candConf_r,map_candConf_c] = ind2sub(size(map_env),map_candConf_ind);
plot(map_candConf_r,map_candConf_c,'xb','MarkerFaceColor','b');
plot(hcs_fusion_r,hcs_fusion_c,'or','MarkerFaceColor','r');

% fixedConfs_ind = [fixedConfs1_ind,fixedConfs2_ind];
[fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),cell2mat(fixedConfs_ind));
plot(fixedConfs_r,fixedConfs_c,'or');

fusionConfs_ind = [fusionConf1,fusionConf2];
[fusionConfs_r,fusionConfs_c] = ind2sub(size(map_env),fusionConfs_ind);
plot(fusionConfs_r,fusionConfs_c,'og');
end

%=======================================================================
% CROSS ANGLES
%=======================================================================

fprintf(1,'----- cross angles \n')

[ Hot2Conf_Angles,...
  conf_crossAngles,...
  conf_crossAngles_G ] = fX( confKept,...
                             fixedConfXG_num,...
                             localNumOfSelectedConfs,...
                             hcs_fusion_ind,...
                             hcs_fusion_num,...
                             map_env,...
                             para_ );


%=======================================================================
% DISTANCE 
%=======================================================================

fprintf(1,'----- traveling distance \n')

[ D ] = fD( refPointsD,...
            fixedConfXG_num,...
            hcs_fusion_ind,...
            confKept,...
            map_env );



%************************************************
%  Conf num and Number of fixed confs
%************************************************ 

% fixedConf1_num = fixedConf1_num';
% fixedConf2_num = fixedConf2_num';
% 
% NumOfFixeConf = numel(unique([fixedConf1_num;fixedConf2_num]));
NumOfFixeConf = numel(unique(cell2mat(fixedConfs_num)));



%**********************************************************************
% Total computation time.
%**********************************************************************
tCompute = 1e-4*(round(toc(tCompute_i)*1e4));

fprintf(1,'---- Computation time: %0.4f sec \n',tCompute)




end



function [ D ] = fD( refPointsD,...
                     fixedConfXG_num,...
                     hcs_fusion_ind,...
                     confKept,...
                     map_env )



% ---- D (distance) matrix is manhattan distance between cand conf position
% and reference position
confKeptCells_ind = confKept; %fix(confKept/f)+(~~mod(confKept,f));
[confKeptCells_r,confKeptCells_c] = ind2sub(size(map_env),confKeptCells_ind);

D = cell(1,numel(hcs_fusion_ind)); %zeros(numel(confKeptCells_ind),1);

for i = 1:numel(hcs_fusion_ind)
    for j = 1:numel(confKeptCells_ind)
        
        %manhattan_distance = abs(confKeptCells_r(j)-refPointsD(i,1))+...
        %    abs(confKeptCells_c(j)-refPointsD(i,2));
        manhattan_distance = abs(confKeptCells_r(j)-refPointsD(1))+...
            abs(confKeptCells_c(j)-refPointsD(2));
        
        D{i}(j,1) = manhattan_distance;
    end
end
   

% Distance values to the fixed conf of the other hotspot is max.
for i = 1:numel(hcs_fusion_ind)
    D{i}(fixedConfXG_num{i},1) = max(D{i}(:));
end
   


end



function [ Hot2Conf_Angles,...
           conf_crossAngles,...
           conf_crossAngles_G ] = fX( confKept,...
                                      fixedConfXG_num,...
                                      localNumOfSelectedConfs,...
                                      hcs_fusion_ind,...
                                      hcs_fusion_num,...
                                      map_env,...
                                      para_ )

%
%


%***********************************************************************
% CONF ANGLES TOWARDS HCs
%***********************************************************************
fprintf('------ angles towards Hcs...\n')

confKeptCell_ind = confKept;
[confKeptCell_r,confKeptCell_c] = ind2sub(size(map_env),confKeptCell_ind);
sensors = [confKeptCell_r',confKeptCell_c'];

[hcs_fusion_r,hcs_fusion_c] = ind2sub(size(map_env),hcs_fusion_ind);

Hot2Conf_Angles = zeros(numel(hcs_fusion_ind),numel(confKept));
for i = 1:numel(hcs_fusion_ind)
    
    %Hot2Conf_Angles(i,:) = rad2deg(angle2Points([hotspots_r(i),hotspots_c(i)],sensors));
    Hot2Conf_Angles(i,:) = angle2Points([hcs_fusion_r(i),hcs_fusion_c(i)],sensors);
end



%***********************************************************************
% CROSS ANGLES
%***********************************************************************
disp('------ cross angles')

conf_crossAngles = cell(1,numel(hcs_fusion_ind));
for i = 1:numel(hcs_fusion_ind)
    for j = 1:numel(confKept)
        for k = 1:numel(confKept)
            conf_crossAngles{i}(j,k) = ...
                rad2deg(angleAbsDiff(Hot2Conf_Angles(i,j),Hot2Conf_Angles(i,k)));
        end
    end
end


%***********************************************************************
% CROSS ANGLES GAIN
%***********************************************************************
disp('------ cross angles gain')


conf_crossAngles_G = cell(1,numel(hcs_fusion_ind));

for h = 1:numel(hcs_fusion_ind)
    
    NumOfAllowedConf = localNumOfSelectedConfs(hcs_fusion_num(h));

    if     NumOfAllowedConf == 2
        gain_meu = para_.xGainMeu2Conf_DEG;
        gain_sigma = para_.xGainSig2Conf_DEG;

    elseif NumOfAllowedConf == 3
        gain_meu = para_.xGainMeu3Conf_DEG;
        gain_sigma = para_.xGainSig3Conf_DEG;

    elseif NumOfAllowedConf == 4
        gain_meu = para_.xGainMeu4Conf_DEG;
        gain_sigma = para_.xGainSig4Conf_DEG;

    elseif NumOfAllowedConf == 5
        gain_meu = para_.xGainMeu5Conf_DEG;
        gain_sigma = para_.xGainSig5Conf_DEG;
    else
        error('Unknown cross angle gains')
    end

    guss_gain = zeros(size(conf_crossAngles{h}));
    for i = 1:numel(gain_meu)
        guss_gain = ...
            max(guss_gain,gaussmf(conf_crossAngles{h},[gain_sigma,gain_meu(i)]));
    end
    conf_crossAngles_G{h} = guss_gain;
    
end




%%***********************************************************************
%
% Special case: gains for the fixed conf
%
%***********************************************************************
disp('------ gain for fixed confs')
% -- gain between fixed confs of hotspot 1 and 2 are zero

for i = 1:numel(hcs_fusion_ind)
    conf_crossAngles_G{i}(:,fixedConfXG_num{i}) = 0;
end



end




function [ map_candConf ] = fConfMapSensingRange( hc_fusion_ind,...
                                                  map_env,...
                                                  SensorRange,...
                                                  OBSconf,...
                                                  refPointDist,...
                                                  para_,...
                                                  visualize)

% [hotspots_r,hotspots_c] = ind2sub(size(map_env),hc_current);
[hc_fusion_r,hc_fusion_c] = ind2sub(size(map_env),hc_fusion_ind);

% initialize map for candidate conf. -- ring only
map_candConf = zeros(size(map_env));    
map_candConf(hc_fusion_ind) = 1;

% -- sample candidate conf within circle radius of sensing range
% candCircleThickness_m = 0 % thickness of candidate conf circle (meters)
% candCircleThickness_cell = candCircleThickness_m/paramPreprocess.MapRes
radius_candCircleOut = ceil(SensorRange.cell);


% grow hotspot to max radius
se = strel('disk',radius_candCircleOut,0); % outer limit
map_candConf = imdilate(map_candConf,se);

% -- plot
if visualize == 1
figure('name','cand conf - max limit and obs are removed');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
map_candConf_ind = find(map_candConf(:)==1);
[map_candConf_r,map_candConf_c] = ind2sub(size(map_env),map_candConf_ind);
plot(map_candConf_r,map_candConf_c,'xb','MarkerFaceColor','b');        
plot(hc_fusion_r,hc_fusion_c,'xr','MarkerFaceColor','r');
end

% -- remove obs cells
map_candConf(OBSconf.ind(ismember(OBSconf.ind,find(map_candConf)))) = 0;

if visualize == 1
% plot
figure('name','cand conf - max limit and obs are removed');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
map_candConf_ind = find(map_candConf(:)==1);
[map_candConf_r,map_candConf_c] = ind2sub(size(map_env),map_candConf_ind);
plot(map_candConf_r,map_candConf_c,'xb','MarkerFaceColor','b');        
plot(hc_fusion_r,hc_fusion_c,'xr','MarkerFaceColor','r');
end


%-- Limit the number of candidate confs before computing visibility (1SE only)
%==============================================================================

% refConfCell = refPointDist(1,1:2);

candConf_i = find(map_candConf);

maxConfs = para_.maxCandidateConf_Fusion_1t+para_.deltaCandConf_Fusion_1t;
% if numel(candConf_i) > (para_.NumOfCandidateConfTomography1SE+20)
if numel(candConf_i) > maxConfs

    
    %[candConf_r,candConf_c] = ind2sub(size(map_env),candConf_i);

    %confdist = pdist2([candConf_r,candConf_c],refConfCell,'euclidean');
        
    % -- sorted indices
    %[~,sorted_ind] = sort(confdist,'ascend');
    %sorted_ind_max = sorted_ind(1:min([para_.NumOfCandidateConfTomography1SE,numel(sorted_ind)]));
    
    rnd_series = randperm(numel(candConf_i),numel(candConf_i));
    [~,sorted_ind] = sort(rnd_series,'ascend');
    
    %sorted_ind_max = sorted_ind(1:min([para_.NumOfCandidateConfTomography1SE+20,...
    %    numel(sorted_ind)]));
    
    sorted_ind_max = sorted_ind(1:min([maxConfs,numel(sorted_ind)]));
        
    candConfKept_i = candConf_i(sorted_ind_max);
    
    map_candConf = zeros(size(map_env));
    map_candConf(candConfKept_i) = 1;
    
end

% % extract edges
% map_candConf_edges = bwmorph(map_candConf,'remove');
% % -- remove edges
% map_candConf(map_candConf_edges(:)==1) = 0;
% 
% 
% % plot
% if visualize == 1
% pause(0.1)
% figure('name','cand conf - extracted edges'); 
% imshow(map_env','InitialMagnification',800); hold on;
% set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
% map_candConf_ind = find(map_candConf(:)==1);
% [map_candConf_r,map_candConf_c] = ind2sub(size(map_env),map_candConf_ind);
% plot(map_candConf_r,map_candConf_c,'xb','MarkerFaceColor','b');        
% plot(hc_r,hc_c,'xr','MarkerFaceColor','r');
% end

end




%--------------------------------- End of the document ----------------------------------


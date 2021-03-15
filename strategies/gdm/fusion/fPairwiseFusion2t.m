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
           conf_crossAngles_G ] = fPairwiseFusion2t( ConfPair_ioh,...
                                                     currentFusedConfs,...
                                                     localNumOfSelectedConfs,...
                                                     refPointDist,...
                                                     hcsSEQ_ind,...
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

% /////////////////////////////////
% -- first conf
% /////////////////////////////////
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

% refPointDist1 = refPointDist(fusionHots1,:);

% //////////////////////////////////
% -- second conf
% //////////////////////////////////
fusionConf2 = ConfPair_ioh(2,1);
fusionHcs2 = find(ConfPair_ioh(2,3:end));

% fixedConfs2_ind = selectedFusedConf_ind(fusionHots2,:);
% fixedConfs2_orn = selectedFusedConf_orn(fusionHots2,:);

fixedConfs2_ind = cell2mat(currentFusedConfs.ind(fusionHcs2));
fixedConfs2_orn = cell2mat(currentFusedConfs.orn(fusionHcs2));

ind2 = find(fixedConfs2_ind==fusionConf2);
fixedConfs2_ind(ind2) = [];
fixedConfs2_orn(ind2) = [];

redundantConf2_ind = fusionConf2;
% redundantConf2_orn = selectedFusedConf_orn(fusionHots2,ind2);
redundantConf2_orn = ConfPair_ioh(2,2);

% refPointDist2 = refPointDist(fusionHots2,:);


% ----------------
redundantConfs_ind = [redundantConf1_ind,redundantConf2_ind];
redundantConfs_orn = [redundantConf1_orn,redundantConf2_orn];


% -- candidate conf map based on sensing range
% hotspots = hotspotSEQ([fusionHots1,fusionHots2]);
hcs_fusion_num = unique([fusionHcs1,fusionHcs2]);
hcs_fusion_ind = hcsSEQ_ind(hcs_fusion_num);


[hcs_fusion_r,hcs_fusion_c] = ind2sub(size(map_env),hcs_fusion_ind);

if numel(hcs_fusion_ind)>2
    disp('more than two hc detected for fusion')
    %hc_fusion_ind
    %pause
end

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
refPointsD = refPointDist(hcs_fusion_num,:);


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
        
        manhattan_distance = abs(confKeptCells_r(j)-refPointsD(i,1))+...
            abs(confKeptCells_c(j)-refPointsD(i,2));
        
        D{i}(j,1) = manhattan_distance;
    end
end
   

% Distance values to the fixed conf of the other hotspot is max.
for i = 1:numel(hcs_fusion_ind)
    D{i}(fixedConfXG_num{i},1) = max(D{i}(:));
end
   


%{                    

refDistPoint1_r = refPointDistance1(1); refDistPoint1_c = refPointDistance1(2);
refDistPoint2_r = refPointDistance2(1); refDistPoint2_c = refPointDistance2(2);


% ---- D (distance) matrix is manhattan distance between cand conf position
% and reference position
confKeptCells_ind = confKept; %fix(confKept/f)+(~~mod(confKept,f));
[confKeptCells_r,confKeptCells_c] = ind2sub(size(map_env),confKeptCells_ind);

D1 = zeros(numel(confKeptCells_ind),1);
D2 = zeros(numel(confKeptCells_ind),1);

for i = 1:numel(confKeptCells_ind)
    manhattan_distance1 = abs(confKeptCells_r(i)-refDistPoint1_r)+...
        abs(confKeptCells_c(i)-refDistPoint1_c);
    manhattan_distance2 = abs(confKeptCells_r(i)-refDistPoint2_r)+...
        abs(confKeptCells_c(i)-refDistPoint2_c);
    
    D1(i) = manhattan_distance1;
    D2(i) = manhattan_distance2;
end


% Distance values to the fixed conf of the other hotspot is max.
% D1(fixedConf2_num) = max(D1(:));
% D2(fixedConf1_num) = max(D2(:));
D1(fixedConf2XG_num) = max(D1(:));
D2(fixedConf1XG_num) = max(D2(:));

%}

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

% conf_theta            := conf orientations [1,c]
% conf_crossAngles      := cross angles between conf [c,c]
% conf_crossAngles_G    := cross angles gain [c,c]


%***********************************************************************
% CONF ANGLES TOWARDS HCs
%***********************************************************************
% fprintf('angles towards Hcs...\n')

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
disp('----- cross angles')

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
disp('----- cross angles gain')



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


%*****************************************************************
%       Gain - cross angles
%*****************************************************************

%{
% ------ cross angles gain 1
% disp('----- gain')
% gain_sigma = deg2rad(10); % standard deviation
% gain_meu   = deg2rad(60); % expected value or center
gain_sigma = para_.xGainSigma_DEG; % standard deviation
gain_meu   = para_.xGainMeu1_DEG; % expected value or center
conf_crossAngles1_G1 = gaussmf(conf_crossAngles1,[gain_sigma,gain_meu]);


% gain_sigma = deg2rad(10); % standard deviation
% gain_meu   = deg2rad(120); % expected value or center
gain_sigma = para_.xGainSigma_DEG; % standard deviation
gain_meu   = para_.xGainMeu2_DEG; % expected value or center
conf_crossAngles1_G2 = gaussmf(conf_crossAngles1,[gain_sigma,gain_meu]);

conf_crossAngles1_G = conf_crossAngles1_G1+conf_crossAngles1_G2;



% ------ cross angles gain 2
% disp('----- gain')
% gain_sigma = deg2rad(10); % standard deviation
% gain_meu   = deg2rad(60); % expected value or center
gain_sigma = para_.xGainSigma_DEG; % standard deviation
gain_meu   = para_.xGainMeu1_DEG; % expected value or center
conf_crossAngles2_G1 = gaussmf(conf_crossAngles2,[gain_sigma,gain_meu]);

% gain_sigma = deg2rad(10); % standard deviation
% gain_meu   = deg2rad(120); % expected value or center
gain_sigma = para_.xGainSigma_DEG; % standard deviation
gain_meu   = para_.xGainMeu2_DEG; % expected value or center
conf_crossAngles2_G2 = gaussmf(conf_crossAngles2,[gain_sigma,gain_meu]);

conf_crossAngles2_G = conf_crossAngles2_G1+conf_crossAngles2_G2;
 
%}


%*****************************************************************
%       Special case: gains for the fixed conf
%*****************************************************************
%{
% -- gain between fixed confs of hotspot 1 and 2 are zero

% fixedConf1_num = find(ismember(confKept,fixedConfs1_ind));
% fixedConf2_num = find(ismember(confKept,fixedConfs2_ind));


% conf_crossAngles1_G(fixedConf1_num,fixedConf2_num) = 0;
% conf_crossAngles1_G(fixedConf2_num,fixedConf1_num) = 0;

% conf_crossAngles1_G(:,fixedConf2_num) = 0;
% conf_crossAngles1_G(fixedConf2_num,:) = 0;

conf_crossAngles1_G(:,fixedConf2XG_num) = 0;
conf_crossAngles1_G(fixedConf2XG_num,:) = 0;


% conf_crossAngles2_G(fixedConf1_num,fixedConf2_num) = 0;
% conf_crossAngles2_G(fixedConf2_num,fixedConf1_num) = 0;

% conf_crossAngles2_G(:,fixedConf1_num) = 0;
% conf_crossAngles2_G(fixedConf1_num,:) = 0;

conf_crossAngles2_G(:,fixedConf1XG_num) = 0;
conf_crossAngles2_G(fixedConf1XG_num,:) = 0;

% conf_crossAngles1_G = conf_crossAngles1_G./(max(conf_crossAngles1_G(:)))
% conf_crossAngles2_G = conf_crossAngles2_G./(max(conf_crossAngles2_G(:)))
%}


%***********************************************************************
%
% Special case: gains for the fixed conf
%
%***********************************************************************
disp('----- gain for fixed confs')
% -- gain between fixed confs of hotspot 1 and 2 are zero

for i = 1:numel(hcs_fusion_ind)
    conf_crossAngles_G{i}(:,fixedConfXG_num{i}) = 0;
end


end




function [ map_candConf ] = fConfMapSensingRange( hc_fusion_ind,...
                                                  map_env,...
                                                  SensorRange,...
                                                  OBSconf,...
                                                  para_,...
                                                  visualize)

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


%                       TMP

%-- Limit the number of candidate confs before computing visibility (1SE only)
%==============================================================================

% refConfCell = refPointDist(1,1:2);

candConf_i = find(map_candConf);

maxConfs = para_.maxCandidateConf_Fusion_2t+para_.deltaCandConf_Fusion_2t;
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
% plot(hc_fusion_r,hc_fusion_c,'xr','MarkerFaceColor','r');
% end

end


% function [ V,...
%            oV,...
%            wV,...
%            Conf_DistOrnGainV,...
%            confKept,...
%            fixedConfs_num,...
%            fixedConfXG_num,... 
%            confRemoved,...
%            cellKept,...   
%            cellRemoved,...
%            oConfOrns,...
%            FoVfaces ] = fV( map_candConf,...
%                             map_coverage,...
%                             map_env,...
%                             map_recon,...
%                             engagedConfs_ind,...
%                             n,...
%                             f,...
%                             FoV,...
%                             FoVfaces,...
%                             SensorRange,...
%                             fixedConfs_ind,...
%                             fixedConfs_orn,...
%                             fusionConf1,...
%                             fusionConf2,...
%                             OBSenv,...
%                             OBSconf,...
%                             hcs_fusion_ind,...
%                             gamma,...
%                             cellsize_env,...
%                             para_,...
%                             visualize )
% %fV generates Visibility matrix V.
% 
% % --- ANGLES:
% FoVfaces.lgt = [1/FoVfaces.num:1/FoVfaces.num:1-(1/FoVfaces.num) 0]'*360;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                                    SensorDomain
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % disp('.. Sensor Domain')
% 
% switch para_.SensorDomainComputation
%     
%     case 'ComputeNew'
%         
%         % For each cell of the map.
%         sD = cell(1,f); % sensor domain (visible).
%         oD = cell(1,f); % sensor's osbtacle domain (used to check the 
%                         % visibility of "sD" for each cell in "oD"]
% 
%         % a single reference cell in the map.
%         for i = 1
%             
%             % sensor position.
%             sensor = [0,0];
%             
%             % For each sensing configuration.
%             for j = 1:f
% 
%                 % zawia := initial, middle and final (opening) angle.
%                 zawia = [FoVfaces.ang(j,1) FoVfaces.lgt(j) FoVfaces.ang(j,2)];
% 
%                 % Finding sensor domain, i.e. visible cells for the
%                 % configuration. 
%                 [ visibleRange ] = fSenDomain( sensor,SensorRange,zawia );
%                 % visibleRange := visible range (row,col).
% 
%                 % post-processing for sorting and removing duplicates
%                 visibleRange = unique(visibleRange,'rows'); % remove repeated numbers
%                 visibleRange = sortrows(visibleRange);      % sort rows
%                 
%                 % here we have jth sensor's domain.
%                 sD{j} = visibleRange;
%                 
%                 
%                 % since, to varify if a cell in sensor's domain is visible,
%                 % we need to check all the cells that effect the
%                 % visibility. Some of them are even not in the sensing
%                 % range, therefore, we need to compute obstacle domain that
%                 % contains the cells in the sensor domain and connecting.
%                 % "obsDomain" will be used to check visibility of the
%                 % sensor's domain.
% 
%                 % obstacle domain
%                 [ obsRange ] = fSenDomainOBS( sensor,SensorRange,zawia );
%                 % obsRange := obstacle range (row,col).
%                 
%                 % post-processing for sorting and removing duplicates
%                 obsRange = unique(obsRange,'rows'); % remove repeated numbers
%                 obsRange = sortrows(obsRange);      % sort rows
%                 
%                 % here we have jth sensor's domain.
%                 oD{j} = obsRange;
%                 
%             end
%         end
% 
%         %--------------------------------------
%         % EFFECTED VISIBILITY BY EACH OBSTACLE
%         %--------------------------------------
%         
%         disp('.. Visibility vs Obstacles')
%         VObs = cell(1,f); % Visibility vs Obstacle
%         
%         for i = 1
%             sensor = [0,0];
%             % For each configuration of a cell.
%             for j = 1:f
%                 
%                 % let take sensor's mask
%                 sensorDom  = sD{j};
%                 % and obstacle domain
%                 obstacDom  = oD{j};
%                 
%                 VObs{j} = zeros(size(obstacDom,1),size(sensorDom,1)); % no more a square matrix
%                 % rows    := obstacle cell,
%                 % columns := visibility effect.
%                 
%                 %for k = 1:size(sensorDom,1)
%                 for k = 1:size(obstacDom,1)
%                     
%                     % kth cell is an obstacle
%                     %obstacDomOBS(k,3)  = 0; 
%                     obs = obstacDom(k,:);
%                     
%                     % visibility effect on all the cells
%                     [ vDom ] = fVisSDvObs( sensor,sensorDom,obs );
%                     
%                     % update visibility against kth obstacle
%                     VObs{j}(k,:) = vDom;
%                    
%                 end
%                
%             end
%         end
%         
%         % save the mask and visibility info to a file.
%         %save(['preprocess/sD_oD_VObs_Tomography_Rng',...
%         %    num2str(SensorRange.cell),'FoV',num2str(FoV)],...
%         %    'sD','oD','VObs');
%         
%         filename = sprintf('preprocess/sD_oD_VObs_Tomography_range%02d_fov%03d_fovfaces%03d.mat',...
%                             SensorRange.cell,FoV,f);
%         save(filename,'sD','oD','VObs');
%         
%                 
%         
%         
%         
%     case 'UseEarlier'
%         
%         %filename = sprintf('preprocess/sD_oD_VObs_Tomography_range%02d_fov%03d_fovfaces%03d.mat',...
%         %                    SensorRange.cell,FoV,f)
%         %load(filename);
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                                       V
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % disp('.. Visibility Matrix')
% % V = sparse(n,c); % initializing
% V = sparse(n,n); % initializing
% % V = zeros(n,c); % initializing
% 
% [hcs_fusion_r,hcs_fusion_c] = ind2sub(size(map_env),hcs_fusion_ind);
% 
% oConfOrns = zeros(1,n);
% % For each cell of the map.
% for i = 1:n
%     if map_candConf(i) % if its not occupied
%         [rr,cc] = ind2sub(size(map_candConf),i);
%         sensor = [rr,cc];
%         
%         %-------------------------------------------------------------
%         %           Orientation of this conf
%         %-------------------------------------------------------------
%         % --
%         %{
%         orientation1 = round(rad2deg(angle2Points(sensor,[hotspots_r(1),hotspots_c(1)])));
%         orientation1 = orientation1 + (orientation1==0)*360;
%         
%         orientation2 = round(rad2deg(angle2Points(sensor,[hotspots_r(2),hotspots_c(2)])));
%         orientation2 = orientation2 + (orientation2==0)*360;
%         
%         %orientation =
%         %round(wrapTo360(meanangle([orientation1,orientation2]))) % this has problem 
%         %with [270,90] or similar cases
%         orientation = round(wrapTo360(rad2deg(circ_mean([deg2rad(orientation1),...
%                                                          deg2rad(orientation2)],[],2))));
%         orientation = orientation + (orientation==0)*360;
%         %}
%         
%         orns = zeros(1,numel(hcs_fusion_ind));
%         for k = 1:numel(hcs_fusion_ind)
%             
%             orn_this = round(rad2deg(angle2Points(sensor,[hcs_fusion_r(k),hcs_fusion_c(k)])));
%             orn_this = orn_this + (orn_this==0)*360;
%             
%             orns(k) = orn_this;
%         end
%         orientation = round(wrapTo360(rad2deg(circ_mean(deg2rad(orns),[],2))));
%         orientation = orientation + (orientation==0)*360;
%                 
%         
%         %orientation = orientation2; % for debugging
%         % For each configuration of a cell.
%         for j = orientation %1%:f
%             
%             sD = []; oD = []; VObs = [];
%             this_dir = sprintf(['sD_oD_VObs/sensing_range_%02dcell/fov_%03ddeg/',...
%                                 'orientation_%03ddeg/'],SensorRange.cell,FoV,j);
%             this_file = sprintf('sD_oD_VObs_rng%02d_fov%03d_orn%03d.mat',...
%                         SensorRange.cell,FoV,j);
%             load([this_dir,this_file],'sD','oD','VObs');
%             
%             % sensorDom := sensor domain translation w.r.t sensor position.
%             sensorDom = [sD(:,1)+sensor(1) sD(:,2)+sensor(2)];
%             
%             % obstacDom := obstalce domain translated w.r.t sensor's
%             % position. 
%             obstacDom = [oD(:,1)+sensor(1) oD(:,2)+sensor(2)];
%             
%             % VOBS := Visibility vs Obstacles matrix for jth FoV.
%             VOBS = VObs;
%             
%             % visibile cells 
%             [ visDom ] = fVisibileDomain( map_env,sensorDom,obstacDom,VOBS );
%             
%             % Update Visibility matrix.
%             %V(visDom,((i-1)*f)+j) = 1;
%             V(visDom,i) = 1;
%             
%         end
%         %{
%         % For each configuration of a cell.
%         for j = orientation %1%:f
%             
%             % sensorDom := sensor domain translation w.r.t sensor position.
%             sensorDom = [sD{j}(:,1)+sensor(1) sD{j}(:,2)+sensor(2)];
%             
%             % obstacDom := obstalce domain translated w.r.t sensor's
%             % position. 
%             obstacDom = [oD{j}(:,1)+sensor(1) oD{j}(:,2)+sensor(2)];
%             
%             % VOBS := Visibility vs Obstacles matrix for jth FoV.
%             VOBS = VObs{j};
%             
%             % visibile cells 
%             [ visDom ] = fVisibileDomain( map_coverage,map_env,sensorDom,obstacDom,VOBS );
%             
%             % Update Visibility matrix.
%             %V(visDom,((i-1)*f)+j) = 1;
%             V(visDom,i) = 1;
%             
%         end
%         %}
%         oConfOrns(i) = orientation;
%         
%     end
% end
% 
% %{
% for i = 1:numel(fixedConfs1_ind)
%             
%     [rr,cc] = ind2sub(size(map_candConf),fixedConfs1_ind(i));
%     sensor  = [rr,cc];
%         
%     j = fixedConfs1_orn(i);
%     
%     %{
%     % sensorDom := sensor domain translation w.r.t sensor position.
%     sensorDom = [sD{j}(:,1)+sensor(1) sD{j}(:,2)+sensor(2)];
% 
%     % obstacDom := obstalce domain translated w.r.t sensor's
%     % position.
%     obstacDom = [oD{j}(:,1)+sensor(1) oD{j}(:,2)+sensor(2)];
% 
%     % VOBS := Visibility vs Obstacles matrix for jth FoV.
%     VOBS = VObs{j};
% 
%     %}
%     
%     sD = []; oD = []; VObs = [];
%     this_dir = sprintf(['sD_oD_VObs/sensing_range_%02dcell/fov_%03ddeg/',...
%                         'orientation_%03ddeg/'],SensorRange.cell,FoV,j);
%     this_file = sprintf('sD_oD_VObs_rng%02d_fov%03d_orn%03d.mat',...
%                 SensorRange.cell,FoV,j);
%     load([this_dir,this_file],'sD','oD','VObs');
%     
%     % sensorDom := sensor domain translation w.r.t sensor position.
%     sensorDom = [sD(:,1)+sensor(1) sD(:,2)+sensor(2)];
% 
%     % obstacDom := obstalce domain translated w.r.t sensor's
%     % position.
%     obstacDom = [oD(:,1)+sensor(1) oD(:,2)+sensor(2)];
% 
%     % VOBS := Visibility vs Obstacles matrix for jth FoV.
%     VOBS = VObs;
%     
%     % visibile cells 
%     [ visDom ] = fVisibileDomain( map_env,sensorDom,obstacDom,VOBS );
%     
%     % Update Visibility matrix.
%     %V(visDom,((i-1)*f)+j) = 1;
%     V(visDom,fixedConfs1_ind(i)) = 1;
%     
%     oConfOrns(fixedConfs1_ind(i)) = j;
% end
%    
% 
% for i = 1:numel(fixedConfs2_ind)
%     
%     [rr,cc] = ind2sub(size(map_candConf),fixedConfs2_ind(i));
%     sensor  = [rr,cc];
%     
%     j = fixedConfs2_orn(i);
% 
%     %{
%     % sensorDom := sensor domain translation w.r.t sensor position.
%     sensorDom = [sD{j}(:,1)+sensor(1) sD{j}(:,2)+sensor(2)];
% 
%     % obstacDom := obstalce domain translated w.r.t sensor's
%     % position.
%     obstacDom = [oD{j}(:,1)+sensor(1) oD{j}(:,2)+sensor(2)];
% 
%     % VOBS := Visibility vs Obstacles matrix for jth FoV.
%     VOBS = VObs{j};
% 
%     %}
%     sD = []; oD = []; VObs = [];
%     this_dir = sprintf(['sD_oD_VObs/sensing_range_%02dcell/fov_%03ddeg/',...
%                         'orientation_%03ddeg/'],SensorRange.cell,FoV,j);
%     this_file = sprintf('sD_oD_VObs_rng%02d_fov%03d_orn%03d.mat',...
%                 SensorRange.cell,FoV,j);
%     load([this_dir,this_file],'sD','oD','VObs');
%     
%     % sensorDom := sensor domain translation w.r.t sensor position.
%     sensorDom = [sD(:,1)+sensor(1) sD(:,2)+sensor(2)];
% 
%     % obstacDom := obstalce domain translated w.r.t sensor's
%     % position.
%     obstacDom = [oD(:,1)+sensor(1) oD(:,2)+sensor(2)];
% 
%     % VOBS := Visibility vs Obstacles matrix for jth FoV.
%     VOBS = VObs;
%     
%     % visibile cells 
%     [ visDom ] = fVisibileDomain( map_env,sensorDom,obstacDom,VOBS );
% 
%     % Update Visibility matrix.
%     %V(visDom,((i-1)*f)+j) = 1;
%     V(visDom,fixedConfs2_ind(i)) = 1;
% 
%     oConfOrns(fixedConfs2_ind(i)) = j;
% end
% %}
% 
% 
% for i = 1:numel(hcs_fusion_ind)
%     
%     for j = 1:numel(fixedConfs_ind{i})
%         
%         [rr,cc] = ind2sub(size(map_candConf),fixedConfs_ind{i}(j));
%         sensor  = [rr,cc];
%         
%         k = fixedConfs_orn{i}(j);
%         
%         sD = []; oD = []; VObs = [];
%         this_dir = sprintf(['sD_oD_VObs/sensing_range_%02dcell/fov_%03ddeg/',...
%                             'orientation_%03ddeg/'],SensorRange.cell,FoV,k);
%         this_file = sprintf('sD_oD_VObs_rng%02d_fov%03d_orn%03d.mat',...
%                     SensorRange.cell,FoV,k);
%         load([this_dir,this_file],'sD','oD','VObs');
%         
% 
%         % sensorDom := sensor domain translation w.r.t sensor position.
%         sensorDom = [sD(:,1)+sensor(1) sD(:,2)+sensor(2)];
% 
%         % obstacDom := obstalce domain translated w.r.t sensor's
%         % position.
%         obstacDom = [oD(:,1)+sensor(1) oD(:,2)+sensor(2)];
% 
%         % VOBS := Visibility vs Obstacles matrix for jth FoV.
%         VOBS = VObs;
%     
%         % visibile cells 
%         [ visDom ] = fVisibileDomain( map_env,sensorDom,obstacDom,VOBS );
%     
%         % Update Visibility matrix.        
%         V(visDom,fixedConfs_ind{i}(j)) = 1;
%         
%         oConfOrns(fixedConfs_ind{i}(j)) = k;
%     end
% end
% 
% 
% % sensing config can not visualize own cell
% for i = 1:numel(map_candConf)    
%     %V(i,(i-1)*f+1:(i-1)*f+f) = 0;
%     V(i,i) = 0;
% end
% 
% oV = V; % oV := visibility with obstacles
% 
% % VISIBILITY MATRIX CORRESPOND TO THE HOTSPOT ONLY:
% % ------------------------------------------------------------------------------
% 
% %////////////////////////////////////////
% %
% % CELLS TO RETAIN
% %
% %////////////////////////////////////////
% 
% %----------------------------------------------------
% % -- cells covered by candidate conf.
% %----------------------------------------------------
% covered_ind = sum(oV,2);
% covered_num = find(covered_ind);
% covered_cell_vector = zeros(1,size(oV,1));
% covered_cell_vector(covered_num) = 1;
% 
% %----------------------------------------------------
% % -- cells need to be covered -- according to the coverage map
% %----------------------------------------------------
% coverage_ind = find(map_coverage);
% coverage_cell_vector = zeros(1,size(oV,1));
% % size(coverage_cell_vector)
% coverage_cell_vector(coverage_ind) = 1;
% % size(coverage_cell_vector)
% 
% %----------------------------------------------------
% % -- unoccupied cells
% %----------------------------------------------------
% free_cell_vector = ones(1,size(oV,1));
% free_cell_vector(OBSenv.ind) = 0;
% 
% 
% % --- selecting cells that are covered and need to be covered and unoccupied.
% %***************************************************************************
% %
% % cells kept to observe are: 
% %       (1) can be covered by cand conf (Check#1) AND 
% %       (2) desired to be covered (Check#2) AND 
% %       (3) are not occupied (Check#3) 
% %
% %***************************************************************************
% cov_cell     = covered_cell_vector.*coverage_cell_vector.*free_cell_vector;
% cellKept     = find(cov_cell>=1);
% cellRemoved  = find(cov_cell==0);
% V            = V(cellKept,:);
% 
% 
% % ////////////////////////////
% % -- cellKept weights
% % ////////////////////////////
% % cellKeptWeights = map_gd(cellKept);
% % cellKeptWeights = cellKeptWeights/(max(cellKeptWeights));
% 
% 
% % //////////////////////////
% % ---- configurations:
% % //////////////////////////
% 
% %----------------------------------------------------
% % -- Check#1: conf corresponding to candidate conf map only
% %----------------------------------------------------
% % cells for the candidate conf.
% candConfCells = find(map_candConf); 
% % conf numbers
% candConfNum = zeros(1,size(oV,2));
% for i = 1:numel(candConfCells)
%     %candConfNum( (candConfCells(i)-1)*f+1:(candConfCells(i)-1)*f+f ) = 1;
%     candConfNum( candConfCells(i) ) = 1;
% end
% % - fixed confs are in the vector
% % candConfNum([fixedConfs1_ind,fixedConfs2_ind]) = 1;
% candConfNum(cell2mat(fixedConfs_ind)) = 1;
% 
% 
% % -- plot 
% if visualize == 1
% figure('name','conf corresponding to candidate conf map only');
% imshow(map_env','InitialMagnification',800); hold on;
% set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
% candConfNum_ind = find(candConfNum(:)==1);
% [candConfNum_r,candConfNum_c] = ind2sub(size(map_env),candConfNum_ind);
% plot(candConfNum_r,candConfNum_c,'xb','MarkerFaceColor','b');        
% plot(hcs_fusion_r,hcs_fusion_c,'or','MarkerFaceColor','r');
% 
% % fixedConfs_ind = [fixedConfs1_ind,fixedConfs2_ind];
% % [fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),fixedConfs_ind);
% [fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),cell2mat(fixedConfs_ind));
% plot(fixedConfs_r,fixedConfs_c,'or');
% 
% fusionConfs_ind = [fusionConf1,fusionConf2];
% [fusionConfs_r,fusionConfs_c] = ind2sub(size(map_env),fusionConfs_ind);
% plot(fusionConfs_r,fusionConfs_c,'og');
% end
% 
% 
% 
% 
% %----------------------------------------------------
% % -- Check#2: confs. that cover all the hotspots
% %----------------------------------------------------
% conf_hotspot = ones(1,size(oV,2));
% for i = 1:numel(hcs_fusion_ind)
%     
%     % conf ind that cover this hotspot
%     conf_ind = find(oV(hcs_fusion_ind(i),:));
%     % binary vector
%     conf_this = zeros(1,size(oV,2));
%     conf_this(conf_ind) = 1;
%     
%     % - update the overall vector
%     conf_hotspot = conf_hotspot.*conf_this;
% end
% % - make sure fixed conf are in the vector
% % conf_hotspot([fixedConfs1_ind,fixedConfs2_ind]) = 1;
% conf_hotspot(cell2mat(fixedConfs_ind)) = 1;
% 
% % -- plot 
% if visualize == 1
% figure('name','conf cover all the hotspots');
% imshow(map_env','InitialMagnification',800); hold on;
% set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
% conf_hotspot_ind = find(conf_hotspot(:)==1);
% [conf_hotspot_r,conf_hotspot_c] = ind2sub(size(map_env),conf_hotspot_ind);
% plot(conf_hotspot_r,conf_hotspot_c,'xb','MarkerFaceColor','b');        
% plot(hcs_fusion_r,hcs_fusion_c,'or','MarkerFaceColor','r');
% 
% % fixedConfs_ind = [fixedConfs1_ind,fixedConfs2_ind];
% % [fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),fixedConfs_ind);
% [fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),cell2mat(fixedConfs_ind));
% plot(fixedConfs_r,fixedConfs_c,'or');
% 
% fusionConfs_ind = [fusionConf1,fusionConf2];
% [fusionConfs_r,fusionConfs_c] = ind2sub(size(map_env),fusionConfs_ind);
% plot(fusionConfs_r,fusionConfs_c,'og');
% end
% 
% %----------------------------------------------------
% % -- Check#3: conf belong to unoccupied cells
% %----------------------------------------------------
% conf_free = ones(1,size(oV,2));
% for i = 1:numel(OBSconf.ind)
%     %conf_free( (OBS.ind(i)-1)*f+1:(OBS.ind(i)-1)*f+f ) = 0;
%     conf_free( OBSconf.ind(i) ) = 0;
% end
% % - make sure fixed conf are in the vector
% % conf_free([fixedConfs1_ind,fixedConfs2_ind]) = 1;
% conf_free(cell2mat(fixedConfs_ind)) = 1;
% 
% % -- plot 
% if visualize == 1
% figure('name','conf belong to unoccupied cells');
% imshow(map_env','InitialMagnification',800); hold on;
% set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
% conf_free_ind = find(conf_free(:)==1);
% [conf_free_r,conf_free_c] = ind2sub(size(map_env),conf_free_ind);
% plot(conf_free_r,conf_free_c,'xb','MarkerFaceColor','b');
% plot(hcs_fusion_r,hcs_fusion_c,'or','MarkerFaceColor','r');
% 
% % fixedConfs_ind = [fixedConfs1_ind,fixedConfs2_ind];
% % [fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),fixedConfs_ind);
% [fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),cell2mat(fixedConfs_ind));
% plot(fixedConfs_r,fixedConfs_c,'or');
% 
% fusionConfs_ind = [fusionConf1,fusionConf2];
% [fusionConfs_r,fusionConfs_c] = ind2sub(size(map_env),fusionConfs_ind);
% plot(fusionConfs_r,fusionConfs_c,'og');
% end
% 
% % ----------------------------------------------------
% % -- Check#4: conf not close enough to the hotspot 
% % -----------------------------------------------------
% 
% % - cose enough radius
% % innerRadius_cell = para_.innerRadius_m/dataYAML.resolution;
% innerRadius_cell = para_.innerRadius_m/cellsize_env;
% 
% % temporary inner limit map
% map_candConf_minLimit = zeros(size(map_env));
% map_candConf_minLimit(hcs_fusion_ind) = 1;
% 
% se = strel('disk',innerRadius_cell,0);
% map_candConf_minLimit = imdilate(map_candConf_minLimit,se);
% 
% % indices of inner limit map
% innerLimit_ind = find(map_candConf_minLimit);
% 
% % -- fixedconf to be removed from inner limit map %setdiff
% % fixedConfs_ind = [fixedConfs1_ind,fixedConfs2_ind]; % indices of fixed conf
% % fixedConfs_cell = fixedConfs_ind; % cells of fixed conf
% fixedConfs_cell = cell2mat(fixedConfs_ind); % cells of fixed conf
% sacred_ind = ismember(innerLimit_ind,fixedConfs_cell); % conflecting indices 
% innerLimit_ind(sacred_ind) = [];
% 
% % -- plot 
% if visualize == 1
% figure('name','conf inner limit');
% imshow(map_env','InitialMagnification',800); hold on;
% set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
% [innerLimit_r,innerLimit_c] = ind2sub(size(map_env),innerLimit_ind);
% plot(innerLimit_r,innerLimit_c,'xb','MarkerFaceColor','b');
% plot(hcs_fusion_r,hcs_fusion_c,'sr');
% 
% % fixedConfs_ind = [fixedConfs1_ind,fixedConfs2_ind];
% % [fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),fixedConfs_ind);
% [fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),cell2mat(fixedConfs_ind));
% plot(fixedConfs_r,fixedConfs_c,'or');
% 
% fusionConfs_ind = [fusionConf1,fusionConf2];
% [fusionConfs_r,fusionConfs_c] = ind2sub(size(map_env),fusionConfs_ind);
% plot(fusionConfs_r,fusionConfs_c,'og');
% end
% 
% confNotClose = ones(1,size(oV,2));
% confNotClose(innerLimit_ind) = 0;
% confNotClose(cell2mat(fixedConfs_ind)) = 1;
% 
% 
% % -- plot 
% if visualize == 1
% figure('name','conf Not Close');
% imshow(map_env','InitialMagnification',800); hold on;
% set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
% 
% [confNotClose_r,confNotClose_c] = ind2sub(size(map_env),find(confNotClose));
% plot(confNotClose_r,confNotClose_c,'xb','MarkerFaceColor','b');
% plot(hcs_fusion_r,hcs_fusion_c,'sr');
% 
% % fixedConfs_ind = [fixedConfs1_ind,fixedConfs2_ind];
% % [fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),fixedConfs_ind);
% [fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),cell2mat(fixedConfs_ind));
% plot(fixedConfs_r,fixedConfs_c,'or');
% 
% fusionConfs_ind = [fusionConf1,fusionConf2];
% [fusionConfs_r,fusionConfs_c] = ind2sub(size(map_env),fusionConfs_ind);
% plot(fusionConfs_r,fusionConfs_c,'og');
% end
% 
% 
% %-----------------------------------------------------------------------------
% % -- Check#5: confs (position) is not already being used for the other hcs
% %-----------------------------------------------------------------------------
% confNotEngaged = ones(1,size(oV,2));
% confNotEngaged(engagedConfs_ind) = 0;
% confNotEngaged(cell2mat(fixedConfs_ind)) = 1;
% 
% %***************************************************************************
% %
% % Candidate conf are: 
% %       (1) candidate conf map (Check#1) AND 
% %       (2) cover the hotspot (Check#2) AND 
% %       (3) are not placed over occupied cells (Check#3) AND 
% %       (4) are not close enough to the hotspots (Check#4)
% %       (5) are not already selected for the other solutions
% %
% %***************************************************************************
% conf_vector  = candConfNum.*conf_hotspot.*conf_free.*confNotClose.*confNotEngaged;
% confKept     = find(conf_vector);
% confRemoved  = find(conf_vector==0);
% % V            = V(:,confKept);
% % V            = full(V);
% 
% 
% % -- plot 
% if visualize == 1
% figure('name','conf meeting all four conditions');
% imshow(map_env','InitialMagnification',800); hold on;
% set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
% [confKept_r,confKept_c] = ind2sub(size(map_env),confKept);
% plot(confKept_r,confKept_c,'xb','MarkerFaceColor','b');
% plot(hcs_fusion_r,hcs_fusion_c,'or','MarkerFaceColor','r');
% 
% % fixedConfs_ind = [fixedConfs1_ind,fixedConfs2_ind];
% % [fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),fixedConfs_ind);
% [fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),cell2mat(fixedConfs_ind));
% plot(fixedConfs_r,fixedConfs_c,'or');
% 
% fusionConfs_ind = [fusionConf1,fusionConf2];
% [fusionConfs_r,fusionConfs_c] = ind2sub(size(map_env),fusionConfs_ind);
% plot(fusionConfs_r,fusionConfs_c,'og');
% end
% 
% %***********************************************************************
% % 
% %    Limit the number of cand conf based on Proportional Coverage
% % 
% %***********************************************************************
% 
% %/////////////////////////////////////// 
% % -- Concentration weights
% %///////////////////////////////////////
% 
% cellKeptWeights = map_recon(cellKept);
% % cellKeptWeights = cellKeptWeights/(max(cellKeptWeights));
% 
% V_this = V(:,confKept);
% V_this = V_this .* repmat(cellKeptWeights',1,numel(confKept));
% wV_concentration_this = sum(V_this,1);
% wV_concentration_this = wV_concentration_this/max(wV_concentration_this(:));
% 
% 
% % -- make up for fixedconf
% % fixedConfs_ind = [fixedConfs1_ind,fixedConfs2_ind];
% fixed_ind = ismember(confKept,cell2mat(fixedConfs_ind));
% wV_concentration_this(fixed_ind) = max(wV_concentration_this(:)+0.5);
% 
% 
% %/////////////////////////////////////// 
% % -- Placement weights
% %/////////////////////////////////////// 
% 
% % --------------------
% % == Orientation
% % --------------------
% % -- conf orientations
% ConfOrns = oConfOrns(confKept);
% 
% % -- desired orinetations w.r.t to hotspot
% % [hotspot_r,hotspot_c] = ind2sub(size(map_env),hotspot);
% [conf_r,conf_c] = ind2sub(size(map_candConf),confKept);
% 
% % desiredOrns1 = round(rad2deg(angle2Points([conf_r',conf_c'],[hotspots_r(1),hotspots_c(1)])));
% % desiredOrns1 = desiredOrns1 + (desiredOrns1==0)*360;
% % 
% % desiredOrns2 = round(rad2deg(angle2Points([conf_r',conf_c'],[hotspots_r(2),hotspots_c(2)])));
% % desiredOrns2 = desiredOrns2 + (desiredOrns2==0)*360;
% 
% desiredOrns = cell(1,numel(hcs_fusion_ind));
% for i = 1:numel(hcs_fusion_ind)
%     desiredOrns{i} = round(rad2deg(angle2Points([conf_r',conf_c'],[hcs_fusion_r(i),hcs_fusion_c(i)])));
%     desiredOrns{i} = desiredOrns{i} + (desiredOrns{i}==0)*360;    
% end
% 
% % -- delta orientations
% % ConfDeltaOrns1 = rad2deg(angleAbsDiff(deg2rad(desiredOrns1(:)),deg2rad(ConfOrns(:))));
% % ConfDeltaOrns2 = rad2deg(angleAbsDiff(deg2rad(desiredOrns2(:)),deg2rad(ConfOrns(:))));
% 
% ConfDeltaOrns = cell(1,numel(hcs_fusion_ind));
% for i = 1:numel(hcs_fusion_ind)
%     ConfDeltaOrns{i} = rad2deg(angleAbsDiff(deg2rad(desiredOrns{i}(:)),deg2rad(ConfOrns(:))));
% end
% 
% % -- angular gain
% gain_sigma = FoV/6;%15; % standard deviation
% gain_meu   = 0; % expected value or center
% % V_angular_gain1_this = gaussmf(ConfDeltaOrns1,[gain_sigma,gain_meu]);
% % V_angular_gain2_this = gaussmf(ConfDeltaOrns2,[gain_sigma,gain_meu]);
% 
% V_angular_gain_this = cell(1,numel(hcs_fusion_ind));
% for i = 1:numel(hcs_fusion_ind)
%     V_angular_gain_this{i} = gaussmf(ConfDeltaOrns{i},[gain_sigma,gain_meu]);
% end
% 
% 
% % ------------------------------
% % == displacement
% % ------------------------------
% 
% % -- dist from conf to hotspot
% % conf2Hot_dist1 = pdist2([conf_r',conf_c'],[hotspots_r(1),hotspots_c(1)],'euclidean');
% % conf2Hot_dist2 = pdist2([conf_r',conf_c'],[hotspots_r(2),hotspots_c(2)],'euclidean');
% 
% conf2Hot_dist = cell(1,numel(hcs_fusion_ind));
% for i = 1:numel(hcs_fusion_ind)
%     conf2Hot_dist{i} = pdist2([conf_r',conf_c'],[hcs_fusion_r(i),hcs_fusion_c(i)],'euclidean');
% end
% 
% % -- distance gain
% gain_sigma       = SensorRange.cell/5; % standard deviation
% gain_meu         = SensorRange.cell/2; % expected value or center
% % V_dist_gain1_this = gaussmf(conf2Hot_dist1,[gain_sigma,gain_meu]);
% % V_dist_gain2_this = gaussmf(conf2Hot_dist2,[gain_sigma,gain_meu]);
% 
% V_dist_gain_this = cell(1,numel(hcs_fusion_ind));
% for i = 1:numel(hcs_fusion_ind)
%     V_dist_gain_this{i} = gaussmf(conf2Hot_dist{i},[gain_sigma,gain_meu]);
% end
% 
% 
% 
% % ------------------------------
% % == orientation + displacement
% % ------------------------------
% % wV_diplacement1_this = (V_angular_gain1_this+V_dist_gain1_this)/2;
% % wV_diplacement2_this = (V_angular_gain2_this+V_dist_gain2_this)/2;
% 
% wV_diplacement_this = cell(1,numel(hcs_fusion_ind));
% for i = 1:numel(hcs_fusion_ind)
%     wV_diplacement_this{i} = (V_angular_gain_this{i}+V_dist_gain_this{i})/2;
% end
% 
% % wV1_this = (gamma*(wV_concentration_this))' + ((1-gamma)* wV_diplacement1_this);
% % wV2_this = (gamma*(wV_concentration_this))' + ((1-gamma)* wV_diplacement2_this);
% % 
% % wV_this = (wV1_this + wV2_this)/2;
% 
% wV_this = cell(1,numel(hcs_fusion_ind));
% for i = 1:numel(hcs_fusion_ind)
%     wV_this{i} = (gamma*(wV_concentration_this))' + ((1-gamma)* wV_diplacement_this{i});
% end
% 
% wV_all = zeros(size(wV_this{1}));
% for i = 1:numel(hcs_fusion_ind)
%     wV_all = wV_all+wV_this{i};
% end
% wV_all = wV_all/numel(hcs_fusion_ind);
% 
% % -- sorted indices
% % [~,sorted_ind] = sort(wV_this,'descend');
% [~,sorted_ind] = sort(wV_all,'descend');
% 
% sorted_ind_max = sorted_ind(1:min([para_.maxCandidateConfFusion2SE,...
%     numel(sorted_ind)]));
% % sorted_ind_max = sorted_ind(1:min([10,numel(sorted_ind)]));
% confToRemove_ind = find(ismember(sorted_ind,sorted_ind_max)==0);
% confRemoved = union(confRemoved,sorted_ind(confToRemove_ind));
% confKept = confKept(sorted_ind_max);
% 
% %--- hardcore fixed confs
% confKept = union(confKept,cell2mat(fixedConfs_ind));
% 
% % -- sorted list
% confKept = sort(confKept);
% 
% % -- plot 
% if visualize == 1
% figure('name','max num of conf');
% imshow(map_env','InitialMagnification',800); hold on;
% set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
% [confKept_r,confKept_c] = ind2sub(size(map_env),confKept);
% plot(confKept_r,confKept_c,'xb','MarkerFaceColor','b');
% plot(hcs_fusion_r,hcs_fusion_c,'sr');
% 
% % fixedConfs_ind = [fixedConfs1_ind,fixedConfs2_ind];
% % [fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),fixedConfs_ind);
% [fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),cell2mat(fixedConfs_ind));
% plot(fixedConfs_r,fixedConfs_c,'or');
% 
% fusionConfs_ind = [fusionConf1,fusionConf2];
% [fusionConfs_r,fusionConfs_c] = ind2sub(size(map_env),fusionConfs_ind);
% plot(fusionConfs_r,fusionConfs_c,'og');
% end
% 
% 
% % ****************************************************************
% %     conf kept visibility matrix
% % ****************************************************************
% 
% V = V(:,confKept);
% V = full(V);
% 
% % size(V)
% 
% % ****************************************************************
% %     weighted visibility matrix
% % ****************************************************************
% 
% 
% % wV = V.*repmat(cellKeptWeights',1,numel(confKept));
% % wV = V.*repmat(map_recon(cellKept)',1,numel(confKept));
% % max_conc_conf = max(sum(wV,1));
% % wV = wV./max_conc_conf;
% WV = V.*repmat(map_recon(cellKept)',1,numel(confKept));
% max_conc_conf = max(sum(WV,1));
% WV = WV./max_conc_conf;
% 
% % -- weighted visibility for fixedCo  nf1 and fixedConf2
% % fixedConf1_num = find(ismember(confKept,fixedConfs1_ind));
% % fixedConf2_num = find(ismember(confKept,fixedConfs2_ind));
% 
% fixedConfs_num = cell(1,numel(hcs_fusion_ind));
% for i = 1:numel(hcs_fusion_ind)
%     fixedConfs_num{i} = find(ismember(confKept,fixedConfs_ind{i}));
% end
% 
% 
% % fixed confs of other local solutions (hcs_not_this) that can not cover
% % this hc, can not be used for the visibility, cross angles, and distance
% % gains.
% % 
% % fixedConfXG_num:= fixed conf num that must not be used for the gain
% % calculations of this hc.
% fixedConfXG_num = cell(1,numel(hcs_fusion_ind));
% wV = cell(1,numel(hcs_fusion_ind));
% num_all = 1:numel(hcs_fusion_ind);
% for i = 1:numel(hcs_fusion_ind)
%     wV{i} = WV;
%     fixedConfOtherThanThisHc_num = cell2mat(fixedConfs_num(num_all~=i));
%     CellNumThisHc = find(cellKept==hcs_fusion_ind(i));
%     %this_ind = wV{i}(cellKept==hc_fusion_ind(i),cell2mat(fixedConfs_num(num_all~=i)))==0
%     this_ind = wV{i}(CellNumThisHc,fixedConfOtherThanThisHc_num)==0;
%     fixedConfXG_num{i} = fixedConfOtherThanThisHc_num(this_ind);
%     wV{i}(:,fixedConfXG_num{i}) = 0;
% end
%     
% 
% %{
% wV1 = wV;
% % wV1(:,fixedConf2_num) = 0;
% % hotspotCellKept = find(cellKept==hotspots(1))
% % fixedConf2XG_num = fixedConf2_num(find(wV1(hotspotCellKept,fixedConf2_num)==0))
% fixedConf2XG_num = fixedConf2_num(wV1(cellKept==hotspots(1),fixedConf2_num)==0);
% wV1(:,fixedConf2XG_num) = 0;
% 
% 
% wV2 = wV;
% % wV2(:,fixedConf1_num) = 0;
% fixedConf1XG_num = fixedConf1_num(wV2(cellKept==hotspots(2),fixedConf1_num)==0);
% wV2(:,fixedConf1XG_num) = 0;
% %}
% 
% % -- fixed confs that do not cover the 2nd hotspot in fusion
% %fixedConf1XG_num 
% %fixedConf2XG_num
% 
% %fixedConf1_num
% %fixedConf2_num
% 
% 
% 
% %//////////////////////////////////////////
% % HOTSPOT PLACEMENT GAIN - ANGULAR
% %//////////////////////////////////////////
% 
% % -- conf orientations
% ConfOrns = oConfOrns(confKept);
% 
% 
% % -- desired orinetations w.r.t to hotspot
% % [hotspots_r,hotspots_c] = ind2sub(size(map_env),hotspots);
% [conf_r,conf_c] = ind2sub(size(map_candConf),confKept);
% 
% % desiredOrns1 = round(rad2deg(angle2Points([conf_r',conf_c'],[hotspots_r(1),hotspots_c(1)])));
% % desiredOrns1 = desiredOrns1 + (desiredOrns1==0)*360;
% % 
% % desiredOrns2 = round(rad2deg(angle2Points([conf_r',conf_c'],[hotspots_r(2),hotspots_c(2)])));
% % desiredOrns2 = desiredOrns2 + (desiredOrns2==0)*360;
% 
% desiredOrns = cell(1,numel(hcs_fusion_ind));
% for i = 1:numel(hcs_fusion_ind)
%     desiredOrns{i} = round(rad2deg(angle2Points([conf_r',conf_c'],[hcs_fusion_r(i),hcs_fusion_c(i)])));
%     desiredOrns{i} = desiredOrns{i} + (desiredOrns{i}==0)*360;
% end
% 
% % -- delta orientations
% % ConfDeltaOrns1 = rad2deg(angleAbsDiff(deg2rad(desiredOrns1(:)),deg2rad(ConfOrns(:))));
% % ConfDeltaOrns2 = rad2deg(angleAbsDiff(deg2rad(desiredOrns2(:)),deg2rad(ConfOrns(:))));
% 
% ConfDeltaOrns = cell(1,numel(hcs_fusion_ind));
% for i = 1:numel(hcs_fusion_ind)
%     ConfDeltaOrns{i} = rad2deg(angleAbsDiff(deg2rad(desiredOrns{i}(:)),deg2rad(ConfOrns(:))));
% end
%    
% % -- angular gain
% gain_sigma = FoV/6;%15; % standard deviation
% gain_meu   = 0; % expected value or center
% % V_angular_gain1 = gaussmf(ConfDeltaOrns1,[gain_sigma,gain_meu]);
% % V_angular_gain2 = gaussmf(ConfDeltaOrns2,[gain_sigma,gain_meu]);
% 
% V_angular_gain = cell(1,numel(hcs_fusion_ind));
% for i = 1:numel(hcs_fusion_ind)
%     V_angular_gain{i} = gaussmf(ConfDeltaOrns{i},[gain_sigma,gain_meu]);
% end
% 
% 
% %//////////////////////////////////////////
% % HOTSPOT PLACEMENT GAIN - DISTANCE
% %//////////////////////////////////////////
% 
% % -- dist from conf to hotspot
% % conf2Hot_dist1 = pdist2([conf_r',conf_c'],[hotspots_r(1),hotspots_c(1)],'euclidean');
% % conf2Hot_dist2 = pdist2([conf_r',conf_c'],[hotspots_r(2),hotspots_c(2)],'euclidean');
% 
% conf2Hot_dist = cell(1,numel(hcs_fusion_ind));
% for i = 1:numel(hcs_fusion_ind)
%     conf2Hot_dist{i} = ...
%         pdist2([conf_r',conf_c'],[hcs_fusion_r(i),hcs_fusion_c(i)],'euclidean');
% end
% 
% 
% 
% % -- distance gain
% gain_sigma = SensorRange.cell/5; % standard deviation
% gain_meu   = SensorRange.cell/2; % expected value or center
% % V_dist_gain1 = gaussmf(conf2Hot_dist1,[gain_sigma,gain_meu]);
% % V_dist_gain2 = gaussmf(conf2Hot_dist2,[gain_sigma,gain_meu]);
% 
% V_dist_gain = cell(1,numel(hcs_fusion_ind));
% for i = 1:numel(hcs_fusion_ind)
%     V_dist_gain{i} = gaussmf(conf2Hot_dist{i},[gain_sigma,gain_meu]);
% end
% 
% 
% %//////////////////////////////////////////////
% % HOTSPOT PLACEMENT GAIN - ANGULAR + DISTANCE
% %//////////////////////////////////////////////
% % Conf_DistOrnGainV1 = (V_angular_gain1+V_dist_gain1)/2;
% % Conf_DistOrnGainV2 = (V_angular_gain2+V_dist_gain2)/2;
% 
% Conf_DistOrnGainV = cell(1,numel(hcs_fusion_ind));
% for i = 1:numel(hcs_fusion_ind)
%     Conf_DistOrnGainV{i} = (V_angular_gain{i}+V_dist_gain{i})/2;
% end
% 
% 
% 
% 
% end
% 

%--------------------------------- End of the document -------------------------------


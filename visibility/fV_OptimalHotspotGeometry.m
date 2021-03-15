function [ V,...
           oV,...
           wV,...
           Conf_DistOrnGainV,...
           confKept,...
           confRemoved,...
           cellKept,...
           cellRemoved,...
           oConfOrns,...
           ConfOrns,...
           executedConfsHot_num,...
           FoVfaces ] = fV_OptimalHotspotGeometry( map_candConf,...
                                                   map_coverage,...
                                                   map_env,...
                                                   map_recon,...
                                                   engagedConfs_ind,...
                                                   executedConfsHot_ind,...
                                                   executedConfsHot_orn,...
                                                   n,...
                                                   FoV,...
                                                   FoVfaces,...
                                                   SensorRange,...
                                                   OBSenv,...
                                                   OBSconf,...
                                                   hc_this_sub,...
                                                   cellsize_env,...
                                                   para_,...
                                                   visualize )
%fV generates Visibility matrix V.

% --- ANGLES:
FoVfaces.lgt = [1/FoVfaces.num:1/FoVfaces.num:1-(1/FoVfaces.num) 0]'*360;
FoVfaces.alf = FoVfaces.lgt-(FoV/2);
FoVfaces.bay = FoVfaces.lgt+(FoV/2);
% compensate negative angles
%FoVfaces.alf(FoVfaces.alf<0) = 360+FoVfaces.alf(FoVfaces.alf<0);
FoVfaces.alf = wrapTo360(FoVfaces.alf);
FoVfaces.bay = wrapTo360(FoVfaces.bay);
FoVfaces.ang = [FoVfaces.alf FoVfaces.bay];
FoVfaces.ang = [FoVfaces.alf FoVfaces.bay];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                V
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disp('----- Visibility Matrix');
% V = sparse(n,c); % initializing
V = sparse(n,n); % initializing
% V = zeros(n,c); % initializing
% [hotspot_r,hotspot_c] = ind2sub(size(map_env),hc_this);
hc_r = hc_this_sub(1);
hc_c = hc_this_sub(2);
hc_i = sub2ind(size(map_env),hc_r,hc_c);

oConfOrns = zeros(1,n);
% For each cell of the map.
for i = 1:n
    if map_candConf(i) % if its not occupied
        [rr,cc] = ind2sub(size(map_candConf),i);
        sensor = [rr,cc];
        % --
        orientation = round(rad2deg(angle2Points(sensor,[hc_r,hc_c])));
        orientation = orientation + (orientation==0)*360;
        
        % For each configuration of a cell.
        for j = orientation %1%:f
            
            sD = []; oD = []; VObs = [];
            %this_dir = sprintf(['sD_oD_VObs/sensing_range_%02dcell/fov_%03ddeg/',...
            %                    'orientation_%03ddeg/'],SensorRange.cell,FoV,j);
            this_dir = sprintf('visibility/templates/sensing_range_%02dcell/fov_%03ddeg/',...
                SensorRange.cell,FoV);
            this_file = sprintf('sD_oD_VObs_rng%02d_fov%03d_orn%03d.mat',...
                        SensorRange.cell,FoV,j);
            load([this_dir,this_file],'sD','oD','VObs');
            
            % sensorDom := sensor domain translation w.r.t sensor position.
            sensorDom = [sD(:,1)+sensor(1) sD(:,2)+sensor(2)];
            
            % obstacDom := obstalce domain translated w.r.t sensor's
            % position. 
            obstacDom = [oD(:,1)+sensor(1) oD(:,2)+sensor(2)];
            
            % VOBS := Visibility vs Obstacles matrix for jth FoV.
            VOBS = VObs;
            
            % visibile cells 
            [ visDom ] = fVisibileDomain( map_env,sensorDom,obstacDom,VOBS );
            
            % Update Visibility matrix.
            %V(visDom,((i-1)*f)+j) = 1;
            V(visDom,i) = 1;
        end
        oConfOrns(i) = orientation;
    end
end


%-------------------------
% V for Executed Confs
%-------------------------
for i = 1:numel(executedConfsHot_ind)
    
    [rr,cc] = ind2sub(size(map_candConf),executedConfsHot_ind(i));
    sensor  = [rr,cc];
        
    j = executedConfsHot_orn(i) + (executedConfsHot_orn(i)==0)*360;
    
    sD = []; oD = []; VObs = [];
    %this_dir = sprintf(['sD_oD_VObs/sensing_range_%02dcell/fov_%03ddeg/',...
    %                    'orientation_%03ddeg/'],SensorRange.cell,FoV,j);
    this_dir = sprintf('visibility/templates/sensing_range_%02dcell/fov_%03ddeg/',...
                SensorRange.cell,FoV);
    this_file = sprintf('sD_oD_VObs_rng%02d_fov%03d_orn%03d.mat',...
                SensorRange.cell,FoV,j);
    load([this_dir,this_file],'sD','oD','VObs');
    
    % sensorDom := sensor domain translation w.r.t sensor position.
    sensorDom = [sD(:,1)+sensor(1) sD(:,2)+sensor(2)];

    % obstacDom := obstalce domain translated w.r.t sensor's
    % position.
    obstacDom = [oD(:,1)+sensor(1) oD(:,2)+sensor(2)];

    % VOBS := Visibility vs Obstacles matrix for jth FoV.
    VOBS = VObs;
        
    % visibile cells 
    [ visDom ] = fVisibileDomain( map_env,sensorDom,obstacDom,VOBS );

    % Update Visibility matrix.
    %V(visDom,((i-1)*f)+j) = 1;
    V(visDom,executedConfsHot_ind(i)) = 1;
    
    oConfOrns(executedConfsHot_ind(i)) = j;
end
   

%-------------------------------------------------
% sensing config can not visualize own cell
%-------------------------------------------------
for i = 1:numel(map_candConf)    
    %V(i,(i-1)*f+1:(i-1)*f+f) = 0;
    V(i,i) = 0;
end

oV = V; % oV := visibility with obstacles

% V = sparse(V);

% VISIBILITY MATRIX CORRESPOND TO THE HOTSPOT ONLY:
% ------------------------------------------------------------------------------

%////////////////////////////////////////
%
% CELLS TO RETAIN
%
%////////////////////////////////////////


%----------------------------------------------------
% -- cells covered by candidate conf.
%----------------------------------------------------
covered_ind = sum(oV,2);
covered_num = find(covered_ind);
covered_cell_vector = zeros(1,size(oV,1));
covered_cell_vector(covered_num) = 1;

%----------------------------------------------------
% -- cells need to be covered -- according to the coverage map
%----------------------------------------------------
coverage_ind = find(map_coverage);
coverage_cell_vector = zeros(1,size(oV,1));
% size(coverage_cell_vector)
coverage_cell_vector(coverage_ind) = 1;
% size(coverage_cell_vector)

%----------------------------------------------------
% -- unoccupied cells
%----------------------------------------------------
free_cell_vector = ones(1,size(oV,1));
free_cell_vector(OBSenv.ind) = 0;


% --- selecting cells that are covered and need to be covered and unoccupied.
%***************************************************************************
%
% cells kept to observe are: 
%       (1) can be covered by cand conf (Check#1) AND 
%       (2) desired to be covered (Check#2) AND 
%       (3) are not occupied (Check#3) 
%
%***************************************************************************
cov_cell     = covered_cell_vector.*coverage_cell_vector.*free_cell_vector;
cellKept     = find(cov_cell>=1);
cellRemoved  = find(cov_cell==0);
V            = V(cellKept,:);

% ////////////////////////////
% -- cellKept weights
% ////////////////////////////
% cellKeptWeights = map_recon(cellKept);
% cellKeptWeights = cellKeptWeights/(max(cellKeptWeights));


% //////////////////////////
% ---- configurations:
% //////////////////////////

%----------------------------------------------------
% -- Check#1: conf corresponding to candidate conf map only
%----------------------------------------------------
% cells for the candidate conf.
candConfCells = find(map_candConf); 
% conf numbers
candConfNum = zeros(1,size(oV,2));
for i = 1:numel(candConfCells)
    %candConfNum( (candConfCells(i)-1)*f+1:(candConfCells(i)-1)*f+f ) = 1;
    candConfNum( candConfCells(i) ) = 1;
end
% - fixed confs are in the vector
candConfNum(executedConfsHot_ind) = 1;

% -- plot 
if visualize == 1
figure('name','conf corresponding to candidate conf map only');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
candConfNum_ind = find(candConfNum(:)==1);
[candConfNum_r,candConfNum_c] = ind2sub(size(map_env),candConfNum_ind);
plot(candConfNum_r,candConfNum_c,'xb','MarkerFaceColor','b');        
plot(hc_r,hc_c,'or','MarkerFaceColor','r');

fixedConfs_ind = executedConfsHot_ind;
[fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),fixedConfs_ind);
plot(fixedConfs_r,fixedConfs_c,'or');
end

%----------------------------------------------------
% -- Check#2: confs. that cover all the hotspots
%----------------------------------------------------
conf_hc = zeros(1,size(oV,2));

% conf ind that cover this hotspot
conf_ind = find(oV(hc_i,:));    
conf_hc(conf_ind) = 1;

% - make sure fixed conf are in the vector
conf_hc(executedConfsHot_ind) = 1;

% -- plot 
if visualize == 1
figure('name','conf cover all the hotspots');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
conf_hotspot_ind = find(conf_hc(:)==1);
[conf_hotspot_r,conf_hotspot_c] = ind2sub(size(map_env),conf_hotspot_ind);
plot(conf_hotspot_r,conf_hotspot_c,'xb','MarkerFaceColor','b');        
plot(hc_r,hc_c,'or','MarkerFaceColor','r');

fixedConfs_ind = executedConfsHot_ind;
[fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),fixedConfs_ind);
plot(fixedConfs_r,fixedConfs_c,'or');

end

%----------------------------------------------------
% -- Check#3: conf belong to unoccupied cells
%----------------------------------------------------
conf_free = ones(1,size(oV,2));
for i = 1:numel(OBSconf.ind)
    %conf_free( (OBS.ind(i)-1)*f+1:(OBS.ind(i)-1)*f+f ) = 0;
    conf_free( OBSconf.ind(i) ) = 0;
end
% - make sure fixed conf are in the vector
conf_free(executedConfsHot_ind) = 1;

% -- plot 
if visualize == 1
figure('name','conf belong to unoccupied cells');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
conf_free_ind = find(conf_free(:)==1);
[conf_free_r,conf_free_c] = ind2sub(size(map_env),conf_free_ind);
plot(conf_free_r,conf_free_c,'xb','MarkerFaceColor','b');
plot(hc_r,hc_c,'or','MarkerFaceColor','r');

fixedConfs_ind = executedConfsHot_ind;
[fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),fixedConfs_ind);
plot(fixedConfs_r,fixedConfs_c,'or');

end


% ----------------------------------------------------
% -- Check#4: conf not close enough to the hotspot 
% -----------------------------------------------------

% - cose enough radius
% innerRadius_cell = para_.innerRadius_m/dataYAML.resolution;
innerRadius_cell = para_.innerRadius_m/cellsize_env;

% temporary inner limit map
map_candConf_minLimit = zeros(size(map_env));
map_candConf_minLimit(hc_i) = 1;

se = strel('disk',innerRadius_cell,0);
map_candConf_minLimit = imdilate(map_candConf_minLimit,se);

% indices of inner limit map
innerLimit_ind = find(map_candConf_minLimit);


% -- fixedconf to be removed from inner limit map %setdiff
fixedConfs_ind = executedConfsHot_ind; % indices of fixed conf
fixedConfs_cell = fixedConfs_ind; % cells of fixed conf
sacred_ind = ismember(innerLimit_ind,fixedConfs_cell); % conflecting indices 
innerLimit_ind(sacred_ind) = [];

% -- plot 
if visualize == 1
figure('name','conf inner limit');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
[innerLimit_r,innerLimit_c] = ind2sub(size(map_env),innerLimit_ind);
plot(innerLimit_r,innerLimit_c,'xb','MarkerFaceColor','b');
plot(hc_r,hc_c,'sr');

fixedConfs_ind = executedConfsHot_ind;
[fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),fixedConfs_ind);
plot(fixedConfs_r,fixedConfs_c,'or');

end

confNotClose = ones(1,size(oV,2));
confNotClose(innerLimit_ind) = 0;
confNotClose(executedConfsHot_ind) = 1;



% -- plot 
if visualize == 1
figure('name','conf Not Close');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])

[confNotClose_r,confNotClose_c] = ind2sub(size(map_env),find(confNotClose));
plot(confNotClose_r,confNotClose_c,'xb','MarkerFaceColor','b');
plot(hc_r,hc_c,'sr');

fixedConfs_ind = executedConfsHot_ind;
[fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),fixedConfs_ind);
plot(fixedConfs_r,fixedConfs_c,'or');

end


%-----------------------------------------------------------------------------
% -- Check#5: confs (position) is not already being used for the other hcs
%-----------------------------------------------------------------------------
confNotEngaged = ones(1,size(oV,2));
confNotEngaged(engagedConfs_ind) = 0;
confNotEngaged(executedConfsHot_ind) = 1;


%***************************************************************************
%
% Candidate conf are: 
%       (1) candidate conf map (Check#1) AND 
%       (2) cover the hotspot (Check#2) AND 
%       (3) are not placed over occupied cells (Check#3) AND 
%       (4) are not close enough to the hotspots (Check#4)
%       (5) are not already selected for the other solutions
%
%***************************************************************************
conf_vector  = candConfNum.*conf_hc.*conf_free.*confNotClose.*confNotEngaged;
confKept     = find(conf_vector);
confRemoved  = find(conf_vector==0);
% V            = V(:,confKept);
% V            = full(V);

% -- plot 
if visualize == 1
figure('name','conf meeting all four conditions');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
[confKept_r,confKept_c] = ind2sub(size(map_env),confKept);
plot(confKept_r,confKept_c,'xb','MarkerFaceColor','b');
plot(hc_r,hc_c,'or','MarkerFaceColor','r');

fixedConfs_ind = executedConfsHot_ind;
[fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),fixedConfs_ind);
plot(fixedConfs_r,fixedConfs_c,'or');
end


%***********************************************************************
% 
%    Limit the number of cand conf based on Proportional Coverage
% 
%***********************************************************************

%/////////////////////////////////////// 
% -- Concentration weights
%///////////////////////////////////////

cellKeptWeights = map_recon(cellKept);
%cellKeptWeights = cellKeptWeights/(max(cellKeptWeights));

V_this = V(:,confKept);
V_this = V_this .* repmat(cellKeptWeights',1,numel(confKept));
wV_concentration_this = sum(V_this,1);
wV_concentration_this = wV_concentration_this/max(wV_concentration_this(:));


% -- make up for fixedconf
fixedConfs_ind = executedConfsHot_ind;
fixed_ind = ismember(confKept,fixedConfs_ind);
wV_concentration_this(fixed_ind) = max(wV_concentration_this(:)+0.5);


%/////////////////////////////////////// 
% -- Placement weights
%/////////////////////////////////////// 

% --------------------
% == Orientation
% --------------------
% -- conf orientations
ConfOrns = oConfOrns(confKept);

% -- desired orinetations w.r.t to hotspot
% [hc_r,hc_c] = ind2sub(size(map_env),hc_this);
[conf_r,conf_c] = ind2sub(size(map_candConf),confKept);

desiredOrns = round(rad2deg(angle2Points([conf_r',conf_c'],[hc_r,hc_c])));
desiredOrns = desiredOrns + (desiredOrns==0)*360;

% -- delta orientations
ConfDeltaOrns = rad2deg(angleAbsDiff(deg2rad(desiredOrns(:)),deg2rad(ConfOrns(:))));

% -- angular gain
gain_sigma = FoV/6;%15; % standard deviation
gain_meu   = 0; % expected value or center
V_angular_gain_this = gaussmf(ConfDeltaOrns,[gain_sigma,gain_meu]);

% ------------------------------
% == displacement
% ------------------------------

% -- dist from conf to hotspot
conf2Hot_dist = pdist2([conf_r',conf_c'],[hc_r,hc_c],'euclidean');

% -- distance gain
gain_sigma       = SensorRange.cell/5; % standard deviation
gain_meu         = SensorRange.cell/2; % expected value or center
V_dist_gain_this = gaussmf(conf2Hot_dist,[gain_sigma,gain_meu]);


% ------------------------------
% == orientation + displacement
% ------------------------------
wV_diplacement_this = (V_angular_gain_this+V_dist_gain_this)/2;

%wV_this_a = ((gamma*(wV_concentration_this)))';
%wV_this_b = ((1-gamma)* wV_diplacement_this);
wV_this_a = ((para_.gamma*(wV_concentration_this)))';
wV_this_b = ((1-para_.gamma)* wV_diplacement_this);

wV_this = wV_this_a+wV_this_b;
%wV_this = (gamma*(wV_concentration_this))' + ((1-gamma)* wV_diplacement_this)


% -- sorted indices
[~,sorted_ind] = sort(wV_this,'descend');
% sorted_ind_max = sorted_ind(1:min([para_.NumOfCandidateConfTomography1SE,...
%     numel(sorted_ind)]));
%sorted_ind_max = sorted_ind(1:min([para_.maxCandidateConfLocal1SE,numel(sorted_ind)]));
switch para_.MissionStrategy
    case '1t-armex'
        para_.maxCandidateConf_Hotspot = para_.maxCandidateConf_Hotspot_1t;
    case '2t-armex'
        para_.maxCandidateConf_Hotspot = para_.maxCandidateConf_Hotspot_2t;
end
sorted_ind_max = sorted_ind(1:min([para_.maxCandidateConf_Hotspot,numel(sorted_ind)]));

% sorted_ind_max = sorted_ind(1:min([10,numel(sorted_ind)]));
confToRemove_ind = find(ismember(sorted_ind,sorted_ind_max)==0);
confRemoved = union(confRemoved,sorted_ind(confToRemove_ind));
confKept = confKept(sorted_ind_max);


%-- make sure that the executed confs are included
confKept = union(confKept,executedConfsHot_ind);

% -- sorted list
confKept = sort(confKept);


%FIXME TODO --- temporary fix column/row vector, debug later
confKept = confKept(:);
confKept = confKept';


% -- plot 
if visualize == 1
figure('name','max num of conf');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
[confKept_r,confKept_c] = ind2sub(size(map_env),confKept);
plot(confKept_r,confKept_c,'xb','MarkerFaceColor','b');
plot(hc_r,hc_c,'sr');

fixedConfs_ind = executedConfsHot_ind;
[fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),fixedConfs_ind);
plot(fixedConfs_r,fixedConfs_c,'or');

end

% ****************************************************************
%     conf kept visibility matrix
% ****************************************************************

V = V(:,confKept);
V = full(V);

% size(V)

% ****************************************************************
%     weighted visibility matrix
% ****************************************************************

% wV = V.*repmat(cellKeptWeights',1,numel(confKept));
wV = V.*repmat(map_recon(cellKept)',1,numel(confKept));
max_conc_conf = max(sum(wV,1));
wV = wV./max_conc_conf;

% sum(map_recon(cellKept))
% wV = wV/sum(map_recon(cellKept));
% pause

% -- weighted visibility for fixedConf1 and fixedConf2
executedConfsHot_num = find(ismember(confKept,executedConfsHot_ind));

%//////////////////////////////////////////
% HOTSPOT PLACEMENT GAIN - ANGULAR
%//////////////////////////////////////////

% -- conf orientations
ConfOrns = oConfOrns(confKept);


% -- desired orinetations w.r.t to hotspot
% [hc_r,hc_c] = ind2sub(size(map_env),hc_this);
[conf_r,conf_c] = ind2sub(size(map_candConf),confKept);

desiredOrns = round(rad2deg(angle2Points([conf_r',conf_c'],[hc_r,hc_c])));
desiredOrns = desiredOrns + (desiredOrns==0)*360;

% -- delta orientations
ConfDeltaOrns = rad2deg(angleAbsDiff(deg2rad(desiredOrns(:)),deg2rad(ConfOrns(:))));

% -- angular gain
gain_sigma = FoV/6;%15; % standard deviation
gain_meu   = 0; % expected value or center
V_angular_gain = gaussmf(ConfDeltaOrns,[gain_sigma,gain_meu]);

%//////////////////////////////////////////
% HOTSPOT PLACEMENT GAIN - DISTANCE
%//////////////////////////////////////////

% -- dist from conf to hotspot
conf2Hot_dist = pdist2([conf_r',conf_c'],[hc_r,hc_c],'euclidean');

% -- distance gain
gain_sigma = SensorRange.cell/5; % standard deviation
gain_meu   = SensorRange.cell/2; % expected value or center
V_dist_gain = gaussmf(conf2Hot_dist,[gain_sigma,gain_meu]);


%//////////////////////////////////////////////
% HOTSPOT PLACEMENT GAIN - ANGULAR + DISTANCE
%//////////////////////////////////////////////
Conf_DistOrnGainV = (V_angular_gain+V_dist_gain)/2;


% ------ publish gussian gain function for test only
%{ 
angles = 0:1:15; % in degrees
gain_sigma = 3; % standard deviation
gain_meu   = 7.5; % expected value or center
guss_gain = gaussmf(angles,[gain_sigma,gain_meu]);
figure; plot(angles,guss_gain ,'-b')
%}        

% sum(wV,1)
% 
% size(Conf_DistOrnGainV)
% % 
% % one = (0.9*wV)
% % 
% % two = 0.1*repmat(Conf_DistOrnGainV',numel(cellKept),1);
% 
% sum(0.9*wV,1)
% 0.1*(Conf_DistOrnGainV')
% 
% wV1 = sum(0.9*wV,1) + 0.1*(Conf_DistOrnGainV')
% % wV = one+two
% % 
% % sum(one,1)
% % 
% % sum(two,1)
% 
% 
% sum(wV,1)
% 
% size(wV)
% 
% Conf_DistOrnGainV(40)
% 
% pause



end


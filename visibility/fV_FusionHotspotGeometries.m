function [ V,...
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
                                                    visualize )
%fV generates Visibility matrix V.


FusionNotPossibleDueToInsufficientCandConf = 0;




%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%       ANGLES
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%       VISIBILITY MATRIX
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% disp('.. Visibility Matrix')
% V = sparse(n,c); % initializing
V = sparse(n,n); % initializing
% V = zeros(n,c); % initializing

% [hotspots_r,hotspots_c] = ind2sub(size(map_env),hcs);
% hcs_fusion_ind = sub2ind(size(map_env),hcs_fusion_r,hcs_fusion_c);
[hcs_fusion_r,hcs_fusion_c] = ind2sub(size(map_env),hcs_fusion_ind);



oConfOrns = zeros(1,n);
% For each cell of the map.
for i = 1:n
    if map_candConf(i) % if its not occupied
        [rr,cc] = ind2sub(size(map_candConf),i);
        sensor = [rr,cc];
        
        %-------------------------------------------------------------
        %           Orientation of this conf
        %-------------------------------------------------------------
        % --
        %{
        orientation1 = round(rad2deg(angle2Points(sensor,[hotspots_r(1),hotspots_c(1)])));
        orientation1 = orientation1 + (orientation1==0)*360;
        
        orientation2 = round(rad2deg(angle2Points(sensor,[hotspots_r(2),hotspots_c(2)])));
        orientation2 = orientation2 + (orientation2==0)*360;
        
        %orientation =
        %round(wrapTo360(meanangle([orientation1,orientation2]))) % this has problem 
        %with [270,90] or similar cases
        orientation = round(wrapTo360(rad2deg(circ_mean([deg2rad(orientation1),...
                                                         deg2rad(orientation2)],[],2))));
        orientation = orientation + (orientation==0)*360;
        %}
        
        orns = zeros(1,numel(hcs_fusion_ind));
        for k = 1:numel(hcs_fusion_ind)
            
            orn_this = round(rad2deg(angle2Points(sensor,[hcs_fusion_r(k),hcs_fusion_c(k)])));
            orn_this = orn_this + (orn_this==0)*360;
            
            orns(k) = orn_this;
        end
        orientation = round(wrapTo360(rad2deg(circ_mean(deg2rad(orns),[],2))));
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
        %{
        % For each configuration of a cell.
        for j = orientation %1%:f
            
            % sensorDom := sensor domain translation w.r.t sensor position.
            sensorDom = [sD{j}(:,1)+sensor(1) sD{j}(:,2)+sensor(2)];
            
            % obstacDom := obstalce domain translated w.r.t sensor's
            % position. 
            obstacDom = [oD{j}(:,1)+sensor(1) oD{j}(:,2)+sensor(2)];
            
            % VOBS := Visibility vs Obstacles matrix for jth FoV.
            VOBS = VObs{j};
            
            % visibile cells 
            [ visDom ] = fVisibileDomain( map_coverage,map_env,sensorDom,obstacDom,VOBS );
            
            % Update Visibility matrix.
            %V(visDom,((i-1)*f)+j) = 1;
            V(visDom,i) = 1;
            
        end
        %}
        oConfOrns(i) = orientation;
        
    end
end


for i = 1:numel(hcs_fusion_ind)
    
    for j = 1:numel(fixedConfs_ind{i})
        
        [rr,cc] = ind2sub(size(map_candConf),fixedConfs_ind{i}(j));
        sensor  = [rr,cc];

        k = fixedConfs_orn{i}(j);

        sD = []; oD = []; VObs = [];
        %this_dir = sprintf(['sD_oD_VObs/sensing_range_%02dcell/fov_%03ddeg/',...
        %                    'orientation_%03ddeg/'],SensorRange.cell,FoV,k);
        this_dir = sprintf('visibility/templates/sensing_range_%02dcell/fov_%03ddeg/',...
                SensorRange.cell,FoV);
        this_file = sprintf('sD_oD_VObs_rng%02d_fov%03d_orn%03d.mat',...
                    SensorRange.cell,FoV,k);
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
        V(visDom,fixedConfs_ind{i}(j)) = 1;
        
        oConfOrns(fixedConfs_ind{i}(j)) = k;
    end
end


% sensing config can not visualize own cell
for i = 1:numel(map_candConf)    
    %V(i,(i-1)*f+1:(i-1)*f+f) = 0;
    V(i,i) = 0;
end

oV = V; % oV := visibility with obstacles




% VISIBILITY MATRIX CORRESPOND TO THE HOTSPOT ONLY:
% ------------------------------------------------------------------------------

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%       CELLS TO KEEP FOR COVERAGE
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
% cellKeptWeights = map_gd(cellKept);
% cellKeptWeights = cellKeptWeights/(max(cellKeptWeights));


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%       CANDIDATE CONFIGURATIONS TO KEEP
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
% candConfNum([fixedConfs1_ind,fixedConfs2_ind]) = 1;
candConfNum(cell2mat(fixedConfs_ind)) = 1;



% -- plot 
if visualize == 1
figure('name','conf corresponding to candidate conf map only');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
candConfNum_ind = find(candConfNum(:)==1);
[candConfNum_r,candConfNum_c] = ind2sub(size(map_env),candConfNum_ind);
plot(candConfNum_r,candConfNum_c,'xb','MarkerFaceColor','b');        
plot(hcs_fusion_r,hcs_fusion_c,'or','MarkerFaceColor','r');

% fixedConfs_ind = [fixedConfs1_ind,fixedConfs2_ind];
% [fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),fixedConfs_ind);
[fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),cell2mat(fixedConfs_ind));
plot(fixedConfs_r,fixedConfs_c,'or');

fusionConfs_ind = [fusionConf1,fusionConf2];
[fusionConfs_r,fusionConfs_c] = ind2sub(size(map_env),fusionConfs_ind);
plot(fusionConfs_r,fusionConfs_c,'og');
end




%----------------------------------------------------
% -- Check#2: confs. that cover all the hotspots
%----------------------------------------------------
conf_hcs = ones(1,size(oV,2));
%for i = 1:size(numel(hcs_fusion_ind,1))
for i = 1:numel(hcs_fusion_ind)
    
    % conf ind that cover this hotspot
    conf_ind = find(oV(hcs_fusion_ind(i),:));
    % binary vector
    conf_this = zeros(1,size(oV,2));
    conf_this(conf_ind) = 1;
    
    % - update the overall vector
    conf_hcs = conf_hcs.*conf_this;
end
% - make sure fixed conf are in the vector
% conf_hcs([fixedConfs1_ind,fixedConfs2_ind]) = 1;
conf_hotspot(cell2mat(fixedConfs_ind)) = 1;

% -- plot 
if visualize == 1
figure('name','conf cover all the hotspots');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
conf_hotspot_ind = find(conf_hcs(:)==1);
[conf_hotspot_r,conf_hotspot_c] = ind2sub(size(map_env),conf_hotspot_ind);
plot(conf_hotspot_r,conf_hotspot_c,'xb','MarkerFaceColor','b');        
plot(hcs_fusion_r,hcs_fusion_c,'or','MarkerFaceColor','r');

% fixedConfs_ind = [fixedConfs1_ind,fixedConfs2_ind];
% [fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),fixedConfs_ind);
[fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),cell2mat(fixedConfs_ind));
plot(fixedConfs_r,fixedConfs_c,'or');

fusionConfs_ind = [fusionConf1,fusionConf2];
[fusionConfs_r,fusionConfs_c] = ind2sub(size(map_env),fusionConfs_ind);
plot(fusionConfs_r,fusionConfs_c,'og');
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
% conf_free([fixedConfs1_ind,fixedConfs2_ind]) = 1;
conf_free(cell2mat(fixedConfs_ind)) = 1;


% -- plot 
if visualize == 1
figure('name','conf belong to unoccupied cells');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
conf_free_ind = find(conf_free(:)==1);
[conf_free_r,conf_free_c] = ind2sub(size(map_env),conf_free_ind);
plot(conf_free_r,conf_free_c,'xb','MarkerFaceColor','b');
plot(hcs_fusion_r,hcs_fusion_c,'or','MarkerFaceColor','r');

% fixedConfs_ind = [fixedConfs1_ind,fixedConfs2_ind];
% [fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),fixedConfs_ind);
[fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),cell2mat(fixedConfs_ind));
plot(fixedConfs_r,fixedConfs_c,'or');

fusionConfs_ind = [fusionConf1,fusionConf2];
[fusionConfs_r,fusionConfs_c] = ind2sub(size(map_env),fusionConfs_ind);
plot(fusionConfs_r,fusionConfs_c,'og');
end

% ----------------------------------------------------
% -- Check#4: conf not close enough to the hotspot 
% -----------------------------------------------------

% - cose enough radius
% innerRadius_cell = para_.innerRadius_m/dataYAML.resolution;
innerRadius_cell = para_.innerRadius_m/cellsize_env;

% temporary inner limit map
map_candConf_minLimit = zeros(size(map_env));
map_candConf_minLimit(hcs_fusion_ind) = 1;

se = strel('disk',innerRadius_cell,0);
map_candConf_minLimit = imdilate(map_candConf_minLimit,se);

% indices of inner limit map
innerLimit_ind = find(map_candConf_minLimit);

% -- fixedconf to be removed from inner limit map %setdiff
% fixedConfs_ind = [fixedConfs1_ind,fixedConfs2_ind]; % indices of fixed conf
% fixedConfs_cell = fixedConfs_ind; % cells of fixed conf
fixedConfs_cell = cell2mat(fixedConfs_ind); % cells of fixed conf
sacred_ind = ismember(innerLimit_ind,fixedConfs_cell); % conflecting indices 
innerLimit_ind(sacred_ind) = [];

% -- plot 
if visualize == 1
figure('name','conf inner limit');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
[innerLimit_r,innerLimit_c] = ind2sub(size(map_env),innerLimit_ind);
plot(innerLimit_r,innerLimit_c,'xb','MarkerFaceColor','b');
plot(hcs_fusion_r,hcs_fusion_c,'sr');

% fixedConfs_ind = [fixedConfs1_ind,fixedConfs2_ind];
% [fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),fixedConfs_ind);
[fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),cell2mat(fixedConfs_ind));
plot(fixedConfs_r,fixedConfs_c,'or');

fusionConfs_ind = [fusionConf1,fusionConf2];
[fusionConfs_r,fusionConfs_c] = ind2sub(size(map_env),fusionConfs_ind);
plot(fusionConfs_r,fusionConfs_c,'og');
end

confNotClose = ones(1,size(oV,2));
confNotClose(innerLimit_ind) = 0;
confNotClose(cell2mat(fixedConfs_ind)) = 1;


% -- plot 
if visualize == 1
figure('name','conf Not Close');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])

[confNotClose_r,confNotClose_c] = ind2sub(size(map_env),find(confNotClose));
plot(confNotClose_r,confNotClose_c,'xb','MarkerFaceColor','b');
plot(hcs_fusion_r,hcs_fusion_c,'sr');

% fixedConfs_ind = [fixedConfs1_ind,fixedConfs2_ind];
% [fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),fixedConfs_ind);
[fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),cell2mat(fixedConfs_ind));
plot(fixedConfs_r,fixedConfs_c,'or');

fusionConfs_ind = [fusionConf1,fusionConf2];
[fusionConfs_r,fusionConfs_c] = ind2sub(size(map_env),fusionConfs_ind);
plot(fusionConfs_r,fusionConfs_c,'og');
end

%-----------------------------------------------------------------------------
% -- Check#5: confs (position) is not already being used for the other hcs
%-----------------------------------------------------------------------------
confNotEngaged = ones(1,size(oV,2));
confNotEngaged(engagedConfs_ind) = 0;
confNotEngaged(cell2mat(fixedConfs_ind)) = 1;

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
conf_vector  = candConfNum.*conf_hcs.*conf_free.*confNotClose.*confNotEngaged;

% candConfNumXX = numel(find(candConfNum))
% conf_hcsXX = numel(find(conf_hcs))
% conf_freeXX = numel(find(conf_free))
% confNotCloseXX = numel(find(confNotClose))
% confNotEngagedXX = numel(find(confNotEngaged))

confKept     = find(conf_vector);
confRemoved  = find(conf_vector==0);
% V            = V(:,confKept);
% V            = full(V);




%------------------------------------------
%   IF NO CANDIDATE CONF LEFT
%------------------------------------------
if numel(confKept) < 1
    FusionNotPossibleDueToInsufficientCandConf = 1;
    V                   = [];
    oV                  = [];
    wV                  = [];
    Conf_DistOrnGainV   = [];
    confKept            = [];
    fixedConfs_num      = [];
    fixedConfXG_num     = [];
    confRemoved         = [];
    cellKept            = [];
    cellRemoved         = [];
    oConfOrns           = [];
    FoVfaces            = [];
    warning('Breaking due to no candidate conf available.');
    return
end

% -- plot 
if visualize == 1
figure('name','conf meeting all four conditions');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
[confKept_r,confKept_c] = ind2sub(size(map_env),confKept);
h_confKept = plot(confKept_r,confKept_c,'xb','MarkerFaceColor','b');
h_hcs4Fusion = plot(hcs_fusion_r,hcs_fusion_c,'or','MarkerFaceColor','r');

% fixedConfs_ind = [fixedConfs1_ind,fixedConfs2_ind];
% [fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),fixedConfs_ind);
[fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),cell2mat(fixedConfs_ind));
h_fixedConf = plot(fixedConfs_r,fixedConfs_c,'or');

fusionConfs_ind = [fusionConf1,fusionConf2];
[fusionConfs_r,fusionConfs_c] = ind2sub(size(map_env),fusionConfs_ind);
h_fusionConf = plot(fusionConfs_r,fusionConfs_c,'og');


hLegend = legend([h_confKept,h_hcs4Fusion,h_fixedConf,h_fusionConf],...
                'conf kept',...
                'hcs for fusion',...
                'fixed conf',...        
                'fusion conf',...
                'location','best','Orientation','vertical');
set(hLegend,'Interpreter','latex');
%legend('boxoff');
set([hLegend,gca],'FontSize',7);
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%
%       LIMIT THE NUMBER OF CONF BASED ON PROPORTIONAL COVERAGE
%
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%/////////////////////////////////////////////////////////////
%
%       PROPORTIONAL COVERAGE WEIGHTS
%
%/////////////////////////////////////////////////////////////

cellKeptWeights = map_recon(cellKept);
% cellKeptWeights = cellKeptWeights/(max(cellKeptWeights));


% confKept
V_this = V(:,confKept);
V_this = V_this .* repmat(cellKeptWeights',1,numel(confKept));
wV_concentration_this = sum(V_this,1);



 wV_concentration_this
 size(wV_concentration_this)
 max(wV_concentration_this(:))
wV_concentration_this = wV_concentration_this/max(wV_concentration_this(:));


% -- make up for fixedconf
% fixedConfs_ind = [fixedConfs1_ind,fixedConfs2_ind];
fixed_ind = ismember(confKept,cell2mat(fixedConfs_ind));
wV_concentration_this(fixed_ind) = max(wV_concentration_this(:)+0.5);


%/////////////////////////////////////////////////////////////
%
%       PERSISTENT COVERAGE WEIGHTS
%
%/////////////////////////////////////////////////////////////

%---------------------------------
%
%       ORIENTATION
%
%---------------------------------
% -- conf orientations
ConfOrns = oConfOrns(confKept);

% -- desired orinetations w.r.t to hotspot
% [hotspot_r,hotspot_c] = ind2sub(size(map_env),hotspot);
[conf_r,conf_c] = ind2sub(size(map_candConf),confKept);

% desiredOrns1 = round(rad2deg(angle2Points([conf_r',conf_c'],[hcs_fusion_r(1),hcs_fusion_c(1)])));
% desiredOrns1 = desiredOrns1 + (desiredOrns1==0)*360;
% 
% desiredOrns2 = round(rad2deg(angle2Points([conf_r',conf_c'],[hcs_fusion_r(2),hcs_fusion_c(2)])));
% desiredOrns2 = desiredOrns2 + (desiredOrns2==0)*360;

desiredOrns = cell(1,numel(hcs_fusion_ind));
for i = 1:numel(hcs_fusion_ind)
    desiredOrns{i} = round(rad2deg(angle2Points([conf_r',conf_c'],[hcs_fusion_r(i),hcs_fusion_c(i)])));
    desiredOrns{i} = desiredOrns{i} + (desiredOrns{i}==0)*360;    
end

% -- delta orientations
% ConfDeltaOrns1 = rad2deg(angleAbsDiff(deg2rad(desiredOrns1(:)),deg2rad(ConfOrns(:))));
% ConfDeltaOrns2 = rad2deg(angleAbsDiff(deg2rad(desiredOrns2(:)),deg2rad(ConfOrns(:))));

ConfDeltaOrns = cell(1,numel(hcs_fusion_ind));
for i = 1:numel(hcs_fusion_ind)
    ConfDeltaOrns{i} = rad2deg(angleAbsDiff(deg2rad(desiredOrns{i}(:)),deg2rad(ConfOrns(:))));
end


% -- angular gain
gain_sigma = FoV/6;%15; % standard deviation
gain_meu   = 0; % expected value or center
% V_angular_gain1_this = gaussmf(ConfDeltaOrns1,[gain_sigma,gain_meu]);
% V_angular_gain2_this = gaussmf(ConfDeltaOrns2,[gain_sigma,gain_meu]);

V_angular_gain_this = cell(1,numel(hcs_fusion_ind));
for i = 1:numel(hcs_fusion_ind)
    V_angular_gain_this{i} = gaussmf(ConfDeltaOrns{i},[gain_sigma,gain_meu]);
end


%---------------------------------
%
%       DISPLACEMENT
%
%---------------------------------

% -- dist from conf to hotspot
% conf2Hot_dist1 = pdist2([conf_r',conf_c'],[hcs_fusion_r(1),hcs_fusion_c(1)],'euclidean');
% conf2Hot_dist2 = pdist2([conf_r',conf_c'],[hcs_fusion_r(2),hcs_fusion_c(2)],'euclidean');

conf2Hot_dist = cell(1,numel(hcs_fusion_ind));
for i = 1:numel(hcs_fusion_ind)
    conf2Hot_dist{i} = pdist2([conf_r',conf_c'],[hcs_fusion_r(i),hcs_fusion_c(i)],'euclidean');
end

% -- distance gain
gain_sigma       = SensorRange.cell/5; % standard deviation
gain_meu         = SensorRange.cell/2; % expected value or center
% V_dist_gain1_this = gaussmf(conf2Hot_dist1,[gain_sigma,gain_meu]);
% V_dist_gain2_this = gaussmf(conf2Hot_dist2,[gain_sigma,gain_meu]);

V_dist_gain_this = cell(1,numel(hcs_fusion_ind));
for i = 1:numel(hcs_fusion_ind)
    V_dist_gain_this{i} = gaussmf(conf2Hot_dist{i},[gain_sigma,gain_meu]);
end




%---------------------------------
%
%   ORIENTATION + DISPLACEMENT
%
%---------------------------------
% wV_diplacement1_this = (V_angular_gain1_this+V_dist_gain1_this)/2;
% wV_diplacement2_this = (V_angular_gain2_this+V_dist_gain2_this)/2;

wV_diplacement_this = cell(1,numel(hcs_fusion_ind));
for i = 1:numel(hcs_fusion_ind)
    wV_diplacement_this{i} = (V_angular_gain_this{i}+V_dist_gain_this{i})/2;
end

% wV1_this = (gamma*(wV_concentration_this))' + ((1-gamma)* wV_diplacement1_this);
% wV2_this = (gamma*(wV_concentration_this))' + ((1-gamma)* wV_diplacement2_this);
% 
% wV_this = (wV1_this + wV2_this)/2;

wV_this = cell(1,numel(hcs_fusion_ind));
for i = 1:numel(hcs_fusion_ind)
    wV_this{i} = (para_.gamma*(wV_concentration_this))' + ...
                 ((1-para_.gamma)* wV_diplacement_this{i});
end

wV_all = zeros(size(wV_this{1}));
for i = 1:numel(hcs_fusion_ind)
    wV_all = wV_all+wV_this{i};
end
wV_all = wV_all/numel(hcs_fusion_ind);

% -- sorted indices
% [~,sorted_ind] = sort(wV_this,'descend');
[~,sorted_ind] = sort(wV_all,'descend');



switch para_.MissionStrategy
    case '1t-armex'
        para_.maxCandidateConf_Fusion = para_.maxCandidateConf_Fusion_1t;
    case '2t-armex'
        para_.maxCandidateConf_Fusion = para_.maxCandidateConf_Fusion_2t;
end
sorted_ind_max = sorted_ind(1:min([para_.maxCandidateConf_Fusion,...
    numel(sorted_ind)]));
% sorted_ind_max = sorted_ind(1:min([10,numel(sorted_ind)]));
confToRemove_ind = find(ismember(sorted_ind,sorted_ind_max)==0);
confRemoved = union(confRemoved,sorted_ind(confToRemove_ind));
confKept = confKept(sorted_ind_max);


%------------------------------------------
%   IF NO CANDIDATE CONF LEFT
%------------------------------------------
if numel(confKept) < 1
    FusionNotPossibleDueToInsufficientCandConf = 1;
    V                   = [];
    oV                  = [];
    wV                  = [];
    Conf_DistOrnGainV   = [];
    confKept            = [];
    fixedConfs_num      = [];
    fixedConfXG_num     = [];
    confRemoved         = [];
    cellKept            = [];
    cellRemoved         = [];
    oConfOrns           = [];
    FoVfaces            = [];
    warning('Breaking due to no candidate conf available.');
    return
end


%--- hardcore fixed confs
confKept = union(confKept,cell2mat(fixedConfs_ind));

% -- sorted list
confKept = sort(confKept);


% -- plot 
if visualize == 1
figure('name','max num of conf');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
[confKept_r,confKept_c] = ind2sub(size(map_env),confKept);
plot(confKept_r,confKept_c,'xb','MarkerFaceColor','b');
plot(hcs_fusion_r,hcs_fusion_c,'sr');

% fixedConfs_ind = [fixedConfs1_ind,fixedConfs2_ind];
% [fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),fixedConfs_ind);
[fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),cell2mat(fixedConfs_ind));
plot(fixedConfs_r,fixedConfs_c,'or');

fusionConfs_ind = [fusionConf1,fusionConf2];
[fusionConfs_r,fusionConfs_c] = ind2sub(size(map_env),fusionConfs_ind);
plot(fusionConfs_r,fusionConfs_c,'og');
end


%****************************************************************
%
%     conf kept visibility matrix
%
%****************************************************************

V = V(:,confKept);
V = full(V);

% size(V)

%****************************************************************
%
%     weighted visibility matrix
%
%****************************************************************


% wV = V.*repmat(cellKeptWeights',1,numel(confKept));
% wV = V.*repmat(map_recon(cellKept)',1,numel(confKept));
% max_conc_conf = max(sum(wV,1));
% wV = wV./max_conc_conf;
WV = V.*repmat(map_recon(cellKept)',1,numel(confKept));
max_conc_conf = max(sum(WV,1));
WV = WV./max_conc_conf;

% -- weighted visibility for fixedCo  nf1 and fixedConf2
% fixedConf1_num = find(ismember(confKept,fixedConfs1_ind));
% fixedConf2_num = find(ismember(confKept,fixedConfs2_ind));

fixedConfs_num = cell(1,numel(hcs_fusion_ind));
for i = 1:numel(hcs_fusion_ind)
    fixedConfs_num{i} = find(ismember(confKept,fixedConfs_ind{i}));
end


% fixed confs of other local solutions (hcs_not_this) that can not cover
% this hc, can not be used for the visibility, cross angles, and distance
% gains.
% 
% fixedConfXG_num:= fixed conf num that must not be used for the gain
% calculations of this hc.
fixedConfXG_num = cell(1,numel(hcs_fusion_ind));
wV = cell(1,numel(hcs_fusion_ind));
num_all = 1:numel(hcs_fusion_ind);

for i = 1:numel(hcs_fusion_ind)
    wV{i} = WV;
    fixedConfOtherThanThisHc_num = cell2mat(fixedConfs_num(num_all~=i));
    CellNumThisHc = find(cellKept==hcs_fusion_ind(i));
    %this_ind = wV{i}(cellKept==hc_fusion_ind(i),cell2mat(fixedConfs_num(num_all~=i)))==0
    this_ind = wV{i}(CellNumThisHc,fixedConfOtherThanThisHc_num)==0;
    fixedConfXG_num{i} = fixedConfOtherThanThisHc_num(this_ind);
    wV{i}(:,fixedConfXG_num{i}) = 0;
end
    

% % wV1 = wV;
% % wV1(:,fixedConf2_num) = 0;
% % hotspotCellKept = find(cellKept==hotspots(1))
% % fixedConf2XG_num = fixedConf2_num(find(wV1(hotspotCellKept,fixedConf2_num)==0))
% fixedConf2XG_num = fixedConf2_num(wV1(cellKept==hcs(1),fixedConf2_num)==0);
% wV1(:,fixedConf2XG_num) = 0;
% 
% 
% wV2 = wV;
% % wV2(:,fixedConf1_num) = 0;
% fixedConf1XG_num = fixedConf1_num(wV2(cellKept==hcs(2),fixedConf1_num)==0);
% wV2(:,fixedConf1XG_num) = 0;

% -- fixed confs that do not cover the 2nd hotspot in fusion
%fixedConf1XG_num 
%fixedConf2XG_num

%fixedConf1_num
%fixedConf2_num



%//////////////////////////////////////////
% HOTSPOT PLACEMENT GAIN - ANGULAR
%//////////////////////////////////////////

% -- conf orientations
ConfOrns = oConfOrns(confKept);


% -- desired orinetations w.r.t to hotspot
% [hotspots_r,hotspots_c] = ind2sub(size(map_env),hotspots);
[conf_r,conf_c] = ind2sub(size(map_candConf),confKept);

% desiredOrns1 = round(rad2deg(angle2Points([conf_r',conf_c'],[hcs_fusion_r(1),hcs_fusion_c(1)])));
% desiredOrns1 = desiredOrns1 + (desiredOrns1==0)*360;
% 
% desiredOrns2 = round(rad2deg(angle2Points([conf_r',conf_c'],[hcs_fusion_r(2),hcs_fusion_c(2)])));
% desiredOrns2 = desiredOrns2 + (desiredOrns2==0)*360;

desiredOrns = cell(1,numel(hcs_fusion_ind));
for i = 1:numel(hcs_fusion_ind)
    desiredOrns{i} = round(rad2deg(angle2Points([conf_r',conf_c'],[hcs_fusion_r(i),hcs_fusion_c(i)])));
    desiredOrns{i} = desiredOrns{i} + (desiredOrns{i}==0)*360;
end

% -- delta orientations
% ConfDeltaOrns1 = rad2deg(angleAbsDiff(deg2rad(desiredOrns1(:)),deg2rad(ConfOrns(:))));
% ConfDeltaOrns2 = rad2deg(angleAbsDiff(deg2rad(desiredOrns2(:)),deg2rad(ConfOrns(:))));

ConfDeltaOrns = cell(1,numel(hcs_fusion_ind));
for i = 1:numel(hcs_fusion_ind)
    ConfDeltaOrns{i} = rad2deg(angleAbsDiff(deg2rad(desiredOrns{i}(:)),deg2rad(ConfOrns(:))));
end

% -- angular gain
gain_sigma = FoV/6;%15; % standard deviation
gain_meu   = 0; % expected value or center
% V_angular_gain1 = gaussmf(ConfDeltaOrns1,[gain_sigma,gain_meu]);
% V_angular_gain2 = gaussmf(ConfDeltaOrns2,[gain_sigma,gain_meu]);

V_angular_gain = cell(1,numel(hcs_fusion_ind));
for i = 1:numel(hcs_fusion_ind)
    V_angular_gain{i} = gaussmf(ConfDeltaOrns{i},[gain_sigma,gain_meu]);
end


%//////////////////////////////////////////
% HOTSPOT PLACEMENT GAIN - DISTANCE
%//////////////////////////////////////////

% -- dist from conf to hotspot
% conf2Hot_dist1 = pdist2([conf_r',conf_c'],[hcs_fusion_r(1),hcs_fusion_c(1)],'euclidean');
% conf2Hot_dist2 = pdist2([conf_r',conf_c'],[hcs_fusion_r(2),hcs_fusion_c(2)],'euclidean');

conf2Hot_dist = cell(1,numel(hcs_fusion_ind));
for i = 1:numel(hcs_fusion_ind)
    conf2Hot_dist{i} = ...
        pdist2([conf_r',conf_c'],[hcs_fusion_r(i),hcs_fusion_c(i)],'euclidean');
end

% -- distance gain
gain_sigma = SensorRange.cell/5; % standard deviation
gain_meu   = SensorRange.cell/2; % expected value or center
% V_dist_gain1 = gaussmf(conf2Hot_dist1,[gain_sigma,gain_meu]);
% V_dist_gain2 = gaussmf(conf2Hot_dist2,[gain_sigma,gain_meu]);

V_dist_gain = cell(1,numel(hcs_fusion_ind));
for i = 1:numel(hcs_fusion_ind)
    V_dist_gain{i} = gaussmf(conf2Hot_dist{i},[gain_sigma,gain_meu]);
end


%//////////////////////////////////////////////
% HOTSPOT PLACEMENT GAIN - ANGULAR + DISTANCE
%//////////////////////////////////////////////
% Conf_DistOrnGainV1 = (V_angular_gain1+V_dist_gain1)/2;
% Conf_DistOrnGainV2 = (V_angular_gain2+V_dist_gain2)/2;

Conf_DistOrnGainV = cell(1,numel(hcs_fusion_ind));
for i = 1:numel(hcs_fusion_ind)
    Conf_DistOrnGainV{i} = (V_angular_gain{i}+V_dist_gain{i})/2;
end





end


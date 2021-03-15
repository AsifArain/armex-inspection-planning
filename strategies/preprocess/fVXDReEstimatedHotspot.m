function [ wV,...
           oV,...
           D,...           
           confKept,...
           map_candConf,...
           Conf_DistOrnGainV,...
           conf_crossAngles,...
           conf_crossAngles_G,...
           oConfOrientations,...
           executedConfsHot_num,...
           tPreprocess ] = fVXDReEstimatedHotspot( map_env,...
                                                   FoV,...
                                                   hotspot_current,...
                                                   executedConfsHot_ind,...
                                                   executedConfsHot_orn,...
                                                   OBSenv,...
                                                   OBSconf,...
                                                   gamma,...
                                                   SensorRange,...
                                                   map_coverage,...
                                                   map_gd_current,...
                                                   refPointDistance,...
                                                   dataYAML,...
                                                   para_ )
% fPreProcess_R4 genereates required parameters for solving a global coverage problem.
% Date: 2014-01-19, Rev 4
% 
% n  := total number of cells in the map.
% nf := number of free cells in the map.
% f  := total number of sensing configurations per cell.
% 
% OUTPUTS:
% ----------------------------------------------------------------------------------------
% cnt: [scalar][ ]   4 or 8 connectivity graph.
% o  : [scalar][ ]   Number of outgoing edges for each conf.
% E  : [nf*f*o][ ]   Edge connectivity. Index number is ID of outgoing edge and array value is ID of connected incoming edge.%                  
% vE : [------][ ]   E (as above) with reduced variables using Corner Reduction Method.
% eE : [nf*f*o][ ] 
% evE: [------][ ]   eE (as above) with reduced variables using Corner Reduction Method.
% T  : [nf*f,o][sec] Travelling cost. Rows are conf number and columns are outgoing edges.
% vT : [------][sec] T (as above) with reduced variables using Corner Reduction Method.
% V  : [nf,nf*f][binary] Visibility matrix, indicates visiblity status (1 for visible, 0 otherwise) of nf cells (rows) from
%                        each conf (columns).
% vV : [-------][binary] V (as above) with reduced variables using Corner Reduction Method.
% vkS: [][]          .cell := selected special vertex cell number.
%                    .conf := configurations that cover selected special vertex cell.
% CRL: [][]          Corner Reduction List, list of configuration selected to be removed using Corner Reduction technique.
% SensorRange: [scalar][integer] .cell sensor range in number of cells.
% map: [l,m][binary] Grid map matrix, 0 for occupied cells and 1 for unoccupied cells. l and m are any positive integer numbers.
% FoVfaces [][]: .alf first/start angle for each FoV orientation.
%                .bay second/final angle for each FoV orientation.
%                .ang set of first/start and second/final angles of each FoV orientation.
%                .num total number of FoV faces/orientations for a sensing position.
% OBS: [][]      .ind Indices of occupied cells.
%                .sub Subscripts (row,col) of occupied cells.
% tPreprocess: [scalar][sec] Total computation time.
% 
% INPUTS:
% --------------------------------------------------------------------------------------------------------------------------------
% FoV          : [][integer] Field of view.
% SensorRange  : [scalar][integer] .m sensor range in meters.
% paramPreprocess     : [][] Solution parameters.
%                     .CRL      := compute list of corner reduction variables (yes/no).
%                     .vkS      := selection of special vertex (yes/no).
%                     .Ct       := compute travelling cost matrix (yes/no).
%                     .GraphCnt := graph connectivity (4/8).
%                     .Res      := required cell size in meter.
%                     .TravCost := travelling cost (sec) for one meter travel.
%                     .PixSize  := original (map) pixel size (meters).
%                     .numFoVs  := number of elementray sensing action per cell.
%                     .A        := generate connectivity matrix A (yes/no).
%                     .AWM      := artificial world map (yes/no).
%                     .RWM      := real world map (yes/no).
%                     .MapFig   := map file e.g. 'frieburg.png'.
%                     .MapTxt   := map file e.g. 'kjula.txt'.
%                     .SenDom   := computing sensor domain. 'ComputeNew' to compute from scratch, and 'UseEarlier' to use
%                                  earlier computed sensor domain parameters.
% map          : [l,m][binary] Grid map matrix, 0 for occupied cells and 1 for unoccupied cells. l and m are any positive
%                              integer numbers. 
%
% NOTES/UPDATES:
% --------------------------------------------------------------------------------------------------------------------------------
% R9 : RAGT::CrossAngles.
% R8 : ...
% R7 : ...
% R6 : ...
% R5 : ...
% R4 : Visibility Method VIS4 (2013-12-13) is used.
% R3 : Visibility Method VIS3 (2013-12-04) is used.
% R2 : ... 
% R1 : ...
% R0 : ...
% 
% oV = Visibility with obstacles
% xV or V = Visibility without obstacles.
% 
%  Map (plot)
%       -------------
%       | 7 | 8 | 9 |
%       -------------
%       | 4 | 5 | 6 |
%       -------------
%       | 1 | 2 | 3 |
%       -------------
%  Map (matrix)
%      =  1 4 7
%         2 5 8
%         3 6 9


tPreprocessi = tic; % initialize pre computational time.
visualize = 0;

%% SENSING PARAMETERS AND CONSTANTS:
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
FoVfaces.num        = para_.numFoVs;   % number of conf per cell.
n                   = numel(map_env);         % number of cells in the map.
f                   = FoVfaces.num;       % for simplification.
c                   = n*f;                % total conf.
% cnt                 = para_.GraphCnt;  % 4 or 8 connected neighbours.

% o := number of outgoing edges for each configuration.



%% CANDIDATE CONF MAP (FINAL) FROM 1ST2ND MAP:
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

hotspot_current = round(hotspot_current);

%////////////////////////////////////////////////////////
% if hotspot is on the obstacle
%////////////////////////////////////////////////////////

if map_env(hotspot_current(1),hotspot_current(2)) ==0
    
    coverage_ind = find(map_env(:)==1);
    [cov_r,cov_c] = ind2sub(size(map_env),coverage_ind);
    
    dist = pdist2(hotspot_current,[cov_r,cov_c],'euclidean');
    [~,ind] = min(dist);
    
    shifted_ind = coverage_ind(ind);
    [shifted_r,shifted_c] = ind2sub(size(map_env),shifted_ind);
    hotspot_current = [shifted_r,shifted_c];
    
end
    


hotspots = sub2ind(size(map_env),hotspot_current(1),hotspot_current(2));

% [hotspot_r,hotspot_c] = ind2sub(size(map_env),hotspot);
hotspots_r = hotspot_current(1);
hotspots_c = hotspot_current(2);



%=======================================================================
% Candidate Conf map based on sensing range
%=======================================================================
[ map_candConf ] = fConfMapSensingRange( hotspots,...
                                         map_env,...
                                         SensorRange,...
                                         OBSconf,...
                                         visualize);


%% VISIBILITY MATRIX:
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

fprintf('Visibility...\n')

[  V,...
   oV,...
   wV,...
   Conf_DistOrnGainV,...
   confKept,...
   confRemoved,...
   cellKept,...           
   cellRemoved,...
   oConfOrientations,...
   executedConfsHot_num,...
   FoVfaces ] = fV( map_candConf,...
                    map_coverage,...
                    map_env,...
                    map_gd_current,...
                    n,...
                    f,...
                    c,...
                    FoV,...
                    FoVfaces,...
                    SensorRange,...
                    executedConfsHot_ind,...
                    executedConfsHot_orn,...
                    OBSenv,...
                    OBSconf,...
                    hotspots,...
                    gamma,...
                    dataYAML,...
                    visualize,...
                    para_ );
               
%=======================================================================
% -- CANDIDATE CONF MAP (FINAL)
%=======================================================================

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
plot(hotspots_r,hotspots_c,'or','MarkerFaceColor','r');

fixedConfs_ind = executedConfsHot_ind;
[fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),fixedConfs_ind);
plot(fixedConfs_r,fixedConfs_c,'or');

end


%=======================================================================
% CROSS ANGLES
%=======================================================================


fprintf('Cross angles...\n')

[ Hot2Conf_Angles,...
  conf_crossAngles,...
  conf_crossAngles_G] = fX( confKept,...
                             hotspots,...
                             map_env );

%=======================================================================
% DISTANCE 
%=======================================================================
fprintf('Distance...\n')

[ D ] = fD( refPointDistance,...
            confKept,...
            map_env );


%% Total computation time.
% ----------------------------------------------------------------------------------------
tPreprocess =  1e-4*(round(toc(tPreprocessi)*1e4));

disp(['tPreprocess: ',num2str(tPreprocess),' sec']);


end



function [ D ] = fD( refPointDistance,...
                     confKept,...
                     map_env )

refDistPoint_r = refPointDistance(1); refDistPoint_c = refPointDistance(2);


% ---- D (distance) matrix is manhattan distance between cand conf position
% and reference position
confKeptCells_ind = confKept; %fix(confKept/f)+(~~mod(confKept,f));
[confKeptCells_r,confKeptCells_c] = ind2sub(size(map_env),confKeptCells_ind);

D = zeros(numel(confKeptCells_ind),1);

for i = 1:numel(confKeptCells_ind)
    manhattan_distance = abs(confKeptCells_r(i)-refDistPoint_r)+...
        abs(confKeptCells_c(i)-refDistPoint_c);
    
    D(i) = manhattan_distance;

end


end





function [ map_candConf ] = fConfMapSensingRange( hotspots,...
                                                  map_env,...
                                                  SensorRange,...
                                                  OBSconf,...
                                                  visualize)

[hotspots_r,hotspots_c] = ind2sub(size(map_env),hotspots);

% initialize map for candidate conf. -- ring only
map_candConf = zeros(size(map_env));    
map_candConf(hotspots) = 1;

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
plot(hotspots_r,hotspots_c,'xr','MarkerFaceColor','r');
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
plot(hotspots_r,hotspots_c,'xr','MarkerFaceColor','r');
end

% extract edges
map_candConf_edges = bwmorph(map_candConf,'remove');
% -- remove edges
map_candConf(map_candConf_edges(:)==1) = 0;


% plot
if visualize == 1
pause(0.1)
figure('name','cand conf - extracted edges'); 
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
map_candConf_ind = find(map_candConf(:)==1);
[map_candConf_r,map_candConf_c] = ind2sub(size(map_env),map_candConf_ind);
plot(map_candConf_r,map_candConf_c,'xb','MarkerFaceColor','b');        
plot(hotspots_r,hotspots_c,'xr','MarkerFaceColor','r');
end

%{
% grow edges with prdefiend thickness 
% circleThickness_m = 2.5
% circleThickness_cell = circleThickness_m/paramPreprocess.MapRes
circleThickness_cell = ceil(SensorRange.cell/4);
se = strel('ball',circleThickness_cell,0); % outer limit

map_candConf = imdilate(map_candConf,se);


% plot
pause(0.1)
figure('name','cand conf - edges are increased to the limit');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
map_candConf_ind = find(map_candConf(:)==1);
[map_candConf_r,map_candConf_c] = ind2sub(size(map_env),map_candConf_ind);
plot(map_candConf_r,map_candConf_c,'xb','MarkerFaceColor','b');        
plot(hotspots_r,hotspots_c,'xr','MarkerFaceColor','r');


% remove obstacle 
%map_candConf_1st2nd(OBS.ind(find(ismember(OBS.ind,find(map_candConf_1st2nd))))) = 0;

% ---- lets have new obs list by growing the obstacles one meter
thikness_cells = ceil(1.5/dataYAML.resolution)
se = strel('square',thikness_cells); % outer limit
map_grown_obs_list = ~map_env;
map_grown_obs_list = imdilate(map_grown_obs_list,se);
map_grown_obs_list = ~map_grown_obs_list;


pause(0.1)
figure('name','map_grown_obs_list'); 
imshow((map_env-map_grown_obs_list)','InitialMagnification',800);
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])


obs_list_grown = find(~map_grown_obs_list);
map_candConf(obs_list_grown(find(ismember(obs_list_grown,find(map_candConf))))) = 0;

% plot
pause(0.1)
figure('name','cand conf - occupied cells are removed');
imshow(map_env','InitialMagnification',800); hold on;
map_candConf_ind = find(map_candConf(:)==1);
[map_candConf_r,map_candConf_c] = ind2sub(size(map_env),map_candConf_ind);
plot(map_candConf_r,map_candConf_c,'xb','MarkerFaceColor','b');        
plot(hotspots_r,hotspots_c,'xr','MarkerFaceColor','r');
%}


end






function [ Hot2Conf_Angles,...
           conf_crossAngles,...
           conf_crossAngles_G] = fX( confKept,...
                                      hotspots,...
                                      map_env )

% conf_theta            := conf orientations [1,c]
% conf_crossAngles      := cross angles between conf [c,c]
% conf_crossAngles_G    := cross angles gain [c,c]



Hot2Conf_Angles = zeros(numel(hotspots),numel(confKept));

confKeptCell_ind = confKept;
[confKeptCell_r,confKeptCell_c] = ind2sub(size(map_env),confKeptCell_ind);
sensors = [confKeptCell_r',confKeptCell_c'];

[hotspots_r,hotspots_c] = ind2sub(size(map_env),hotspots);


for i = 1:numel(hotspots)
    
    %Hot2Conf_Angles(i,:) = rad2deg(angle2Points([hotspots_r(i),hotspots_c(i)],sensors));
    Hot2Conf_Angles(i,:) = angle2Points([hotspots_r(i),hotspots_c(i)],sensors);
end


% ----- cross angles 1
disp('----- cross angles')
conf_crossAngles = zeros(size(Hot2Conf_Angles,2));

for i = 1:numel(confKept)
    for j = 1:numel(confKept)
        
        conf_crossAngles(i,j) = rad2deg(angleAbsDiff(Hot2Conf_Angles(1,i),Hot2Conf_Angles(1,j)));

    end
end

% ------ cross angles gain 1
disp('----- gain')
% gain_sigma = deg2rad(10); % standard deviation
% gain_meu   = deg2rad(60); % expected value or center
gain_sigma = 10; % standard deviation
gain_meu   = 60; % expected value or center
conf_crossAngles_G1 = gaussmf(conf_crossAngles,[gain_sigma,gain_meu]);


% gain_sigma = deg2rad(10); % standard deviation
% gain_meu   = deg2rad(120); % expected value or center
gain_sigma = 10; % standard deviation
gain_meu   = 120; % expected value or center
conf_crossAngles_G2 = gaussmf(conf_crossAngles,[gain_sigma,gain_meu]);

conf_crossAngles_G = conf_crossAngles_G1+conf_crossAngles_G2;




end







function [ V,...
           oV,...
           wV,...
           Conf_DistOrnGainV,...
           confKept,...
           confRemoved,...
           cellKept,...           
           cellRemoved,...
           oConfOrns,...
           executedConfsHot_num,...
           FoVfaces ] = fV( map_candConf,...
                            map_coverage,...
                            map_env,...
                            map_recon,...
                            n,...
                            f,...
                            c,...
                            FoV,...
                            FoVfaces,...
                            SensorRange,...
                            executedConfsHot_ind,...
                            executedConfsHot_orn,...
                            OBSenv,...
                            OBSconf,...
                            hotspots,...
                            gamma,...
                            dataYAML,...
                            visualize,...
                            para_ )
%fV generates Visibility matrix V.

% --- ANGLES:
FoVfaces.lgt = [1/FoVfaces.num:1/FoVfaces.num:1-(1/FoVfaces.num) 0]'*360;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    SensorDomain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('.. Sensor Domain')

switch para_.SensorDomainComputation
    
    case 'ComputeNew'
        
        % For each cell of the map.
        sD = cell(1,f); % sensor domain (visible).
        oD = cell(1,f); % sensor's osbtacle domain (used to check the 
                        % visibility of "sD" for each cell in "oD"]

        % a single reference cell in the map.
        for i = 1
            
            % sensor position.
            sensor = [0,0];
            
            % For each sensing configuration.
            for j = 1:f

                % zawia := initial, middle and final (opening) angle.
                zawia = [FoVfaces.ang(j,1) FoVfaces.lgt(j) FoVfaces.ang(j,2)];

                % Finding sensor domain, i.e. visible cells for the
                % configuration. 
                [ visibleRange ] = fSenDomain( sensor,SensorRange,zawia );
                % visibleRange := visible range (row,col).

                % post-processing for sorting and removing duplicates
                visibleRange = unique(visibleRange,'rows'); % remove repeated numbers
                visibleRange = sortrows(visibleRange);      % sort rows
                
                % here we have jth sensor's domain.
                sD{j} = visibleRange;
                
                
                % since, to varify if a cell in sensor's domain is visible,
                % we need to check all the cells that effect the
                % visibility. Some of them are even not in the sensing
                % range, therefore, we need to compute obstacle domain that
                % contains the cells in the sensor domain and connecting.
                % "obsDomain" will be used to check visibility of the
                % sensor's domain.

                % obstacle domain
                [ obsRange ] = fSenDomainOBS( sensor,SensorRange,zawia );
                % obsRange := obstacle range (row,col).
                
                % post-processing for sorting and removing duplicates
                obsRange = unique(obsRange,'rows'); % remove repeated numbers
                obsRange = sortrows(obsRange);      % sort rows
                
                % here we have jth sensor's domain.
                oD{j} = obsRange;
                
            end
        end

        %--------------------------------------
        % EFFECTED VISIBILITY BY EACH OBSTACLE
        %--------------------------------------
        
        disp('.. Visibility vs Obstacles')
        VObs = cell(1,f); % Visibility vs Obstacle
        
        for i = 1
            sensor = [0,0];
            % For each configuration of a cell.
            for j = 1:f
                
                % let take sensor's mask
                sensorDom  = sD{j};
                % and obstacle domain
                obstacDom  = oD{j};
                
                VObs{j} = zeros(size(obstacDom,1),size(sensorDom,1)); % no more a square matrix
                % rows    := obstacle cell,
                % columns := visibility effect.
                
                %for k = 1:size(sensorDom,1)
                for k = 1:size(obstacDom,1)
                    
                    % kth cell is an obstacle
                    %obstacDomOBS(k,3)  = 0; 
                    obs = obstacDom(k,:);
                    
                    % visibility effect on all the cells
                    [ vDom ] = fVisSDvObs( sensor,sensorDom,obs );
                    
                    % update visibility against kth obstacle
                    VObs{j}(k,:) = vDom;
                   
                end
               
            end
        end
        
        % save the mask and visibility info to a file.
        save(['preprocess/sD_oD_VObs_Tomography_Rng',...
            num2str(SensorRange.cell),'FoV',num2str(FoV)],...
            'sD','oD','VObs');
        
        pause
        
    case 'UseEarlier'
        
        load(['preprocess/sD_oD_VObs_Tomography_Rng',...
            num2str(SensorRange.cell),'FoV',num2str(FoV)]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       V
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('.. Visibility Matrix')
% V = sparse(n,c); % initializing
V = sparse(n,n); % initializing
% V = zeros(n,c); % initializing

[hotspots_r,hotspots_c] = ind2sub(size(map_env),hotspots);

oConfOrns = zeros(1,n);
% For each cell of the map.
for i = 1:n
    if map_candConf(i) % if its not occupied
        [rr,cc] = ind2sub(size(map_candConf),i);
        sensor = [rr,cc];
        
        % --
        orientation1 = round(rad2deg(angle2Points(sensor,[hotspots_r(1),hotspots_c(1)])));
        orientation1 = orientation1 + (orientation1==0)*360;
        
        %orientation2 = round(rad2deg(angle2Points(sensor,[hotspots_r(2),hotspots_c(2)])));
        %orientation2 = orientation2 + (orientation2==0)*360;
        
        %orientation =
        %round(wrapTo360(meanangle([orientation1,orientation2]))) % this has problem 
        %with [270,90] or similar cases
        %orientation = round(wrapTo360(rad2deg(circ_mean([deg2rad(orientation1),...
        %                                                 deg2rad(orientation2)],[],2))));
        orientation = orientation1 + (orientation1==0)*360;
        
        %orientation = orientation2; % for debugging
        
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
        
        oConfOrns(i) = orientation;
        
    end
end


for i = 1:numel(executedConfsHot_ind)
    
        
    [rr,cc] = ind2sub(size(map_candConf),executedConfsHot_ind(i));
    sensor  = [rr,cc];
        
    j = executedConfsHot_orn(i);
    
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
    V(visDom,executedConfsHot_ind(i)) = 1;
    
    oConfOrns(executedConfsHot_ind(i)) = j;
end
   

% sensing config can not visualize own cell
for i = 1:numel(map_candConf)    
    %V(i,(i-1)*f+1:(i-1)*f+f) = 0;
    V(i,i) = 0;
end

oV = V; % oV := visibility with obstacles

% VISIBILITY MATRIX CORRESPOND TO THE HOTSPOT ONLY:
% ------------------------------------------------------------------------------

% //////////////////////////
% ---- cells:
% //////////////////////////

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
cov_cell     = covered_cell_vector.*coverage_cell_vector.*free_cell_vector;
cellKept     = find(cov_cell>=1);
cellRemoved  = find(cov_cell==0);
V            = V(cellKept,:);


% ////////////////////////////
% -- cellKept weights
% ////////////////////////////
% cellKeptWeights = map_gd(cellKept);
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
plot(hotspots_r,hotspots_c,'or','MarkerFaceColor','r');

fixedConfs_ind = executedConfsHot_ind;
[fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),fixedConfs_ind);
plot(fixedConfs_r,fixedConfs_c,'or');

end




%----------------------------------------------------
% -- Check#2: confs. that cover all the hotspots
%----------------------------------------------------
conf_hotspot = ones(1,size(oV,2));
for i = 1:numel(hotspots)
    
    % conf ind that cover this hotspot
    conf_ind = find(oV(hotspots(i),:));
    % binary vector
    conf_this = zeros(1,size(oV,2));
    conf_this(conf_ind) = 1;
    
    % - update the overall vector
    conf_hotspot = conf_hotspot.*conf_this;
end
% - make sure fixed conf are in the vector
conf_hotspot(executedConfsHot_ind) = 1;

% -- plot 
if visualize == 1
figure('name','conf cover all the hotspots');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
conf_hotspot_ind = find(conf_hotspot(:)==1);
[conf_hotspot_r,conf_hotspot_c] = ind2sub(size(map_env),conf_hotspot_ind);
plot(conf_hotspot_r,conf_hotspot_c,'xb','MarkerFaceColor','b');        
plot(hotspots_r,hotspots_c,'or','MarkerFaceColor','r');

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
plot(hotspots_r,hotspots_c,'or','MarkerFaceColor','r');

fixedConfs_ind = executedConfsHot_ind;
[fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),fixedConfs_ind);
plot(fixedConfs_r,fixedConfs_c,'or');

end

% ----------------------------------------------------
% -- Check#4: conf not close enough to the hotspot 
% -----------------------------------------------------

% - cose enough radius
innerRadius_cell = para_.innerRadius_m/dataYAML.resolution;

% temporary inner limit map
map_candConf_minLimit = zeros(size(map_env));
map_candConf_minLimit(hotspots) = 1;

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
plot(hotspots_r,hotspots_c,'sr');

fixedConfs_ind = executedConfsHot_ind;
[fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),fixedConfs_ind);
plot(fixedConfs_r,fixedConfs_c,'or');

end

confNotClose = ones(1,size(oV,2));
confNotClose(innerLimit_ind) = 0;

% -- plot 
if visualize == 1
figure('name','conf Not Close');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])

[confNotClose_r,confNotClose_c] = ind2sub(size(map_env),find(confNotClose));
plot(confNotClose_r,confNotClose_c,'xb','MarkerFaceColor','b');
plot(hotspots_r,hotspots_c,'sr');

fixedConfs_ind = executedConfsHot_ind;
[fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),fixedConfs_ind);
plot(fixedConfs_r,fixedConfs_c,'or');

end

% ///////////////////////////////////////////////////////////////////////////
% Candidate conf: candidate conf map (Check#1) AND cover the hotspot
% (Check#2) AND are not placed over occupied cells (Check#3) AND are not
% close enough to the hotspots (Check#4)
% ///////////////////////////////////////////////////////////////////////////

conf_vector  = candConfNum.*conf_hotspot.*conf_free.*confNotClose;
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
plot(hotspots_r,hotspots_c,'or','MarkerFaceColor','r');

fixedConfs_ind = executedConfsHot_ind;
[fixedConfs_r,fixedConfs_c] = ind2sub(size(map_env),fixedConfs_ind);
plot(fixedConfs_r,fixedConfs_c,'or');
end

% **********************************************************
% -- Limit Number of Cand Conf based on Proportional Coverage
% **********************************************************



%/////////////////////////////////////// 
% -- Concentration weights
%///////////////////////////////////////

cellKeptWeights = map_recon(cellKept);
% cellKeptWeights = cellKeptWeights/(max(cellKeptWeights));

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
[hotspot_r,hotspot_c] = ind2sub(size(map_env),hotspots);
[conf_r,conf_c] = ind2sub(size(map_candConf),confKept);

desiredOrns = round(rad2deg(angle2Points([conf_r',conf_c'],[hotspot_r,hotspot_c])));
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
conf2Hot_dist = pdist2([conf_r',conf_c'],[hotspot_r,hotspot_c],'euclidean');

% -- distance gain
gain_sigma       = SensorRange.cell/5; % standard deviation
gain_meu         = SensorRange.cell/2; % expected value or center
V_dist_gain_this = gaussmf(conf2Hot_dist,[gain_sigma,gain_meu]);



% ------------------------------
% == orientation + displacement
% ------------------------------
wV_diplacement_this = (V_angular_gain_this+V_dist_gain_this)/2;

wV_this_a = ((gamma*(wV_concentration_this)))';
wV_this_b = ((1-gamma)* wV_diplacement_this);

wV_this = wV_this_a+wV_this_b;
%wV_this = (gamma*(wV_concentration_this))' + ((1-gamma)* wV_diplacement_this)



% -- sorted indices
[~,sorted_ind] = sort(wV_this,'descend');

sorted_ind_max = sorted_ind(1:min([para_.NumOfCandConf,numel(sorted_ind)]));
% sorted_ind_max = sorted_ind(1:min([10,numel(sorted_ind)]));
confToRemove_ind = find(ismember(sorted_ind,sorted_ind_max)==0);
confRemoved = union(confRemoved,sorted_ind(confToRemove_ind));
confKept = confKept(sorted_ind_max);

% -- sorted list
confKept = sort(confKept);

% -- plot 
if visualize == 1
figure('name','max num of conf');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
[confKept_r,confKept_c] = ind2sub(size(map_env),confKept);
plot(confKept_r,confKept_c,'xb','MarkerFaceColor','b');
plot(hotspots_r,hotspots_c,'sr');para_.numFoVs

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


% -- weighted visibility for fixedConf1 and fixedConf2
executedConfsHot_num = find(ismember(confKept,executedConfsHot_ind));




%//////////////////////////////////////////
% HOTSPOT PLACEMENT GAIN - ANGULAR
%//////////////////////////////////////////

% -- conf orientations
ConfOrns = oConfOrns(confKept);


% -- desired orinetations w.r.t to hotspot
% [hotspots_r,hotspots_c] = ind2sub(size(map_env),hotspots);
[conf_r,conf_c] = ind2sub(size(map_candConf),confKept);

desiredOrns = round(rad2deg(angle2Points([conf_r',conf_c'],[hotspots_r(1),hotspots_c(1)])));
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
conf2Hot_dist = pdist2([conf_r',conf_c'],[hotspots_r(1),hotspots_c(1)],'euclidean');
% -- distance gain
gain_sigma = SensorRange.cell/5; % standard deviation
gain_meu   = SensorRange.cell/2; % expected value or center
V_dist_gain = gaussmf(conf2Hot_dist,[gain_sigma,gain_meu]);


%//////////////////////////////////////////////
% HOTSPOT PLACEMENT GAIN - ANGULAR + DISTANCE
%//////////////////////////////////////////////
Conf_DistOrnGainV = (V_angular_gain+V_dist_gain)/2;




end




function [ visibleRange ] = fSenDomain( sensor,SensorRange,zawia )
%fSenDomain returns visible cells (domain/mask) for a given sensing
%configuration. 
%

% INTERPOLATING INTERMEDIATE ANGLES:
% -------------------------------------
DELTA_ANGLE = 1;
% --- Intermediate angles of visibile window. ---
% First half:
angl_1h = zawia(1):DELTA_ANGLE:zawia(2);
if zawia(1) > zawia(2) % if final angle is greater than 360 degrees.
    angl_1h = zawia(1):DELTA_ANGLE:360;
end
% Second half:
angl_2h = zawia(2):DELTA_ANGLE:zawia(3);
if zawia(2) > zawia(3) % if final angle is greater than 360 degrees.
    angl_2h = [zawia(2):DELTA_ANGLE:360 0:DELTA_ANGLE:zawia(3)];
end
angl = unique([angl_1h angl_2h]);


visibleRange      = [];        %
% subCellANGspecial = cell(1,4); % subject cells of special angles.

for i = 1:numel(angl)
    ALF = angl(i);
    % Find visible cells along one edge point.
    [ visCells ] = fRayTracSD(sensor,SensorRange,ALF);
    
    % Updating visible range.
    visibleRange = unique([visibleRange; visCells],'rows');   

end

end

function [ effective_visCells ] = fRayTracSD(sensor,SensorRange,ALF)
% fRayTracSD is ray tracing for sensor domain.

% CONSTANTS:
stepSize = 0.1; % step size of the ray length while sampling from origion to max. 

% initialize some high number of indices
effective_visCells = zeros(round(2*length(0.5:stepSize:SensorRange.cell)),2);
num = 0;
for d = [0.5:stepSize:SensorRange.cell SensorRange.cell]
    dr = 0.001*round((d*cosd(ALF))/0.001);  % row
    dc = 0.001*round((d*sind(ALF))/0.001);  % column
    
    vis_r = round(sensor(1)+round(dr)); % visible cell if validated
    vis_c = round(sensor(2)+round(dc)); % visible cell if validated
    
    if abs((dc-0.5)-fix(dc-0.5))>0.5 && abs((dr-0.5)-fix(dr-0.5))>0.5
        num = num+1;
        effective_visCells(num,:) = [vis_r vis_c];
    end
end
effective_visCells(num+1:end,:) = [];
effective_visCells = unique(effective_visCells,'rows');


end

function [ obsRange ] = fSenDomainOBS( sensor,SensorRange,zawia )
%fSenDomainOBS returns sensor's domain to check visibility against given
%obstacles. It contains more cells than the sensor's domain (mask).
%

% INTERPOLATING INTERMEDIATE ANGLES:
% ----------------------------------
DELTA_ANGLE = 1;
% --- Intermediate angles of visibile window. ---
% First half:
angl_1h = zawia(1):DELTA_ANGLE:zawia(2);
if zawia(1) > zawia(2) % if final angle is greater than 360 degrees.
    angl_1h = zawia(1):DELTA_ANGLE:360;
end
% Second half:
angl_2h = zawia(2):DELTA_ANGLE:zawia(3);
if zawia(2) > zawia(3) % if final angle is greater than 360 degrees.
    angl_2h = [zawia(2):DELTA_ANGLE:360 0:DELTA_ANGLE:zawia(3)];
end
angl = unique([angl_1h angl_2h]);

% initialize the obstacle range.
obsRange = []; 

% find the list of cells (for obstacles) in the sensing range.
for i = 1:numel(angl)
    
    % angle.
    ALF = angl(i);
    
    % Finding cells that are intersected along the ray.
    [ obsCells ] = fRayTracOD(sensor,SensorRange,ALF);
    
    % Updating the list of cells to be checked for obstacles.
    obsRange = unique([obsRange; obsCells],'rows');

end


end

function [ obsCells ] = fRayTracOD(sensor,SensorRange,ALF)
% fRayTracOD is ray tracing to find sensor domain for visibility check
% against listed obstacles

% point corresponding to the sensor range and ray angle.
subPoint_x = sensor(1) + (SensorRange.cell*cosd(ALF));
subPoint_y = sensor(2) + (SensorRange.cell*sind(ALF));

subPoint = [subPoint_x,subPoint_y];


% conditions.
cond_x1 = ALF>=0&&ALF<=45;
cond_x2 = ALF>=135&&ALF<=225;
cond_x3 = ALF>=315&&ALF<=360;

cond_y1 = ALF>=45&&ALF<=135;
cond_y2 = ALF>=225&&ALF<=315;

% sampling length.
delta = 0.01;
% Note: fine/coarse sampling changes the results for intersection points
% and hence changes the list of cells interseted. TODO: Find the optimal
% value for delta.


if cond_x1 || cond_x3
    
    % intermediate x and y points.
    x = [sensor(1):delta:subPoint(1) subPoint(1)];
    y = x*tand(ALF);
    
    % 4 decimal precision
    x = round(x*1e+4)/1e+4;
    y = round(y*1e+4)/1e+4;
    
    % corresponding cell numbers.
    cell_xy = round([x',y']);
    
    % cells that belong to the intersection point (x.5,y.5)
    ind_x  = find(abs(fix(x)-x)==0.5);
    ind_yy = abs(fix(y(ind_x))-y(ind_x))==0.5;
    ind_y  = ind_x(ind_yy);
    
    for i = 1:numel(ind_y)
        cell_xy = [cell_xy;...
                   ceil(x(ind_y(i))),  ceil(y(ind_y(i)))  ;...
                   ceil(x(ind_y(i))),  floor(y(ind_y(i))) ;...
                   floor(x(ind_y(i))), ceil(y(ind_y(i)))  ;...
                   floor(x(ind_y(i))), floor(y(ind_y(i)))];
    end    
    
elseif cond_x2
    
    % intermediate x and y points.    
    x = [sensor(1):-delta:subPoint(1) subPoint(1)];
    y = x*tand(ALF);
    
    % 4 decimal precision
    x = round(x*1e+4)/1e+4;
    y = round(y*1e+4)/1e+4;
    
    % corresponding cell numbers.
    cell_xy = round([x',y']);
    
    % cells that belong to the intersection point (x.5,y.5)
    ind_x  = find(abs(fix(x)-x)==0.5);
    ind_yy = find(abs(fix(y(ind_x))-y(ind_x))==0.5);
    ind_y  = ind_x(ind_yy);
    
    for i = 1:numel(ind_y)
        cell_xy = [cell_xy;...
                   ceil(x(ind_y(i))),  ceil(y(ind_y(i)))  ;...
                   ceil(x(ind_y(i))),  floor(y(ind_y(i))) ;...
                   floor(x(ind_y(i))), ceil(y(ind_y(i)))  ;...
                   floor(x(ind_y(i))), floor(y(ind_y(i)))];
    end
    
    
elseif cond_y1
    
    % intermediate x and y points.
    y = [sensor(2):delta:subPoint(2) subPoint(2)];
    x = y/tand(ALF);
    
    % 4 decimal precision
    x = round(x*1e+4)/1e+4;
    y = round(y*1e+4)/1e+4;
    
    % corresponding cell numbers.
    cell_xy = round([x',y']);
    
    % cells that belong to the intersection point (x.5,y.5)
    ind_x  = find(abs(fix(x)-x)==0.5);
    ind_yy = find(abs(fix(y(ind_x))-y(ind_x))==0.5);
    ind_y  = ind_x(ind_yy);
    
    for i = 1:numel(ind_y)
        cell_xy = [cell_xy;...
                   ceil(x(ind_y(i))),  ceil(y(ind_y(i)))  ;...
                   ceil(x(ind_y(i))),  floor(y(ind_y(i))) ;...
                   floor(x(ind_y(i))), ceil(y(ind_y(i)))  ;...
                   floor(x(ind_y(i))), floor(y(ind_y(i)))];
    end
    
    
elseif cond_y2
    
    % intermediate x and y points.
    y = [sensor(2):-delta:subPoint(2) subPoint(2)];
    x = y/tand(ALF);
    
    % 4 decimal precision
    x = round(x*1e+4)/1e+4;
    y = round(y*1e+4)/1e+4;
    
    % corresponding cell numbers.
    cell_xy = round([x',y']);
    
    % cells that belong to the intersection point (x.5,y.5)
    ind_x  = find(abs(fix(x)-x)==0.5);
    ind_yy = find(abs(fix(y(ind_x))-y(ind_x))==0.5);
    ind_y  = ind_x(ind_yy);
    
    for i = 1:numel(ind_y)
        cell_xy = [cell_xy;...
                   ceil(x(ind_y(i))),  ceil(y(ind_y(i)))  ;...
                   ceil(x(ind_y(i))),  floor(y(ind_y(i))) ;...
                   floor(x(ind_y(i))), ceil(y(ind_y(i)))  ;...
                   floor(x(ind_y(i))), floor(y(ind_y(i)))];
    end
    
end


 obsCells = unique(cell_xy,'rows');

end

function [ vDom ] = fVisSDvObs( sensor,sensorDom,obs )
%fVisibileArea returns visible cells for a given sensing configuration.
%

% initialize the list of visible cells.
vDom = zeros(size(sensorDom,1),1);

% for all cells in the sensor domain, check visibility against obs.
for i = 1:size(sensorDom,1)
    %i
    % cell to be validated for visibility
    subCell = [sensorDom(i,1) sensorDom(i,2)];

    % validatation for the visibility
    %[ cellValid ] = fCellValidSD( sensor,subCell,obstacDomOBS );
    [ validation ] = fCellValidSD( sensor,subCell,obs );

    % update visDom
    vDom(i) = validation;
    %pause
end

end


function [ validation ] = fCellValidSD( sensor,subCell,obs )
% 
% 

% angle to the center of the cell
ang = atan2(((subCell(2))-(sensor(2))),((subCell(1))-(sensor(1))))*(180/pi);
ang = ang+(360*(ang<0));

% conditions to project the line along (+/-x,+/-y).
cond_x1 = ang>=0&&ang<=45;
cond_x2 = ang>=135&&ang<=225;
cond_x3 = ang>=315&&ang<=360;

cond_y1 = ang>=45&&ang<=135;
cond_y2 = ang>=225&&ang<=315;

% sampling length.
delta = 0.01;
% Note: fine/coarse sampling changes the results for intersection points
% and hence changes the list of cells interseted. TODO: Find the optimal
% value for delta.

if cond_x1 || cond_x3
    
    % intermediate x and y points.
    x = sensor(1):delta:subCell(1);
    y = x*tand(ang);
    
    % 4 decimal precision
    x = round(x*1e+4)/1e+4;
    y = round(y*1e+4)/1e+4;
    
    % corresponding cell numbers.
    cell_xy = round([x',y']);
    
    % cells that belong to the intersection point (x.5,y.5)
    ind_x  = find(abs(fix(x)-x)==0.5);
    ind_yy = find(abs(fix(y(ind_x))-y(ind_x))==0.5);
    ind_y  = ind_x(ind_yy);
    
    for i = 1:numel(ind_y)
        cell_xy = [cell_xy;...
                   ceil(x(ind_y(i))),  ceil(y(ind_y(i)))  ;...
                   ceil(x(ind_y(i))),  floor(y(ind_y(i))) ;...
                   floor(x(ind_y(i))), ceil(y(ind_y(i)))  ;...
                   floor(x(ind_y(i))), floor(y(ind_y(i)))];
    end    
    
elseif cond_x2
    
    % intermediate x and y points.
    x = sensor(1):-delta:subCell(1);
    y = x*tand(ang);
    
    % 4 decimal precision
    x = round(x*1e+4)/1e+4;
    y = round(y*1e+4)/1e+4;
    
    % corresponding cell numbers.
    cell_xy = round([x',y']);
    
    % cells that belong to the intersection point (x.5,y.5)
    ind_x  = find(abs(fix(x)-x)==0.5);
    ind_yy = find(abs(fix(y(ind_x))-y(ind_x))==0.5);
    ind_y  = ind_x(ind_yy);
    
    for i = 1:numel(ind_y)
        cell_xy = [cell_xy;...
                   ceil(x(ind_y(i))),  ceil(y(ind_y(i)))  ;...
                   ceil(x(ind_y(i))),  floor(y(ind_y(i))) ;...
                   floor(x(ind_y(i))), ceil(y(ind_y(i)))  ;...
                   floor(x(ind_y(i))), floor(y(ind_y(i)))];
    end
    
    
elseif cond_y1
    
    % intermediate x and y points.
    y = sensor(2):delta:subCell(2);
    x = y/tand(ang);
    
    % 4 decimal precision
    x = round(x*1e+4)/1e+4;
    y = round(y*1e+4)/1e+4;
    
    % corresponding cell numbers.
    cell_xy = round([x',y']);
    
    % cells that belong to the intersection point (x.5,y.5)
    ind_x  = find(abs(fix(x)-x)==0.5);
    ind_yy = find(abs(fix(y(ind_x))-y(ind_x))==0.5);
    ind_y  = ind_x(ind_yy);
    
    for i = 1:numel(ind_y)
        cell_xy = [cell_xy;...
                   ceil(x(ind_y(i))),  ceil(y(ind_y(i)))  ;...
                   ceil(x(ind_y(i))),  floor(y(ind_y(i))) ;...
                   floor(x(ind_y(i))), ceil(y(ind_y(i)))  ;...
                   floor(x(ind_y(i))), floor(y(ind_y(i)))];
    end
    
    
elseif cond_y2
    
    % intermediate x and y points.
    y = sensor(2):-delta:subCell(2);
    x = y/tand(ang);
    
    % 4 decimal precision
    x = round(x*1e+4)/1e+4;
    y = round(y*1e+4)/1e+4;
    
    % corresponding cell numbers.
    cell_xy = round([x',y']);
    
    % cells that belong to the intersection point (x.5,y.5)
    ind_x  = find(abs(fix(x)-x)==0.5);
    ind_yy = find(abs(fix(y(ind_x))-y(ind_x))==0.5);
    ind_y  = ind_x(ind_yy);
    
    for i = 1:numel(ind_y)
        cell_xy = [cell_xy;...
                   ceil(x(ind_y(i))),  ceil(y(ind_y(i)))  ;...
                   ceil(x(ind_y(i))),  floor(y(ind_y(i))) ;...
                   floor(x(ind_y(i))), ceil(y(ind_y(i)))  ;...
                   floor(x(ind_y(i))), floor(y(ind_y(i)))];
    end
    
end

% remove repeated cells
cell_xy = unique(cell_xy,'rows');

% validation for the visibility
ind  = find(cell_xy(:,1)==obs(1)&cell_xy(:,2)==obs(2));
% if size(ind,1) > 0
%     validation = 0;
% end
validation = (size(ind,1) < 1);

end

function [ visDom ] = fVisibileDomain( map_coverage,map_env,sensorDom,obstacDom,VOBS )
%fVisibileDomain returns visible cells for a given sensing configuration.
%It takes the mask (sensor domain), finds the non-visible cells due to
%obstacles in the map and remove them from the mask.
%
% INPUTS:
% ------
% map       := [x,y] map
% sensorDom := [n,2] mask for the sensing configuration in question.
% obstacDom := [m,2] mask for the sensing configuration in question (used
%                    to check visibility against obstacles).
% VOBS      := [m,n] effect on visibility (columns) against each obstacle 
%                    (rows)
% OUTPUTS:
% ------
% visDom    := [k] visibile cells for the configuration in question 
%              (subscripts)

% rows and columns for the map
[map_row,map_col] = size(map_coverage);

% temporary domain/masks to handle intermediate computation quickly.
visDomX = sensorDom;
obsDomX = obstacDom;

% preRemoveList:= cells outside the map (sensor domain).
% along x-axis:
   i = find(visDomX(:,1)<=0);      % less than row number
   j = find(visDomX(:,1)>map_row); % greater than row number
   preRemoveList = [i;j];
% along y-axis:
   i = find(visDomX(:,2)<=0);      % less than col number
   j = find(visDomX(:,2)>map_col); % greater than col number
   preRemoveList = [preRemoveList;i;j];

% remove preRemoveList from temporary visible domain/mask.
% visDomX(preRemoveList,:) = [];

% preRemoveListOD:= cells outside the map (obstacle domain).
% along x-axis:
   i = find(obsDomX(:,1)<=0);      % less than row number
   j = find(obsDomX(:,1)>map_row); % greater than row number
   preRemoveListOD = [i;j];
% along y-axis:
   i = find(obsDomX(:,2)<=0);      % less than col number
   j = find(obsDomX(:,2)>map_col); % greater than col number
   preRemoveListOD = [preRemoveListOD;i;j];

% remove preRemoveList from temporary visible domain/mask.
obsDomX(preRemoveListOD,:) = [];


% convert visDomX from ind to sub.
% visDomX_ind = sub2ind(size(map),visDomX(:,1),visDomX(:,2));
obsDomX_ind = sub2ind(size(map_coverage),obsDomX(:,1),obsDomX(:,2));

% find the list of obstacles in visDomX;
% obs_list = find(map(visDomX_ind)==0);
obs_listOD = find(map_env(obsDomX_ind)==0);

% VOBSx:= VOBS with outside cells are removed.
% VOBSx = VOBS; 
% VOBSx(preRemoveList,:) = [];

VOBSx = VOBS; 
VOBSx(preRemoveListOD,:) = [];

% non-visible cells due to obstacles (logical values).
% nVisCells = (VOBSx(obs_list',:)==0);
nVisCells = (VOBSx(obs_listOD',:)==0);

% convert logical values to a list.
list = find(nVisCells'==1);

% translate the list to the column numbers of VOBS, which is remove list.
[removeList,~] = ind2sub(size(VOBSx'),list);

% update the RemoveList with preRemoveList and removeList.
RemoveList = [preRemoveList;removeList];

% finally, we know non-visible cells in the domain/mask.
visDom = sensorDom;
visDom(RemoveList,:) = [];

% convert indices into subscripts and sort them.
if visDom  
    visDom = sort(sub2ind(size(map_coverage),visDom(:,1),visDom(:,2)));
end

end

% ---------------------------------------------------- End of the document -------------------------------------------------------

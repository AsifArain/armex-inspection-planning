function [ conf_theta,...
           conf_crossAngles,...
           conf_crossAngles_G,...
           SensorRange,...
           OBS,...
           V,...
           oV,...
           confKept,...
           confRemoved,...
           cellKept,...
           cellRemoved,...
           FoVfaces,...
           map_candConf,...
           map_coverage,...
           map_env,... 
           BestFocusPointGain,...
           BestFocusPointAng,...
           refDistPoint,...
           D,...
           tPreprocess ] = fPreprocessHotspotVCrossAngles( FoV,...
                                                           OBS,...
                                                           SensorRange,...
                                                           map_coverage,...
                                                           map_env,...
                                                           hotspot,...
                                                           FocusTSPAng,...
                                                           dataYAML,...
                                                           para_,...
                                                           visualize )
% fPreProcess_R4 genereates required parameters for solving a global coverage problem.
% Date: 2014-01-19, Rev 4
% 
% n  := total number of cells in the map.
% nf := number of free cells in the map.
% f  := total number of sensing configurations per cell.
% 
% OUTPUTS:
% --------------------------------------------------------------------------------------------------------------------------------
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


%% SENSING PARAMETERS AND CONSTANTS:
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
FoVfaces.num        = para_.numFoVs;   % number of conf per cell.
n                   = numel(map_env);         % number of cells in the map.
f                   = FoVfaces.num;       % for simplification.
c                   = n*f;                % total conf.
cnt                 = para_.GraphCnt;  % 4 or 8 connected neighbours.

% o := number of outgoing edges for each configuration.
if cnt == 4
    if f == 1
        o = 4;
        oO = [1 1 1 1;...    % T - 1st 1
              0 0 0 0;...    % T - 2nd 1.4
              0 0 0 0;...    % T - 3rd inf
              0 0 0 0;...    % R - 1st 0.1
              0 0 0 0];      % R - 2nd 0.2          
         eO = oO;
          
    elseif f == 4
        o = 3;
        oO = [1 0 0;...    % T - 1st 1
              0 0 0;...    % T - 2nd 1.4
              0 0 0;...    % T - 3rd inf
              0 0 0;...    % R - 1st 0.1
              0 1 1];      % R - 2nd 0.2
        eO = oO;        
    end   
    
elseif cnt == 8
    if f == 1
        o = 8;
        %     1 2 3 4 5 6 7 8 -     IDs of outgoing edges
        oO = [1 0 1 0 1 0 1 0;...    % T - 1st 1
              0 1 0 1 0 1 0 1;...    % T - 2nd 1.4
              0 0 0 0 0 0 0 0;...    % T - 3rd inf
              0 0 0 0 0 0 0 0;...    % R - 1st 0.1
              0 0 0 0 0 0 0 0];      % R - 2nd 0.2
        eO = oO;
        
    elseif f == 8
        
        o = 3;
        oO = [1 0 0;...    % T - 1st 1
              0 0 0;...    % T - 2nd 1.4
              0 0 0;...    % T - 3rd inf
              0 1 1;...    % R - 1st 0.1
              0 0 0];      % R - 2nd 0.2
          
        eO = [0 0 0;...    % T - 1st 1
              1 0 0;...    % T - 2nd 1.4
              0 0 0;...    % T - 3rd inf
              0 1 1;...    % R - 1st 0.1
              0 0 0];      % R - 2nd 0.2
    end
end
e = c*o; % # of edges



%**********************************************************************
% CANDIDATE CONF MAP BASED ON SENSING RANGE:
%**********************************************************************

fprintf(1,'\n::::: candidate conf based on sensing range \n')

[hotspot_r,hotspot_c] = ind2sub(size(map_env),hotspot);

[ map_candConf ] = fConfMapSensingRange( hotspot,...
                                         map_env,...
                                         SensorRange,...
                                         OBS,...
                                         visualize );
                                     
% %% CANDIDATE CONF MAP (1ST2ND) LATEST VERSION:
% % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 
% [hotspot_r,hotspot_c] = ind2sub(size(map_env),hotspot);
% 
% % initialize map for candidate conf. -- ring only
% map_candConf_1st2nd = zeros(size(map_env));    
% map_candConf_1st2nd(hotspot) = 1;
% 
% % max radius to grow hotspot
% % candCircleThickness_m = 0 % thickness of candidate conf circle (meters)
% % candCircleThickness_cell = candCircleThickness_m/paramPreprocess.MapRes
% % radius_candCircleOut = (SensorRange.cell/4)+(candCircleThickness_cell/2)
% radius_candCircleOut = ceil(SensorRange.cell/4)
% % radius_candCircleOut = round(SensorRange.cell/4)
% 
% % grow hotspot to max radius
% se = strel('disk',radius_candCircleOut,0); % outer limit
% map_candConf_1st2nd = imdilate(map_candConf_1st2nd,se);
% 
% % plot
% figure('name','cand conf - max limit and obs are removed');
% imshow(map_env','InitialMagnification',800); hold on;
% set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
% for i = 1:size(map_candConf_1st2nd,1)
%     for j = 1:size(map_candConf_1st2nd,2)
%         if map_candConf_1st2nd(i,j)
%             %plot(j,i,'xb','MarkerFaceColor','b');
%             plot(i,j,'xb','MarkerFaceColor','b');
%         end
%     end
% end
% % plot(hotspot_c,hotspot_r,'xr','MarkerFaceColor','r');
% plot(hotspot_r,hotspot_c,'xr','MarkerFaceColor','r');
% 
% % pause
% 
% % remove obs cells
% map_candConf_1st2nd(OBS.ind(find(ismember(OBS.ind,find(map_candConf_1st2nd))))) = 0;
% 
% % plot
% figure('name','cand conf - max limit and obs are removed');
% imshow(map_env','InitialMagnification',800); hold on;
% set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
% for i = 1:size(map_candConf_1st2nd,1)
%     for j = 1:size(map_candConf_1st2nd,2)
%         if map_candConf_1st2nd(i,j)
%             %plot(j,i,'xb','MarkerFaceColor','b');
%             plot(i,j,'xb','MarkerFaceColor','b');
%         end
%     end
% end
% % plot(hotspot_c,hotspot_r,'xr','MarkerFaceColor','r');
% plot(hotspot_r,hotspot_c,'xr','MarkerFaceColor','r');
% 
% % extract edges
% map_candConf_1st2nd = bwmorph(map_candConf_1st2nd,'remove');
% 
% % plot
% pause(0.1)
% figure('name','cand conf - extracted edges'); 
% imshow(map_env','InitialMagnification',800); hold on;
% set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
% for i = 1:size(map_candConf_1st2nd,1)
%     for j = 1:size(map_candConf_1st2nd,2)
%         if map_candConf_1st2nd(i,j)
%             %plot(j,i,'xb','MarkerFaceColor','b');
%             plot(i,j,'xb','MarkerFaceColor','b');
%         end
%     end
% end
% %plot(hotspot_c,hotspot_r,'xr','MarkerFaceColor','r');
% plot(hotspot_r,hotspot_c,'xr','MarkerFaceColor','r');
% 
% % grow edges with prdefiend thickness 
% % circleThickness_m = 2.5
% % circleThickness_cell = circleThickness_m/paramPreprocess.MapRes
% circleThickness_cell = ceil(SensorRange.cell/4);
% se = strel('ball',circleThickness_cell,0); % outer limit
% 
% map_candConf_1st2nd = imdilate(map_candConf_1st2nd,se);
% 
% 
% % plot
% pause(0.1)
% figure('name','cand conf - edges are increased to the limit');
% imshow(map_env','InitialMagnification',800); hold on;
% set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
% for i = 1:size(map_candConf_1st2nd,1)
%     for j = 1:size(map_candConf_1st2nd,2)
%         if map_candConf_1st2nd(i,j)
%             %plot(j,i,'xb','MarkerFaceColor','b');
%             plot(i,j,'xb','MarkerFaceColor','b');
%         end
%     end
% end
% % plot(hotspot_c,hotspot_r,'xr','MarkerFaceColor','r');
% plot(hotspot_r,hotspot_c,'xr','MarkerFaceColor','r');
% 
% % remove obstacle 
% %map_candConf_1st2nd(OBS.ind(find(ismember(OBS.ind,find(map_candConf_1st2nd))))) = 0;
% 
% % ---- lets have new obs list by growing the obstacles one meter
% thikness_cells = ceil(1.5/dataYAML.resolution)
% se = strel('square',thikness_cells); % outer limit
% map_grown_obs_list = ~map_env;
% map_grown_obs_list = imdilate(map_grown_obs_list,se);
% map_grown_obs_list = ~map_grown_obs_list;
% 
% 
% pause(0.1)
% figure('name','map_grown_obs_list'); 
% imshow((map_env-map_grown_obs_list)','InitialMagnification',800);
% set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
% 
% 
% obs_list_grown = find(~map_grown_obs_list);
% map_candConf_1st2nd(obs_list_grown(find(ismember(obs_list_grown,find(map_candConf_1st2nd))))) = 0;
% 
% % plot
% pause(0.1)
% figure('name','cand conf - occupied cells are removed');
% imshow(map_env','InitialMagnification',800); hold on;
% set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
% for i = 1:size(map_candConf_1st2nd,1)
%     for j = 1:size(map_candConf_1st2nd,2)
%         if map_candConf_1st2nd(i,j)
%             %plot(j,i,'xb','MarkerFaceColor','b');
%             plot(i,j,'xb','MarkerFaceColor','b');
%         end
%     end
% end
% % plot(hotspot_c,hotspot_r,'xr','MarkerFaceColor','r');
% plot(hotspot_r,hotspot_c,'xr','MarkerFaceColor','r');
% % figure; imshow(map_candConf_1st2nd)
% 
% 
% % ---- now, no need to remove cand conf from inner part. it is implemented
% % in visibility function
% %{
% % min radius where sensing conf are not allowed.
% % innerRadius_m = 2
% % paramPreprocess.innerRadius_m = 2
% innerRadius_cell = paramPreprocess.innerRadius_m/paramPreprocess.MapRes;
% 
% % temporary inner limit map
% map_candConf_minLimit = zeros(size(map_env));
% map_candConf_minLimit(hotspot) = 1;
% 
% se = strel('disk',innerRadius_cell,0);
% map_candConf_minLimit = imdilate(map_candConf_minLimit,se);
% % figure; imshow(map_candConf_minLimit)
% 
% % indices of inner limit map
% % innerLimit_ind = find(map_candConf_minLimit);
% 
% % remove cand conf from inner limit
% map_candConf_1st2nd(find(map_candConf_minLimit)) = 0;
% 
% % plot
% figure('name','cand conf - inner radius is removed'); imshow(map_env); hold on;
% for i = 1:size(map_candConf_1st2nd,1)
%     for j = 1:size(map_candConf_1st2nd,2)
%         if map_candConf_1st2nd(i,j)
%             plot(j,i,'xb','MarkerFaceColor','b');
%         end
%     end
% end
% plot(hotspot_c,hotspot_r,'xr','MarkerFaceColor','r');
% % figure; imshow(map_candConf_1st2nd)
% 
% numel(find(map_candConf_1st2nd))
% %}
% 
% % ---------------------
% 

%% VISIBILITY MATRIX:
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

disp('V...')
[ V,...
  oV,...
  confKept,...
  confRemoved,...
  cellKept,...
  cellRemoved,...
  FoVfaces ] = fV( map_candConf,...
                   map_coverage,...
                   map_env,...
                   n,...
                   f,...
                   c,...
                   FoV,...
                   FoVfaces,...
                   SensorRange,...
                   OBS,...
                   hotspot,...
                   dataYAML,...
                   para_ );




%% CANDIDATE CONF MAP (FINAL) FROM 1ST2ND MAP:
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% cand conf that do not cover hotspot are removed from the second map.

% initialize final candidate conf map
map_candConf = zeros(size(map_env));
% corresponding cell numbers for conf kept.
confCELL = fix(confKept/f)+(~~mod(confKept,f)); % cell num for each conf num.
% update the map
map_candConf(confCELL) = 1;
% figure; imshow(map_candConf)
% numel(find(map_candConf))*4




% plot on environment map
if visualize == 1
figure('name','cand conf - final');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
% [candCell_r,candCell_c] = find(map_candConf);
for i = 1:size(map_candConf,1)
    for j = 1:size(map_candConf,2)        
        if map_candConf(i,j)
            %plot(j,i,'xb','MarkerFaceColor','b');
            plot(i,j,'xb','MarkerFaceColor','b');
        end
    end
end
% plot(hotspot_c,hotspot_r,'xr','MarkerFaceColor','r');
plot(hotspot_r,hotspot_c,'xr','MarkerFaceColor','r');
end




%% CROSS ANGLES:
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

disp('Cross angles...')
[ conf_theta,...
  conf_crossAngles,...
  conf_crossAngles_G] = fCrossAngles( oV,map_env,map_candConf,f,confKept,hotspot,para_ );


%% REFERENCE POINT FOR DISTANCE:
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

disp('Reference point for distance...')

[map_row,map_col] = size(map_candConf);

% cell numbers of all valid candidate conf
AllCandConfCell = fix(confKept/f)+(~~mod(confKept,f)); % cell num for each conf num.

% --- light angles of possible areas for distance point
% LightAngles = [0,90,180,270]
% deltaAng = 10;
% paramPreprocess.AngularStepForDistanceFocusPoint = 10;
LightAngles = 0:para_.AngularStepForDistanceFocusPoint:360-para_.AngularStepForDistanceFocusPoint;
LightAngles = LightAngles+FocusTSPAng;
LightAngles = LightAngles+(360*(LightAngles<0))-(360*(LightAngles>360));

% --- cand conf domain for each focus point
sensor = [hotspot_r,hotspot_c];
% DomainRadius.cell = radius_candCircleOut;
DomainRadius.cell = SensorRange.cell/2;

% --- initialize variables
AllFcousPoints_NumOfCandConfs = zeros(numel(LightAngles),1);
AllFcousPoints_NumOfFreeCells = zeros(numel(LightAngles),1);
AllFcousPoints_CircmCandConfs = zeros(numel(LightAngles),1);
for i = 1:numel(LightAngles)
    
    % -- domain of focus point
    zawia = [LightAngles(i)-90,LightAngles(i),LightAngles(i)+90];
    zawia+(360*(zawia<0))-(360*(zawia>360));
    [ FocusPointDomain ] = fSenDomain( sensor,DomainRadius,zawia );
    
    % -- finding domain cells outside map
    idxR1 = find(FocusPointDomain(:,1)<=0);      % less than row number
    idxR2 = find(FocusPointDomain(:,1)>map_row); % greater than row number    
    idxC1 = find(FocusPointDomain(:,2)<=0);      % less than col number
    idxC2 = find(FocusPointDomain(:,2)>map_col); % greater than col number
    
    invalidList = [idxR1;idxR2;idxC1;idxC2];
    invalidList = unique(invalidList);
    
    % -- valid domain
    FocusPointDomainValid = FocusPointDomain;
    FocusPointDomainValid(invalidList,:) = [];        
    FocusPointDomainValid_ind = sub2ind(size(map_candConf),...
                                        FocusPointDomainValid(:,1),FocusPointDomainValid(:,2));
    
    % -- number of valid candidate conf in the domain of focus point
    FocusPointCandConfs = numel(find(ismember(AllCandConfCell,FocusPointDomainValid_ind)));
    
    % -- number of obstacles in the domain
    %FocusPointObstacles = numel(find(ismember(OBS.ind,FocusPointDomainValid_ind)));
    %FocusPointObstacles = FocusPointObstacles+numel(invalidList);
    % instead, lets find number of unoccupied cells    
    FocusPointFreeCells = size(FocusPointDomain,1)-...
                          numel(find(ismember(OBS.ind,FocusPointDomainValid_ind)))-...
                          numel(invalidList);
    
    % -- angles between valid domain confs
    CandConfInDomain      = find(ismember(AllCandConfCell,FocusPointDomainValid_ind));
    CrossAnglesDomainConf = conf_crossAngles(CandConfInDomain,CandConfInDomain);
    FocusPointConfCircm   = sum(sum(CrossAnglesDomainConf*1,1),2); % cirum = radians*radius
    
    % -- combined results for all focus points
    AllFcousPoints_NumOfCandConfs(i) = FocusPointCandConfs;
    AllFcousPoints_NumOfFreeCells(i) = FocusPointFreeCells;
    AllFcousPoints_CircmCandConfs(i) = FocusPointConfCircm;
    
    
    
end

% normalized - number of cand conf, number of free cells and circumference
AllFcousPoints_NumOfCandConfsNorm = AllFcousPoints_NumOfCandConfs/max(AllFcousPoints_NumOfCandConfs);
AllFcousPoints_NumOfFreeCellsNorm = AllFcousPoints_NumOfFreeCells/max(AllFcousPoints_NumOfFreeCells);
AllFcousPoints_CircmCandConfsNorm = AllFcousPoints_CircmCandConfs/max(AllFcousPoints_CircmCandConfs);


% overall gain
AllFcousPoints_Gain = AllFcousPoints_NumOfCandConfsNorm.*...
                      AllFcousPoints_NumOfFreeCellsNorm.*...
                      AllFcousPoints_CircmCandConfsNorm;
% best region
[BestFocusPointGain,BestFocusPointInd] = max(AllFcousPoints_Gain);
BestFocusPointAng = LightAngles(BestFocusPointInd);



for i = 1:numel(LightAngles)
    %LightAngles(i)
    % point corresponding to the sensor range and ray angle.
    %endPoint_r = hotspot_r + (DomainRadius.cell*cosd(LightAngles(i)));
    %endPoint_c = hotspot_c + (DomainRadius.cell*sind(LightAngles(i)));
    endPoint_r = hotspot_r + (SensorRange.cell/4*cosd(LightAngles(i)));
    endPoint_c = hotspot_c + (SensorRange.cell/4*sind(LightAngles(i)));
    %plot(endPoint_c,endPoint_r,'+g');
    %plot([hotspot_c,endPoint_c],[hotspot_r,endPoint_r],'-g'); 
    
    if visualize == 1
    plot(endPoint_r,endPoint_c,'+g');
    plot([hotspot_r,endPoint_r],[hotspot_c,endPoint_c],'-g'); 
    end
end

% refDistPoint_r = hotspot_r + (DomainRadius.cell*cosd(BestFocusPointAng))
% refDistPoint_c = hotspot_c + (DomainRadius.cell*sind(BestFocusPointAng))
% 
% refDistPoint = [refDistPoint_r,refDistPoint_c]

% --- lets go back to TSP angles
% refDistPoint_r = hotspot_r + (DomainRadius.cell*cosd(FocusTSPAng))
% refDistPoint_c = hotspot_c + (DomainRadius.cell*sind(FocusTSPAng))
refDistPoint_r = hotspot_r + (SensorRange.cell/4*cosd(FocusTSPAng));
refDistPoint_c = hotspot_c + (SensorRange.cell/4*sind(FocusTSPAng));

refDistPoint = [refDistPoint_r,refDistPoint_c];

if visualize == 1
% plot(refDistPoint_c,refDistPoint_r,'*r');
plot(refDistPoint_r,refDistPoint_c,'*r');
end

%% DISTANCE MATRIX:
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
disp('Distance matrix...')

% ---- D (distance) matrix is manhattan distance between cand conf position
% and reference position
confKeptCells_ind = fix(confKept/f)+(~~mod(confKept,f));
[confKeptCells_r,confKeptCells_c] = ind2sub(size(map_env),confKeptCells_ind);

D = zeros(numel(confKeptCells_ind),1);
for i = 1:numel(confKeptCells_ind)
    manhattan_distance = abs(confKeptCells_r(i)-refDistPoint_r)+abs(confKeptCells_c(i)-refDistPoint_c);
    D(i) = manhattan_distance;
end


%% Total computation time.
% ---------------------------------------------------------------------------------------------------------------------------
tPreprocess =  1e-4*(round(toc(tPreprocessi)*1e4));

disp(['tPreprocess: ',num2str(tPreprocess),' sec']);


end




function [ map_candConf ] = fConfMapSensingRange( hotspot,...
                                                  map_env,...
                                                  SensorRange,...
                                                  OBS,...
                                                  visualize )

[hotspot_r,hotspot_c] = ind2sub(size(map_env),hotspot);

% initialize map for candidate conf. -- ring only
map_candConf = zeros(size(map_env));    
map_candConf(hotspot) = 1;

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
plot(hotspot_r,hotspot_c,'xr','MarkerFaceColor','r');
end

% -- remove obs cells
map_candConf(OBS.ind(ismember(OBS.ind,find(map_candConf)))) = 0;

% plot
if visualize == 1
figure('name','cand conf - max limit and obs are removed');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
map_candConf_ind = find(map_candConf(:)==1);
[map_candConf_r,map_candConf_c] = ind2sub(size(map_env),map_candConf_ind);
plot(map_candConf_r,map_candConf_c,'xb','MarkerFaceColor','b');        
plot(hotspot_r,hotspot_c,'xr','MarkerFaceColor','r');
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
plot(hotspot_r,hotspot_c,'xr','MarkerFaceColor','r');
end



end




function [conf_theta,...
          conf_crossAngles,...
          conf_crossAngles_G] =...
                    fCrossAngles( oV,map_env,map_candConf,f,confKept,hotspot,para_ )

% conf_theta            := conf orientations [1,c]
% conf_crossAngles      := cross angles between conf [c,c]
% conf_crossAngles_G    := cross angles gain [c,c]

%{
center_of_mass = paramPreprocess.CoM; % reference cell
[center_of_mass_y,center_of_mass_x] = ind2sub(size(map_candConf),center_of_mass);

% ----- conf orientations
conf_theta = 1000*ones(1,size(oV,2)); % vector contains orientations of candidate conf.

for conf = 1:size(oV,2)
    conf_cell_num = fix(conf/f)+(~~mod(conf,f));   % cell num for each conf num.
    [conf_cell_y,conf_cell_x] = ind2sub(size(map_candConf),conf_cell_num); % row and col of cell num.
    
    if map_candConf(conf_cell_y,conf_cell_x)
        deltaX = conf_cell_x - center_of_mass_x;
        deltaY = conf_cell_y - center_of_mass_y;
        conf_theta(conf) = atan2(deltaY,deltaX);
    end
    
end
%}

center_of_mass = para_.CoM; % reference cell
[center_of_mass_r,center_of_mass_c] = ind2sub(size(map_candConf),center_of_mass);

% ----- conf orientations
conf_theta = 1000*ones(1,size(oV,2)); % vector contains orientations of candidate conf.

for conf = 1:size(oV,2)
    conf_cell_num = fix(conf/f)+(~~mod(conf,f));   % cell num for each conf num.
    [conf_cell_r,conf_cell_c] = ind2sub(size(map_candConf),conf_cell_num); % row and col of cell num.
    
    if map_candConf(conf_cell_r,conf_cell_c)
        deltaR = conf_cell_r - center_of_mass_r;
        deltaC = conf_cell_c - center_of_mass_c;
        conf_theta(conf) = atan2(deltaC,deltaR);
    end
    
end



% ----- cross angles
disp('----- cross angles')
% conf_crossAngles = inf(size(oV,2)); % matrix of cross angles between candidate conf.
conf_crossAngles = sparse(size(oV,2)); % matrix of cross angles between candidate conf.

% for conf_first = 1:size(oV,2)
%     for conf_second = 1:size(oV,2)
%         if map_candConf(fix(conf_first/f)+(~~mod(conf_first,f)))
%             if map_candConf(fix(conf_second/f)+(~~mod(conf_second,f)))
%                 
%                 conf_crossAngles(conf_first,conf_second) =...
%                     angleAbsDiff(conf_theta(conf_first),conf_theta(conf_second));
%             end
%         end
%     end
% end

% plot for debugging
% figure('name','cross angles debug'); imshow(map_env); hold on;
% [hotspot_r,hotspot_c] = ind2sub(size(map_candConf),hotspot);
% plot(hotspot_c,hotspot_r,'xr','MarkerFaceColor','r');

for i = 1:numel(confKept)
    for j = 1:numel(confKept)
        
        conf_crossAngles(confKept(i),confKept(j)) =...
            angleAbsDiff(conf_theta(confKept(i)),conf_theta(confKept(j)));
        
        % -- for debugging
        %{
        cross_angles_degree = rad2deg(angleAbsDiff(conf_theta(confKept(i)),conf_theta(confKept(j))))
        
        conf1_cell_ind = fix(confKept(i)/f)+(~~mod(confKept(i),f));
        conf2_cell_ind = fix(confKept(j)/f)+(~~mod(confKept(j),f));
        
        [conf1_r,conf1_c] = ind2sub(size(map_env),conf1_cell_ind);
        [conf2_r,conf2_c] = ind2sub(size(map_env),conf2_cell_ind);
        
        [conf1_r,conf1_c]
        [conf2_r,conf2_c]
        
        h1 = plot(conf1_c,conf1_r,'+b','MarkerFaceColor','b');
        h2 = plot(conf2_c,conf2_r,'+g','MarkerFaceColor','g');
        
        pause;
        delete(h1,h2);
        clc;
        %}
        
    end
end

% ------- removing conf sampled over occupied cells
%{
disp('----- removing conf sampled over occupied cells')

% obs-conf nums
confOBS = zeros(numel(OBS_candConf.ind)*f,1);
for i = 1:numel(OBS_candConf.ind)
    confOBS(((i-1)*f)+1:((i-1)*f)+f) = (((OBS_candConf.ind(i)-1)*f)+1):(((OBS_candConf.ind(i)-1)*f)+f);
end

conf_theta(:,confOBS)  = [];
conf_crossAngles(:,confOBS)        = [];
% conf_cross_angles_gain(:,confOBS)   = [];

conf_crossAngles(confOBS,:)        = [];
% conf_cross_angles_gain(confOBS,:)   = [];
%}

% conf variables corresponding to candidate conf map only:
% ------------------------------------------------------------------------------
%{
% --- selecting conf number corresponding to candidate conf map only
% cells for the candidate conf.
candConfCells = find(map_candConf); 
% conf numbers
candConfNum = zeros(numel(map_candConf)*f,1);
for i = 1:numel(candConfCells)
    candConfNum( (candConfCells(i)-1)*f+1:(candConfCells(i)-1)*f+f ) = 1;
end

% update conf orientation vector
conf_theta = conf_theta(find(candConfNum));
% update conf cross angles matrix
conf_crossAngles = conf_crossAngles(:,find(candConfNum));
conf_crossAngles = conf_crossAngles(find(candConfNum),:);
%}

% --- instead, use info from fV:
% update conf orientation vector
conf_theta = conf_theta(confKept);
% update conf cross angles matrix
conf_crossAngles = conf_crossAngles(:,confKept);
conf_crossAngles = conf_crossAngles(confKept,:);
conf_crossAngles = full(conf_crossAngles);


% ------ cross angles gain
disp('----- gain')
gain_sigma = deg2rad(para_.xGainSigma_DEG); % standard deviation
gain_meu   = deg2rad(para_.xGainMeu1_DEG); % expected value or center
conf_crossAngles_G1 = gaussmf(conf_crossAngles,[gain_sigma,gain_meu]);


gain_sigma = deg2rad(para_.xGainSigma_DEG); % standard deviation
gain_meu   = deg2rad(para_.xGainMeu2_DEG); % expected value or center
conf_crossAngles_G2 = gaussmf(conf_crossAngles,[gain_sigma,gain_meu]);

conf_crossAngles_G = conf_crossAngles_G1+conf_crossAngles_G2;
 


% ------ publish gussian gain function for test only
%{ 
angles = 0:1:180; % in degrees
gain_sigma = 10; % standard deviation
% gain-1
gain_meu   = 60; % expected value or center
guss_gain1 = gaussmf(angles,[gain_sigma,gain_meu]);
% gain-2
gain_meu   = 120;% expected value or center
guss_gain2 = gaussmf(angles,[gain_sigma,gain_meu]);
% final
guss_gain = guss_gain1+guss_gain2;
% plot
figure; plot(angles,guss_gain1,'-r')
figure; plot(angles,guss_gain2,'-g')
figure; plot(angles,guss_gain ,'-b')


hXLabel = xlabel('Cross angles (degree)');
hYLabel = ylabel('Gain');
xlim([0,180])
ylim([0,1.2])
set([hXLabel,hYLabel],'Interpreter','latex')
set(gca, ...
      'Box'         , 'off'     , ...
      'TickDir'     , 'out'     , ...
      'TickLength'  , [.01 .01] , ...
      'XMinorTick'  , 'off'     , ...
      'YMinorTick'  , 'off'      , ...
      'YGrid'       , 'off'      , ...
      'XColor'      , [.3 .3 .3], ...
      'YColor'      , [.3 .3 .3], ...
      'XTick'       , (0:20:180),   ...  
      'YTick'       , (0:.2:1.2),   ...
      'LineWidth'   , 1         );
set([gca,hXLabel,hYLabel],'FontSize',14);
set(gcf,'PaperPositionMode','auto','PaperSize',[5.85 4.7]);
eval(['matlab2tikz(''fig/cross_angles_gain.tikz'',''height'',''\figureheight'',''width'',''\figurewidth'',''extraAxisOptions'',''xlabel near ticks'',''extraAxisOptions'',''ylabel near ticks'')']);

%}


end

function [ V,...
           oV,...
           confKept,...
           confRemoved,...
           cellKept,...
           cellRemoved,...
           FoVfaces ] = fV( map_candConf,...
                            map_coverage,...
                            map_env,...
                            n,...
                            f,...
                            c,...
                            FoV,...
                            FoVfaces,...
                            SensorRange,...
                            OBS,...
                            hotspot,...
                            dataYAML,...
                            para_ )
%fV generates Visibility matrix V.

% --- ANGLES:
FoVfaces.lgt = [1/FoVfaces.num:1/FoVfaces.num:1-(1/FoVfaces.num) 0]'*360;
FoVfaces.alf = FoVfaces.lgt-(FoV/2);
FoVfaces.bay = FoVfaces.lgt+(FoV/2);
% compensate negative alf
FoVfaces.alf(FoVfaces.alf<0) = 360+FoVfaces.alf(FoVfaces.alf<0);
FoVfaces.ang = [FoVfaces.alf FoVfaces.bay];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SensorDomain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('.. Sensor Domain')

switch para_.SenDom
    
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
        save(['preprocess/SenDom_Rng',num2str(SensorRange.cell),'FoV',num2str(FoV)],...
            'sD','oD','VObs');
        
    case 'UseEarlier'
        
        load(['preprocess/SenDom_Rng',num2str(SensorRange.cell),'FoV',num2str(FoV)]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('.. Visibility Matrix')
V = sparse(n,c); % initializing
% V = zeros(n,c); % initializing

% For each cell of the map.
for i = 1:n
    if map_candConf(i) % if its not occupied
        [rr,cc] = ind2sub(size(map_candConf),i);
        sensor = [rr,cc];
        
        % For each configuration of a cell.
        for j = 1:f
            
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
            V(visDom,((i-1)*f)+j) = 1;
        end
    end
end



% sensing config can not visualize own cell
for i = 1:numel(map_candConf)    
    V(i,(i-1)*f+1:(i-1)*f+f) = 0;
end

oV = V; % oV := visibility with obstacles

% V = sparse(V);

%{
% OBSTACLE FREE VISIBILITY:
% ------------------------------------------------------------------------------
% -- removing cells
V(OBS.ind,:) = [];

% -- removing conf.
confOBS = zeros(numel(OBS.ind)*f,1);
for i = 1:numel(OBS.ind)
    % removing sensing configurations.
    %V(:,((OBS.ind(i)-1)*f+1)-(f*(i-1)):((OBS.ind(i)-1)*f+f)-(f*(i-1))) = [];
    confOBS(((i-1)*f)+1:((i-1)*f)+f) = (((OBS.ind(i)-1)*f)+1):(((OBS.ind(i)-1)*f)+f);
end
V(:,confOBS) = [];
%}

% VISIBILITY MATRIX CORRESPOND TO THE HOTSPOT ONLY:
% ------------------------------------------------------------------------------

% ---- cells:

% cells covered by candidate conf.
% covered_ind = ( spotV*ones(size(spotV,2),1) ); % or
covered_ind = sum(oV,2);
covered_num = find(covered_ind);
covered_cell_vector = zeros(1,size(oV,1));
covered_cell_vector(covered_num) = 1;

% cells need to be covered -- according to the coverage map
coverage_ind = find(map_coverage);
coverage_cell_vector = zeros(1,size(oV,1));
% size(coverage_cell_vector)
coverage_cell_vector(coverage_ind) = 1;
% size(coverage_cell_vector)

% unoccupied cells
free_cell_vector = ones(1,size(oV,1));
free_cell_vector(OBS.ind) = 0;

% size(coverage_cell_vector)

% pause
% --- selecting cells that are covered and need to be covered and unoccupied.
% size(covered_cell_vector)
% size(coverage_cell_vector)
% size(free_cell_vector)
cov_cell     = covered_cell_vector.*coverage_cell_vector.*free_cell_vector;
cellKept     = find(cov_cell);
cellRemoved  = find(cov_cell==0);
V            = V(cellKept,:);

% ---- configurations:

% --- selecting conf number corresponding to candidate conf map only
% cells for the candidate conf.
candConfCells = find(map_candConf); 
% conf numbers
% candConfNum = zeros(numel(map_candConf)*f,1);
candConfNum = zeros(1,size(oV,2));
for i = 1:numel(candConfCells)
    candConfNum( (candConfCells(i)-1)*f+1:(candConfCells(i)-1)*f+f ) = 1;
end

% conf that cover hotspot
tmpV = oV;
tmpV = tmpV(hotspot,:);
conf_hotspot_num = find(tmpV);
conf_hotspot = zeros(1,size(oV,2));
conf_hotspot(conf_hotspot_num) = 1;

% conf belong to unoccupied cells
conf_free = ones(1,size(oV,2));
for i = 1:numel(OBS.ind)
    conf_free( (OBS.ind(i)-1)*f+1:(OBS.ind(i)-1)*f+f ) = 0;
end

% --- keeping candidate conf in the list that are from candidate conf map
% AND cover the hotspot AND are placed over unoccupied cells 
conf_vector  = candConfNum.*conf_hotspot.*conf_free;
confKept     = find(conf_vector);
confRemoved  = find(conf_vector==0);
% V            = V(:,confKept);
% V            = full(V);

% ------------- limited number of varibales to x --------------------------
% removing cand conf close to the hotspot to keep a limited number of cand
% conf.

% --- lets remove cand conf very close to the hotspot (say 1m on 0.5m
% resolution map).

innerRadius_cell = para_.innerRadius_m/dataYAML.resolution;

% temporary inner limit map
map_candConf_minLimit = zeros(size(map_env));
map_candConf_minLimit(hotspot) = 1;

se = strel('disk',innerRadius_cell,0);
map_candConf_minLimit = imdilate(map_candConf_minLimit,se);
% figure; imshow(map_candConf_minLimit)

% indices of inner limit map
innerLimit_ind = find(map_candConf_minLimit);

%  conf kept cell numbers
confKeptCells = fix(confKept/f)+(~~mod(confKept,f));

% conf.s that are on the inner limit cells
confToRemove_ind = find(ismember(confKeptCells,innerLimit_ind));

% update confRemoved list
confRemoved = union(confRemoved,confKept(confToRemove_ind));

% remove cand conf
confKept(confToRemove_ind) = [];


% --- now, lets keep limited number of cand conf by removing excessive conf
% from outer ring, iteratively.

i = 0;
while numel(confKept)>para_.NumOfCandidateConfTomography    
    i = i+1;
    
    outRadius_cell = SensorRange.cell-(i-1);
    inRadius_cell  = SensorRange.cell-(i-0);
    
    
    %innerRadius_cell = i;
    %innerRadius_cell = paramPreprocess.innerRadius_m/paramPreprocess.MapRes;

    % temporary inner limit map
    map_candConf_minLimit = zeros(size(map_env));
    map_candConf_minLimit(hotspot) = 1;

    se = strel('disk',inRadius_cell,0);
    map_candConf_minLimit = imdilate(map_candConf_minLimit,se);
    %figure('name','min limit'); imshow(map_candConf_minLimit)
        
    % temporary outer limit map
    map_candConf_maxLimit = zeros(size(map_env));
    map_candConf_maxLimit(hotspot) = 1;

    se = strel('disk',outRadius_cell,0);
    map_candConf_maxLimit = imdilate(map_candConf_maxLimit,se);
    %figure('name','max limit'); imshow(map_candConf_maxLimit)
    
    % differential map
    map_candConf_outerRing = map_candConf_maxLimit-map_candConf_minLimit;
    %figure('name','outer ring'); imshow(map_candConf_outerRing)

    % indices of inner limit map
    outerRing_ind = find(map_candConf_outerRing);
    
    %  conf kept cell numbers
    confKeptCells = fix(confKept/f)+(~~mod(confKept,f));
        
    % conf.s that are on the inner limit cells
    confToRemove_ind = find(ismember(confKeptCells,outerRing_ind));
        
    % update confRemoved list
    confRemoved = union(confRemoved,confKept(confToRemove_ind));
    
    % remove cand conf
    confKept(confToRemove_ind) = [];
    
    %size(confToRemove_ind)   
    %size(confKept)   
    %pause
end
%}

% --- update visibility matrix for the conf kept only.
V = V(:,confKept);
V = full(V);



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

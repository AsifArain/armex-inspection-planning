function [ oConfOrns,...
           ConfOrns,...
           SensorRange,...
           OBSenv,...
           OBSconf,...
           V,...
           oV,...
           wV,...
           Conf_DistOrnGainV,...
           confKept,...
           confRemoved,...
           cellKept,...
           cellRemoved,...
           FoVfaces,...
           map_candConf,...
           map_coverage,...
           refDistPoint,...
           conf_crossAngles,...
           Angles_FromHc2Conf,...
           D,...
           tPreprocess ] = fVDHc2t( FoV,...
                                    OBSenv,...
                                    OBSconf,...
                                    SensorRange,...
                                    map_coverage,...
                                    engagedConfs_ind,...
                                    map_recon,...
                                    hc_ind,...
                                    FocusTSPAng,...
                                    para_,...
                                    dir_,...
                                    visualize )
% fPreProcess_R4 genereates required parameters for solving a global coverage problem.
% Date: 2014-01-19, Rev 4
% 
% n  := total number of cells in the map.
% nf := number of free cells in the map.
% f  := total number of sensing configurations per cell.
% 
% OUTPUTS:
%-----------------------------------------------------------------------------------------
% cnt: [scalar][ ]   4 or 8 connectivity graph.
% o  : [scalar][ ]   Number of outgoing edges for each conf.
% E  : [nf*f*o][ ]   Edge connectivity. Index number is ID of outgoing edge
% and array value is ID of connected incoming edge.%
% vE : [------][ ]   E (as above) with reduced variables using Corner Reduction Method.
% eE : [nf*f*o][ ] 
% evE: [------][ ]   eE (as above) with reduced variables using Corner Reduction Method.
% T  : [nf*f,o][sec] Travelling cost. Rows are conf number and columns are outgoing edges.
% vT : [------][sec] T (as above) with reduced variables using Corner Reduction Method.
% V  : [nf,nf*f][binary] Visibility matrix, indicates visiblity status (1
% for visible, 0 otherwise) of nf cells (rows) from each conf (columns).
% vV : [-------][binary] V (as above) with reduced variables using Corner
% Reduction Method. 
% vkS: [][]          .cell := selected special vertex cell number.
%                    .conf := configurations that cover selected special vertex cell.
% CRL: [][]          Corner Reduction List, list of configuration selected
% to be removed using Corner Reduction technique. 
% SensorRange: [scalar][integer] .cell sensor range in number of cells.
% map: [l,m][binary] Grid map matrix, 0 for occupied cells and 1 for
% unoccupied cells. l and m are any positive integer numbers. 
% FoVfaces [][]: .alf first/start angle for each FoV orientation.
%                .bay second/final angle for each FoV orientation.
%                .ang set of first/start and second/final angles of each FoV orientation.
%                .num total number of FoV faces/orientations for a sensing position.
% OBS: [][]      .ind Indices of occupied cells.
%                .sub Subscripts (row,col) of occupied cells.
% tPreprocess: [scalar][sec] Total computation time.
% 
% INPUTS:
%-----------------------------------------------------------------------------------------
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
%                     .SenDom   := computing sensor domain. 'ComputeNew' to
%                     compute from scratch, and 'UseEarlier' to use earlier
%                     computed sensor domain parameters. 
% map          : [l,m][binary] Grid map matrix, 0 for occupied cells and 1
% for unoccupied cells. l and m are any positive integer numbers. 
%
% NOTES/UPDATES:
% ----------------------------------------------------------------------------------------
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




%***************************************
%                 MAP:
%***************************************
[ map_env,...
  map_conf,...
  origin_env,...
  origin_conf,...
  cellsize_env,...
  cellsize_conf,...
  cellsize_org ] = fMapStuff( para_,dir_ );




%**********************************************************************
% SENSING PARAMETERS AND CONSTANTS:
%**********************************************************************
FoVfaces.num        = para_.numFoVs;   % number of conf per cell.
n                   = numel(map_env);         % number of cells in the map.
f                   = FoVfaces.num;       % for simplification.
c                   = n*f;                % total conf.
% cnt                 = para_.GraphCnt;  % 4 or 8 connected neighbours.

% o := number of outgoing edges for each configuration.
%{
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
%}

%**********************************************************************
% CANDIDATE CONF MAP BASED ON SENSING RANGE:
%**********************************************************************

fprintf(1,'---- candidate conf based on sensing range \n')

[hc_r,hc_c] = ind2sub(size(map_env),hc_ind);

[ map_candConf ] = fConfMapSensingRange( hc_ind,...
                                         map_conf,...
                                         SensorRange,...
                                         OBSconf,...
                                         para_,...
                                         visualize );

%**********************************************************************                                     
% VISIBILITY MATRIX:
%**********************************************************************
fprintf(1,'---- visibility matrix \n')

executedConfsHot_ind = [];
executedConfsHot_orn = [];
[hc_r,hc_c] = ind2sub(size(map_env),hc_ind);
hc_this_sub = [hc_r,hc_c];

[ V,...
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
                                          visualize );

%**********************************************************************
% CANDIDATE CONF MAP BASED ON SENSING RANGE AND VISIBILITY (FINAL):
%**********************************************************************
fprintf(1,'---- candidate confs map \n');

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
plot(hc_r,hc_c,'or','MarkerFaceColor','r');
end


%**********************************************************************
% CROSS ANGLES:
%**********************************************************************
fprintf(1,'---- cross angles \n');

[ conf_crossAngles,...
  Angles_FromHc2Conf ] = fX( hc_ind,...
                             map_env,...
                             confKept );
                         
%**********************************************************************
% REFERENCE POINT FOR DISTANCE:
%**********************************************************************
fprintf(1,'---- reference points for traveling distances \n');

[ BestFocusPointGain,...
  BestFocusPointAng,...
  refDistPoint ] = fReferencePointDistance( map_candConf,...
                                            confKept,...
                                            FocusTSPAng,...
                                            hc_r,...
                                            hc_c,...
                                            SensorRange,...
                                            conf_crossAngles,...
                                            OBSconf,...
                                            para_,...
                                            visualize );


%**********************************************************************
% DISTANCE MATRIX:
%**********************************************************************

fprintf(1,'---- distance matrix \n')


[ D ] = fD( refDistPoint,confKept,map_env );


%**********************************************************************
% Total computation time.
%**********************************************************************
tPreprocess = 1e-4*(round(toc(tPreprocessi)*1e4));

fprintf(1,'---- Computation time: %0.4f sec \n',tPreprocess)




end



function [ map_candConf ] = fConfMapSensingRange( hc_ind,...
                                                  map_env,...
                                                  SensorRange,...
                                                  OBSconf,...
                                                  para_,...
                                                  visualize )

[hc_r,hc_c] = ind2sub(size(map_env),hc_ind);

% initialize map for candidate conf. -- ring only
map_candConf = zeros(size(map_env));    
map_candConf(hc_ind) = 1;

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
plot(hc_r,hc_c,'xr','MarkerFaceColor','r');
end

% -- remove obs cells
map_candConf(OBSconf.ind(ismember(OBSconf.ind,find(map_candConf)))) = 0;

% plot
if visualize == 1
figure('name','cand conf - max limit and obs are removed');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
map_candConf_ind = find(map_candConf(:)==1);
[map_candConf_r,map_candConf_c] = ind2sub(size(map_env),map_candConf_ind);
plot(map_candConf_r,map_candConf_c,'xb','MarkerFaceColor','b');        
plot(hc_r,hc_c,'xr','MarkerFaceColor','r');
end


%               TMP

%-- Limit the number of candidate confs before computing visibility (1SE only)
%==============================================================================

% refConfCell = confs_executed_hc(end,1:2);
candConf_i = find(map_candConf);

maxConfs = para_.maxCandidateConf_Hotspot_2t+para_.deltaCandConf_Hotspot_2t;
% if numel(candConf_i) > (para_.NumOfCandidateConfTomography1SE+20)
if numel(candConf_i) > maxConfs

    %[candConf_r,candConf_c] = ind2sub(size(map_env),candConf_i);
    %confdist = pdist2([candConf_r,candConf_c],refConfCell,'euclidean');    
    % -- sorted indices
    %[~,sorted_ind] = sort(confdist,'ascend');
    %sorted_ind_max = sorted_ind(1:min([para_.NumOfCandidateConfTomography1SE+20,...
    %    numel(sorted_ind)]));
    
    rnd_series = randperm(numel(candConf_i),numel(candConf_i));
    [~,sorted_ind] = sort(rnd_series,'ascend');
    
    %sorted_ind_max = sorted_ind(1:min([para_.NumOfCandidateConfTomography1SE+20,...
    %    numel(sorted_ind)]));
    
    sorted_ind_max = sorted_ind(1:min([maxConfs,numel(sorted_ind)]));    
    
    candConfKept_i = candConf_i(sorted_ind_max);
    
    map_candConf = zeros(size(map_env));
    map_candConf(candConfKept_i) = 1;
    
end




%--- removing conf from edges 
%-----------------------------------
% warning('The edges will not be removed since we have already extracted a conf map.')
%{
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
%}


end


function [ D ] = fD( refDistPoint,confKept,map_env )

%
%
%

refDistPoint_r = refDistPoint(1);
refDistPoint_c = refDistPoint(2);


% ---- D (distance) matrix is manhattan distance between cand conf position
% and reference position
confKeptCells_ind = confKept; %fix(confKept/f)+(~~mod(confKept,f));
[confKeptCells_r,confKeptCells_c] = ind2sub(size(map_env),confKeptCells_ind);

D = zeros(numel(confKeptCells_ind),1);
for i = 1:numel(confKeptCells_ind)
    manhattan_distance = ...
        abs(confKeptCells_r(i)-refDistPoint_r)+abs(confKeptCells_c(i)-refDistPoint_c);
    D(i) = manhattan_distance;
end



end



function [ BestFocusPointGain,...
           BestFocusPointAng,...
           refDistPoint ] = fReferencePointDistance( map_candConf,...
                                                     confKept,...
                                                     FocusTSPAng,...
                                                     hc_r,...
                                                     hc_c,...
                                                     SensorRange,...
                                                     conf_crossAngles,...
                                                     OBSconf,...
                                                     para_,...
                                                     visualize )



[map_row,map_col] = size(map_candConf);

% cell numbers of all valid candidate conf
% AllCandConfCell = fix(confKept/f)+(~~mod(confKept,f)); % cell num for each conf num.
AllCandConfCell = confKept; % cell num for each conf num.

% --- light angles of possible areas for distance point
% LightAngles = [0,90,180,270]
% deltaAng = 10;
% paramPreprocess.AngularStepForDistanceFocusPoint = 10;
LightAngles = 0:para_.AngularStepForDistanceFocusPoint:360-para_.AngularStepForDistanceFocusPoint;
LightAngles = LightAngles+FocusTSPAng;
LightAngles = LightAngles+(360*(LightAngles<0))-(360*(LightAngles>360));

% --- cand conf domain for each focus point
sensor = [hc_r,hc_c];
% DomainRadius.cell = radius_candCircleOut;
DomainRadius.cell = SensorRange.cell/2;

% --- initialize variables
AllFcousPoints_NumOfCandConfs = zeros(numel(LightAngles),1);
AllFcousPoints_NumOfFreeCells = zeros(numel(LightAngles),1);
AllFcousPoints_CircmCandConfs = zeros(numel(LightAngles),1);
for i = 1:numel(LightAngles)
    
    % -- domain of focus point
    %zawia = [LightAngles(i)-90,LightAngles(i),LightAngles(i)+90];
    %zawia+(360*(zawia<0))-(360*(zawia>360));
        
    deltaZawia = 1;
    zawia = LightAngles(i)-90:deltaZawia:LightAngles(i)+90;
    zawia = wrapTo360(zawia);
    
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
                          numel(find(ismember(OBSconf.ind,FocusPointDomainValid_ind)))-...
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

% AllFcousPoints_NumOfCandConfs

% normalized - number of cand conf, number of free cells and circumference
AllFcousPoints_NumOfCandConfsNorm = AllFcousPoints_NumOfCandConfs/max(AllFcousPoints_NumOfCandConfs);
AllFcousPoints_NumOfFreeCellsNorm = AllFcousPoints_NumOfFreeCells/max(AllFcousPoints_NumOfFreeCells);
AllFcousPoints_CircmCandConfsNorm = AllFcousPoints_CircmCandConfs/max(AllFcousPoints_CircmCandConfs);

% AllFcousPoints_CircmCandConfs

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
    endPoint_r = hc_r + (SensorRange.cell/4*cosd(LightAngles(i)));
    endPoint_c = hc_c + (SensorRange.cell/4*sind(LightAngles(i)));
    
    if visualize == 1
    plot(endPoint_r,endPoint_c,'+g');
    plot([hc_r,endPoint_r],[hc_c,endPoint_c],'-g'); 
    end
    
end

% refDistPoint_r = hotspot_r + (DomainRadius.cell*cosd(BestFocusPointAng))
% refDistPoint_c = hotspot_c + (DomainRadius.cell*sind(BestFocusPointAng))
% 
% refDistPoint = [refDistPoint_r,refDistPoint_c]

% --- lets go back to TSP angles
% refDistPoint_r = hotspot_r + (DomainRadius.cell*cosd(FocusTSPAng))
% refDistPoint_c = hotspot_c + (DomainRadius.cell*sind(FocusTSPAng))
refDistPoint_r = hc_r + (SensorRange.cell/4*cosd(FocusTSPAng));
refDistPoint_c = hc_c + (SensorRange.cell/4*sind(FocusTSPAng));

refDistPoint = [refDistPoint_r,refDistPoint_c];
% DomainRadius.cell

if visualize == 1
% plot(refDistPoint_c,refDistPoint_r,'*r');
plot(refDistPoint_r,refDistPoint_c,'*r');
end






end


function [ conf_crossAngles,...
           Angles_FromHc2Conf ] = fX( hc_ind,...
                                      map_env,...
                                      confKept )


                                 
%***********************************************************************
%
% CONF ANGLES TOWARDS HC
%
%***********************************************************************
% fprintf('angles towards Hc...\n')

confKeptCell_ind = confKept;
[confKeptCell_r,confKeptCell_c] = ind2sub(size(map_env),confKeptCell_ind);
sensor_placements = [confKeptCell_r',confKeptCell_c'];

[hc_r,hc_c] = ind2sub(size(map_env),hc_ind);

Angles_FromHc2Conf = angle2Points([hc_r,hc_c],sensor_placements);


%***********************************************************************
%
% CROSS ANGLES
%
%***********************************************************************
% disp('----- cross angles')

conf_crossAngles = zeros(numel(confKept));
for i = 1:numel(confKept)
    for j = 1:numel(confKept)
        conf_crossAngles(i,j) = ...
            rad2deg(angleAbsDiff(Angles_FromHc2Conf(i),Angles_FromHc2Conf(j)));
    end
end

                                  

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



% function [ V,...
%            oV,...
%            wV,...
%            Conf_DistOrnGainV,...
%            confKept,...
%            confRemoved,...
%            cellKept,...
%            cellRemoved,...
%            oConfOrns,...
%            ConfOrns,...
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
%                             OBSenv,...
%                             OBSconf,...
%                             hc_ind,...
%                             gamma,...
%                             cellsize_env,...
%                             para_,...
%                             visualize )
% %fV generates Visibility matrix V.
% 
% % --- ANGLES:
% FoVfaces.lgt = [1/FoVfaces.num:1/FoVfaces.num:1-(1/FoVfaces.num) 0]'*360;
% FoVfaces.alf = FoVfaces.lgt-(FoV/2);
% FoVfaces.bay = FoVfaces.lgt+(FoV/2);
% % compensate negative angles
% %FoVfaces.alf(FoVfaces.alf<0) = 360+FoVfaces.alf(FoVfaces.alf<0);
% FoVfaces.alf = wrapTo360(FoVfaces.alf);
% FoVfaces.bay = wrapTo360(FoVfaces.bay);
% FoVfaces.ang = [FoVfaces.alf FoVfaces.bay];
% FoVfaces.ang = [FoVfaces.alf FoVfaces.bay];
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %               V
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % disp('.. Visibility Matrix')
% % V = sparse(n,c); % initializing
% V = sparse(n,n); % initializing
% % V = zeros(n,c); % initializing
% [hc_r,hc_c] = ind2sub(size(map_env),hc_ind);
% 
% oConfOrns = zeros(1,n);
% % For each cell of the map.
% for i = 1:n
%     if map_candConf(i) % if its not occupied
%         [rr,cc] = ind2sub(size(map_candConf),i);
%         sensor = [rr,cc];
%         % --
%         orientation = round(rad2deg(angle2Points(sensor,[hc_r,hc_c])));
%         orientation = orientation + (orientation==0)*360;
% 
%         for j = orientation %1%:f
%             
%             sD = []; oD = []; VObs = [];
%             %this_dir = sprintf(['sD_oD_VObs/sensing_range_%02dcell/fov_%03ddeg/',...
%             %                    'orientation_%03ddeg/'],SensorRange.cell,FoV,j);
%             this_dir = sprintf('visibility/templates/sensing_range_%02dcell/fov_%03ddeg/',...
%                 SensorRange.cell,FoV);
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
%         end
%         oConfOrns(i) = orientation;
%     end
% end
% 
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
% % V = sparse(V);
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
% % ////////////////////////////
% % -- cellKept weights
% % ////////////////////////////
% % cellKeptWeights = map_recon(cellKept);
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
% 
% % -- plot 
% if visualize == 1
% figure('name','conf corresponding to candidate conf map only');
% imshow(map_env','InitialMagnification',800); hold on;
% set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
% candConfNum_ind = find(candConfNum(:)==1);
% [candConfNum_r,candConfNum_c] = ind2sub(size(map_env),candConfNum_ind);
% plot(candConfNum_r,candConfNum_c,'xb','MarkerFaceColor','b');        
% plot(hc_r,hc_c,'or','MarkerFaceColor','r');
% end
% 
% %----------------------------------------------------
% % -- Check#2: confs. that cover all the hotspots
% %----------------------------------------------------
% conf_hotspot = zeros(1,size(oV,2));
% 
% % conf ind that cover this hotspot
% conf_ind = find(oV(hc_ind,:));    
% conf_hotspot(conf_ind) = 1;
%     
% 
% % -- plot 
% if visualize == 1
% figure('name','conf cover all the hotspots');
% imshow(map_env','InitialMagnification',800); hold on;
% set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
% conf_hotspot_ind = find(conf_hotspot(:)==1);
% [conf_hotspot_r,conf_hotspot_c] = ind2sub(size(map_env),conf_hotspot_ind);
% plot(conf_hotspot_r,conf_hotspot_c,'xb','MarkerFaceColor','b');        
% plot(hc_r,hc_c,'or','MarkerFaceColor','r');
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
% 
% 
% % -- plot 
% if visualize == 1
% figure('name','conf belong to unoccupied cells');
% imshow(map_env','InitialMagnification',800); hold on;
% set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
% conf_free_ind = find(conf_free(:)==1);
% [conf_free_r,conf_free_c] = ind2sub(size(map_env),conf_free_ind);
% plot(conf_free_r,conf_free_c,'xb','MarkerFaceColor','b');
% plot(hc_r,hc_c,'or','MarkerFaceColor','r');
% end
% 
% 
% % ----------------------------------------------------
% % -- Check#4: conf not close enough to the hotspot 
% % -----------------------------------------------------
% 
% % - cose enough radius
% %innerRadius_cell = para_.innerRadius_m/dataYAML.resolution;
% innerRadius_cell = para_.innerRadius_m/cellsize_env;
% % temporary inner limit map
% map_candConf_minLimit = zeros(size(map_env));
% map_candConf_minLimit(hc_ind) = 1;
% 
% se = strel('disk',innerRadius_cell,0);
% map_candConf_minLimit = imdilate(map_candConf_minLimit,se);
% 
% % indices of inner limit map
% innerLimit_ind = find(map_candConf_minLimit);
% 
% % -- plot 
% if visualize == 1
% figure('name','conf inner limit');
% imshow(map_env','InitialMagnification',800); hold on;
% set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
% [innerLimit_r,innerLimit_c] = ind2sub(size(map_env),innerLimit_ind);
% plot(innerLimit_r,innerLimit_c,'xb','MarkerFaceColor','b');
% plot(hc_r,hc_c,'sr');
% end
% 
% confNotClose = ones(1,size(oV,2));
% confNotClose(innerLimit_ind) = 0;
% 
% 
% % -- plot 
% if visualize == 1
% figure('name','conf Not Close');
% imshow(map_env','InitialMagnification',800); hold on;
% set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
% [confNotClose_r,confNotClose_c] = ind2sub(size(map_env),find(confNotClose));
% plot(confNotClose_r,confNotClose_c,'xb','MarkerFaceColor','b');
% plot(hc_r,hc_c,'sr');
% end
% 
% 
% 
% %-----------------------------------------------------------------------------
% % -- Check#5: confs (position) is not already being used for the other hcs
% %-----------------------------------------------------------------------------
% confNotEngaged = ones(1,size(oV,2));
% confNotEngaged(engagedConfs_ind) = 0;
% 
% 
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
% 
% conf_vector  = candConfNum.*conf_hotspot.*conf_free.*confNotClose.*confNotEngaged;
% confKept     = find(conf_vector);
% confRemoved  = find(conf_vector==0);
% % V            = V(:,confKept);
% % V            = full(V);
% 
% % -- plot 
% if visualize == 1
% figure('name','conf meeting all four conditions');
% imshow(map_env','InitialMagnification',800); hold on;
% set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
% [confKept_r,confKept_c] = ind2sub(size(map_env),confKept);
% plot(confKept_r,confKept_c,'xb','MarkerFaceColor','b');
% plot(hc_r,hc_c,'or','MarkerFaceColor','r');
% end
% 
% 
% % **********************************************************
% % -- Limit Number of Cand Conf based on Proportional Coverage
% % **********************************************************
% 
% %/////////////////////////////////////// 
% % -- Concentration weights
% %///////////////////////////////////////
% 
% cellKeptWeights = map_recon(cellKept);
% %cellKeptWeights = cellKeptWeights/(max(cellKeptWeights));
% 
% V_this = V(:,confKept);
% V_this = V_this .* repmat(cellKeptWeights',1,numel(confKept));
% wV_concentration_this = sum(V_this,1);
% wV_concentration_this = wV_concentration_this/max(wV_concentration_this(:));
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
% [hc_r,hc_c] = ind2sub(size(map_env),hc_ind);
% [conf_r,conf_c] = ind2sub(size(map_candConf),confKept);
% 
% desiredOrns = round(rad2deg(angle2Points([conf_r',conf_c'],[hc_r,hc_c])));
% desiredOrns = desiredOrns + (desiredOrns==0)*360;
% 
% % -- delta orientations
% ConfDeltaOrns = rad2deg(angleAbsDiff(deg2rad(desiredOrns(:)),deg2rad(ConfOrns(:))));
% 
% % -- angular gain
% gain_sigma = FoV/6;%15; % standard deviation
% gain_meu   = 0; % expected value or center
% V_angular_gain_this = gaussmf(ConfDeltaOrns,[gain_sigma,gain_meu]);
% 
% % ------------------------------
% % == displacement
% % ------------------------------
% 
% % -- dist from conf to hotspot
% conf2Hot_dist = pdist2([conf_r',conf_c'],[hc_r,hc_c],'euclidean');
% 
% % -- distance gain
% gain_sigma       = SensorRange.cell/5; % standard deviation
% gain_meu         = SensorRange.cell/2; % expected value or center
% V_dist_gain_this = gaussmf(conf2Hot_dist,[gain_sigma,gain_meu]);
% 
% 
% % ------------------------------
% % == orientation + displacement
% % ------------------------------
% wV_diplacement_this = (V_angular_gain_this+V_dist_gain_this)/2;
% 
% wV_this_a = ((gamma*(wV_concentration_this)))';
% wV_this_b = ((1-gamma)* wV_diplacement_this);
% 
% wV_this = wV_this_a+wV_this_b;
% %wV_this = (gamma*(wV_concentration_this))' + ((1-gamma)* wV_diplacement_this)
% 
% 
% % -- sorted indices
% [~,sorted_ind] = sort(wV_this,'descend');
% sorted_ind_max = sorted_ind(1:min([para_.maxCandidateConfLocal2SE,numel(sorted_ind)]));
% % sorted_ind_max = sorted_ind(1:min([10,numel(sorted_ind)]));
% confToRemove_ind = find(ismember(sorted_ind,sorted_ind_max)==0);
% confRemoved = union(confRemoved,sorted_ind(confToRemove_ind));
% confKept = confKept(sorted_ind_max);
% confKept = sort(confKept);
% 
% 
% % -- plot 
% if visualize == 1
% figure('name','max num of conf');
% imshow(map_env','InitialMagnification',800); hold on;
% set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
% [confKept_r,confKept_c] = ind2sub(size(map_env),confKept);
% plot(confKept_r,confKept_c,'xb','MarkerFaceColor','b');
% plot(hc_r,hc_c,'sr');
% end
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
% % wV = V.*repmat(cellKeptWeights',1,numel(confKept));
% wV = V.*repmat(map_recon(cellKept)',1,numel(confKept));
% max_conc_conf = max(sum(wV,1));
% wV = wV./max_conc_conf;
% 
% % sum(map_recon(cellKept))
% % wV = wV/sum(map_recon(cellKept));
% % pause
% 
% 
% % wV = wV/3;
% % sum(wV,1)
% % pause
% % numel(cellKept)
% % size(wV)
% % size(sum(wV,1))
% % max(sum(wV,1))
% % pause
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
% [hc_r,hc_c] = ind2sub(size(map_env),hc_ind);
% [conf_r,conf_c] = ind2sub(size(map_candConf),confKept);
% 
% desiredOrns = round(rad2deg(angle2Points([conf_r',conf_c'],[hc_r,hc_c])));
% desiredOrns = desiredOrns + (desiredOrns==0)*360;
% 
% % -- delta orientations
% ConfDeltaOrns = rad2deg(angleAbsDiff(deg2rad(desiredOrns(:)),deg2rad(ConfOrns(:))));
% 
% % -- angular gain
% gain_sigma = FoV/6;%15; % standard deviation
% gain_meu   = 0; % expected value or center
% V_angular_gain = gaussmf(ConfDeltaOrns,[gain_sigma,gain_meu]);
% 
% %//////////////////////////////////////////
% % HOTSPOT PLACEMENT GAIN - DISTANCE
% %//////////////////////////////////////////
% 
% % -- dist from conf to hotspot
% conf2Hot_dist = pdist2([conf_r',conf_c'],[hc_r,hc_c],'euclidean');
% 
% % -- distance gain
% gain_sigma = SensorRange.cell/5; % standard deviation
% gain_meu   = SensorRange.cell/2; % expected value or center
% V_dist_gain = gaussmf(conf2Hot_dist,[gain_sigma,gain_meu]);
% 
% 
% %//////////////////////////////////////////////
% % HOTSPOT PLACEMENT GAIN - ANGULAR + DISTANCE
% %//////////////////////////////////////////////
% Conf_DistOrnGainV = (V_angular_gain+V_dist_gain)/2;
% 
% 
% % ------ publish gussian gain function for test only
% %{ 
% angles = 0:1:15; % in degrees
% gain_sigma = 3; % standard deviation
% gain_meu   = 7.5; % expected value or center
% guss_gain = gaussmf(angles,[gain_sigma,gain_meu]);
% figure; plot(angles,guss_gain ,'-b')
% %}        
% 
% % sum(wV,1)
% % 
% % size(Conf_DistOrnGainV)
% % % 
% % % one = (0.9*wV)
% % % 
% % % two = 0.1*repmat(Conf_DistOrnGainV',numel(cellKept),1);
% % 
% % sum(0.9*wV,1)
% % 0.1*(Conf_DistOrnGainV')
% % 
% % wV1 = sum(0.9*wV,1) + 0.1*(Conf_DistOrnGainV')
% % % wV = one+two
% % % 
% % % sum(one,1)
% % % 
% % % sum(two,1)
% % 
% % 
% % sum(wV,1)
% % 
% % size(wV)
% % 
% % Conf_DistOrnGainV(40)
% % 
% % pause
% 
% 
% 
% end
% 



function [ visibleRange ] = fSenDomain( sensor,SensorRange,zawia )
%fSenDomain returns visible cells (domain/mask) for a given sensing
%configuration. 
%

% INTERPOLATING INTERMEDIATE ANGLES:
% -------------------------------------
%{
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
%}

DELTA_ANGLE = 1;
if     zawia(2) > zawia(1)
    angl_1h = zawia(1):DELTA_ANGLE:zawia(2);
elseif zawia(2) < zawia(1) % if final angle is greater than 360 degrees.
    angl_1h = [zawia(1):DELTA_ANGLE:360 0:DELTA_ANGLE:zawia(2)];
end

% Second half:
if     zawia(3) > zawia(2)
    angl_2h = zawia(2):DELTA_ANGLE:zawia(3);    
elseif zawia(3) < zawia(2) % if final angle is greater than 360 degrees.
    angl_2h = [zawia(2):DELTA_ANGLE:360 0:DELTA_ANGLE:zawia(3)];
end

angl = unique([angl_1h angl_2h]);
angl = wrapTo360(angl);

% angl = zawia
visibleRange = [];        %
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
%{
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
%}

DELTA_ANGLE = 1;
if     zawia(2) > zawia(1)
    angl_1h = zawia(1):DELTA_ANGLE:zawia(2);
elseif zawia(2) < zawia(1) % if final angle is greater than 360 degrees.
    angl_1h = [zawia(1):DELTA_ANGLE:360 0:DELTA_ANGLE:zawia(2)];
end

% Second half:
if     zawia(3) > zawia(2)
    angl_2h = zawia(2):DELTA_ANGLE:zawia(3);    
elseif zawia(3) < zawia(2) % if final angle is greater than 360 degrees.
    angl_2h = [zawia(2):DELTA_ANGLE:360 0:DELTA_ANGLE:zawia(3)];
end

angl = unique([angl_1h angl_2h]);
angl = wrapTo360(angl);

% angl = zawia;

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
cond_x1 = ang>=000&&ang<=045;
cond_x2 = ang>=135&&ang<=225;
cond_x3 = ang>=315&&ang<=360;

cond_y1 = ang>=045&&ang<=135;
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

function [ visDom ] = fVisibileDomain( map_env,sensorDom,obstacDom,VOBS )
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
[map_row,map_col] = size(map_env);

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
obsDomX_ind = sub2ind(size(map_env),obsDomX(:,1),obsDomX(:,2));

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
    visDom = sort(sub2ind(size(map_env),visDom(:,1),visDom(:,2)));
end

end

% --------------------------- End of the document --------------------------------------

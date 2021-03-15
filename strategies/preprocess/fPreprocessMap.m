
function [ map_env,...
           map_recon,...
           map_coverage,...
           OBSenv,...
           OBSconf,...
           o,...
           E,...
           T,...
           hotspotSEQ,...
           InOutTSPAnglesAll,...
           InOutTSPAnglesMean,...
           FocusTSPAngles,...
           FoVfaces,...
           SensorRange,...
           tPreprocess ] = fPreprocessMap( FoV,...
                                           SensorRange,...
                                           cell_coords_x,...
                                           cell_coords_y,...
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
% ----------------------------------------------------------------------------------------
% cnt: [scalar][ ]   4 or 8 connectivity graph.
% o  : [scalar][ ]   Number of outgoing edges for each conf.
% E  : [nf*f*o][ ]   Edge connectivity. Index number is ID of outgoing edge
% and array value is ID of connected incoming edge.%
% vE : [------][ ]   E (as above) with reduced variables using Corner Reduction Method.
% eE : [nf*f*o][ ] 
% evE: [------][ ]   eE (as above) with reduced variables using Corner Reduction Method.
% T  : [nf*f,o][sec] Travelling cost. Rows are conf number and columns are
% outgoing edges. 
% vT : [------][sec] T (as above) with reduced variables using Corner Reduction Method.
% V  : [nf,nf*f][binary] Visibility matrix, indicates visiblity status (1
% for visible, 0 otherwise) of nf cells (rows) from 
%                        each conf (columns).
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
% ---------------------------------------------------------------------------------------
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
%                     compute from scratch, and 'UseEarlier' to use 
%                                  earlier computed sensor domain parameters.
% map          : [l,m][binary] Grid map matrix, 0 for occupied cells and 1
% for unoccupied cells. l and m are any positive 
%                              integer numbers. 
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

% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             RESOLUTION PARAMETERS:
% 
% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤

% % Number of pixels per meter.
% pix_m = 1/dataYAML.resolution;
% 
% % Number of pixels in a cell for desired cell size.
% pix_cell = para_.EnvironmentCellSizeM*pix_m;
% 
% % Cell length in meter.
% cell_length_m = pix_cell/pix_m;
% 
% % Sensor range in number of cells.
% SensorRange.cell = (SensorRange.m/cell_length_m); % sensor range in cells


%**********************************************************************
% 
%   ENVIRONMENT MAP:
% 
%**********************************************************************
% fprintf(1,'::::: environment map \n\n')

% [ map_env,...
%   mapCut,...
%   robotStartCellInMap ] = fMapEnv( para_,dir_,visualize );


%**********************************************************************
% MAP:
%**********************************************************************

fprintf(1,'---> reterive environment map \n')


file__map_env   = sprintf('%s_map_coverage.dat',para_.environment);
file__map_conf  = sprintf('%s_map_conf.dat',para_.environment);

file__org_env   = sprintf('%s_origin_coverage.dat',para_.environment);
file__org_conf  = sprintf('%s_origin_conf.dat',para_.environment);

file__cell_env  = sprintf('%s_cellsize_coverage.dat',para_.environment);
file__cell_conf = sprintf('%s_cellsize_conf.dat',para_.environment);
file__cell_org  = sprintf('%s_cellsize_original.dat',para_.environment);

file__cut_org   = sprintf('%s_mapcut_original.dat',para_.environment);

map_env         = load([dir_.env,file__map_env]);
map_conf        = load([dir_.env,file__map_conf]);
origin_env      = load([dir_.env,file__org_env]);
origin_conf     = load([dir_.env,file__org_conf]);
cellsize_env    = load([dir_.env,file__cell_env]);
cellsize_conf   = load([dir_.env,file__cell_conf]);
cellsize_org    = load([dir_.env,file__cell_org]);
mapCut          = load([dir_.env,file__cut_org]);




% Sensor range in number of cells.
SensorRange.cell = (SensorRange.m/cellsize_env); % sensor range in cells

%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             OBSTACLE LIST:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% Finding subscripts of the obstacles.
OBSenv.ind            = find(~map_env);     % linear indices of obstacles.
[OBSenv_r,OBSenv_c]   = find(~map_env);     % rows and columns of obstacles.
OBSenv.sub            = [OBSenv_r OBSenv_c];  % subscripts of obstacles.

OBSconf.ind           = find(~map_conf);     % linear indices of obstacles.
[OBSconf_r,OBSconf_c] = find(~map_conf);     % rows and columns of obstacles.
OBSconf.sub           = [OBSconf_r OBSconf_c];  % subscripts of obstacles.

% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             GAS DISTRIBUTION MAP:
% 
% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
fprintf(1,'---> coarse reconstruction map \n')

[ map_recon ] = fMapRecon( map_env,...
                           OBSenv,...
                           origin_env,...
                           cell_coords_x,...
                           cell_coords_y,...
                           para_,...
                           dir_,...
                           visualize );


%****************************************************************************
% 
%                        GAS COVERAGE MAP
% 
%****************************************************************************
fprintf(1,'---> coverage map (positive concentration) \n')

map_coverage = zeros(size(map_recon));
map_coverage(map_recon(:)>0) = 1;
map_coverage(OBSenv.ind) = 0;


if visualize==1
% -- plot 
figure('name','gas coverage map');
imshow(map_env','InitialMagnification',750); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])

[cov_r,cov_c] = ind2sub(size(map_env),find(map_coverage));
plot(cov_r+0.0,cov_c+0.0,'sb','MarkerFaceColor','b');
       
end

%****************************************************************************
% 
%             GAS DISTRIBUTION MAP - HOTSPOTS (HIGH CONCENTRATION):
% 
%****************************************************************************
fprintf(1,'---> hotspot centers (hcs) \n')

para_.HotspotsType

if strcmp(para_.HotspotsType,'manual')
    
    if strcmp(para_.ExperimentTitle,'corridor03')
        hotspots = [86,34;...
                    60,30;...
                    45,15;...
                    25,29;...
                    07,26];
    end
    if strcmp(para_.ExperimentTitle,'corridor05')
        hotspots = [17,29;...
                    09,17;...
                    27,26;...
                    44,33;...
                    47,28;...
                    45,10;...
                    64,31;...                    
                    82,33]; %pause
    end
    if strcmp(para_.ExperimentTitle,'corridor06')
        hotspots = [16,29;...
                    44,32;...
                    44,08;...
                    63,31]; %pause
    end
    if strcmp(para_.ExperimentTitle,'corridor08')
        hotspots = [08,17;...
                    21,29;...
                    49,28;...
                    46,05;...
                    73,30;...
                    81,35]; %pause
    end
    if strcmp(para_.ExperimentTitle,'corridor10')
        hotspots = [12,14;...
                    15,30;...
                    28,31;...
                    42,20;...
                    68,30;...
                    85,34]; %pause
    end
    if strcmp(para_.ExperimentTitle,'corridor16')
        hotspots = [08,12;...
                    13,31;...
                    29,28;...
                    42,07;...
                    55,29;...
                    83,36]; %pause
    end
    if strcmp(para_.ExperimentTitle,'corridor17')
        hotspots = [08,16;...
                    03,23;...
                    23,26;...
                    42,22;...
                    54,29;...
                    82,34]; %pause
    end
    if strcmp(para_.ExperimentTitle,'corridor18')
        hotspots = [11,13;...
                    29,27;...
                    40,06;...
                    53,32;...
                    84,34]; %pause
    end
    if strcmp(para_.ExperimentTitle,'corridor19')
        hotspots = [09,02;...
                    05,23;...
                    35,32;...
                    51,05;...
                    67,30;...
                    93,35]; %pause
    end
    if strcmp(para_.ExperimentTitle,'corridor20')
        hotspots = [10,12;...
                    24,31;...
                    44,09;...
                    61,34;...
                    97,34]; %pause
    end
    
    if strcmp(para_.ExperimentTitle,'johnsson-metal04')
        hotspots = [17,24;...
                    09,11;...
                    34,14;...
                    34,25]; %pause
    end
    if strcmp(para_.ExperimentTitle,'johnsson-metal05')
        hotspots = [03,21;...
                    03,12;...
                    34,24;...
                    40,15]; %pause
    end
    
    
    if strcmp(para_.ExperimentTitle,'sample-environment02')
        hotspots = [09,08;...
                    16,31;...
                    39,40;...
                    63,31]; %pause
    end
    
    if strcmp(para_.ExperimentTitle,'sample-reconstruction-03')
        hotspots = [25,25]; %pause
    end
    
    
elseif strcmp(para_.HotspotsType,'computed')
    
    if strcmp(para_.MissionStrategy,'2t-armex') ||...
            strcmp(para_.MissionStrategy,'2t-armex-icra16-newHc')
        
        
        %--- to print old procedure -----------------------
        %[ hotspots ] = fHotspotsICRA16( map_env,...
        %                                 map_recon,...
        %                                 SensorRange,...
        %                                 OBSenv,...
        %                                 para_,...
        %                                 cellsize_env,...
        %                                 visualize );
        %--------------------------------------------------
    
        [ hotspots ] = fEstimatedHcs2t( map_env,...
                                        map_recon,...
                                        SensorRange,...
                                        OBSenv,...
                                        cellsize_env,...
                                        para_,...
                                        dir_,...
                                        visualize );
        %{                                          
        [ hotspots ] = fHotspots( map_env,...
                                  map_recon,...
                                  SensorRange,...
                                  OBSenv,...
                                  para_,...
                                  dataYAML,...
                                  dir_,...
                                  visualize );
        %}
                                  
                              
    elseif strcmp(para_.MissionStrategy,'2t-armex-icra16')
        
        visualize = 1;
        [ hotspots ] = fHotspotsICRA16( map_env,...
                                         map_recon,...
                                         SensorRange,...
                                         OBS,...
                                         para_,...
                                         cellsize_env,...
                                         visualize );
    end
    
elseif strcmp(para_.HotspotsType,'manual-GUI')
    
    
    hotspots = para_.HotspotsManualGUI
        
end

% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                       SENSING PARAMETERS AND CONSTANTS:
% 
% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
FoVfaces.num  = para_.numFoVs;   % number of conf per cell.
n             = numel(map_env);         % number of cells in the map.
f             = FoVfaces.num;       % for simplification.
c             = n*f;                % total conf.
cnt           = para_.GraphCnt;  % 4 or 8 connected neighbours.

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


% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                        EDGES AND COST RELATED MATRICS:
% 
% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
fprintf(1,'---> ETAC \n')

[ T,E,eE,A,oA,Ct,P ] = fETAC( map_env,...
                              cnt,...
                              n,...
                              f,...
                              c,...
                              o,...
                              e,...
                              oO,...
                              eO,...
                              cellsize_env,...
                              OBSenv,...
                              para_ );


% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                   HOTSPOTS WITH IN/OUT ANGLES USING TSP
% 
% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
fprintf(1,'---> tsp in/out directions \n')

[ hotspotSEQ,...
  InOutTSPAnglesAll,...
  InOutTSPAnglesMean,...
  FocusTSPAngles ] = fHotspotInOutAngles( hotspots,...
                                          map_env,...
                                          o,...
                                          c,...
                                          f,...
                                          E,...
                                          T,...
                                          OBSenv,...
                                          FoV,...
                                          FoVfaces,...
                                          origin_env,...
                                          SensorRange,...
                                          para_,...
                                          visualize );
%

%% Total computation time.
% ---------------------------------------------------------------------------------------
tPreprocess =  1e-4*(round(toc(tPreprocessi)*1e4));

disp(['---> computation time: ',num2str(tPreprocess),' sec']);

end





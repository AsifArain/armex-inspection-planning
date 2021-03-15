%
% clear all; close all; %clc;
disp('||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||')
disp('||                         VISIBILITY TEST                            ||')
disp('||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||')
disp(' ')


MapSizes       = 50%:1:100;  % map sizes. start frm 95,9:10
MapNums        = 1%:1:10;   % map nums.

SensorRanges_m = 15;       % sensor range in meters.
FoVs           = 90;       % field of views.



%% PRE-PROCESSING PARAMETERS:
% --------------------------------------------------------------------------------------------------------------------------------
ParamPreprocess.CRL      = 0;              % list of corner reduction variables (yes/no).
ParamPreprocess.vkS      = 0;              % selection of special vertex (yes/no).
ParamPreprocess.Ct       = 0;              % travelling cost matrix (yes/no).
ParamPreprocess.GraphCnt = 4;              % graph connectivity (4/8).
% ParamPreprocess.Res      = Resolutions;    % required cell size in meter.
ParamPreprocess.TravCost = 10;             % travelling cost (sec) for one meter travel.
ParamPreprocess.PixSize  = 0.5;           % original (map) pixel size (meters).
ParamPreprocess.numFoVs  = 4;              % nummber of elementray sensing action per cell.
ParamPreprocess.A        = 0;              % generate connectivity matrix A (yes/no).
ParamPreprocess.AWM      = 0;              % artificial world map (yes/no).
ParamPreprocess.RWM      = 1;              % real world map (yes/no).
% ParamPreprocess.MapFig   = [expTitle,'/maps/MapCoverage_HANDMADE.png'];
% ParamPreprocess.MapFig   = [expTitle,'/',mapDate,'/map/map.png'];
% ParamPreprocess.MapTxt   = [];
% ParamPreprocess.SenDom   = 'ComputeNew';
% ParamPreprocess.SenDom   = 'UseEarlier';
% ParamPreprocess.FileName = 'AASS_StudentArea_20140804_t1';
% ParamPreprocess.FileName = [expTitle,mapDate];
% ParamPreprocess.do_nTCell  = 1;
% ParamPreprocess.nTCellSize = 0.5;

map_size = MapSizes    
SensorRange.m = SensorRanges_m
FoV = FoVs
map_num = MapNums

load(['../Solution/MapSize_',num2str(map_size),'x',num2str(map_size),'cells/MapNum_',...
    num2str(map_num),'/SensorRange_',num2str(SensorRange.m),'cells/FoV_',...
    num2str(FoV),'deg/Preprocess_R8.mat']);            

                    
                    
FoVfaces.num        = ParamPreprocess.numFoVs;   % number of conf per cell.
n                   = numel(map);         % number of cells in the map.
f                   = FoVfaces.num;       % for simplification.
c                   = n*f;                % total conf.
cnt                 = ParamPreprocess.GraphCnt;  % 4 or 8 connected neighbours.
            
Fo = ones(FoVfaces.num*numel(map),3);
[map_row,map_col] = size(map);
n = numel(map);
o = cnt-1; % orientations

for i = 1:numel(OBS.ind)
    
    % Visibility
    V = [V(:,1:(OBS.ind(i)-1)*f) zeros(length(V(:,1)),f) ...
        V(:,(OBS.ind(i)-1)*f+1:end)];
    V = [V(1:(OBS.ind(i)-1),:); zeros(1,length(V(1,:))); ...
        V((OBS.ind(i)):end,:)];
end

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% % 1. PLOT THE MAP:
% ------------------------------------------------------------------------------
figure;
hold on;

for i = 1:map_row
    for j = 1:map_col
        % if obstacle/free
        if ~map(i,j) % occupied
            %col = [0.5 0.5 0.5];
            col = [0.1 0.1 0.1];
        elseif map(i,j) % free
            %col = 'y';
            col = [0.9 0.9 0.9];
        end
        plot(i-0.5,j-0.5,'S', 'color',col,'MarkerFaceColor',col,'MarkerSize',20);
    end
end
% xlim([-1 5])
% ylim([-1 5])
axis equal
grid on
% grid minor

pause

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% % 4. DRAW FIELD OF VIEWS:
% ------------------------------------------------------------------------------
% Draw FoVs and visible cells.
radius = SensorRange.cell;
for i = 1:c
    if any(Fo(i,:))
        %i
        %full(Fo(i,:))
        % Finding cell number.
        cell_num = round(i/f);
        if mod(i,f)/f < 0.5 && mod(i,f)/f > 0
            cell_num = cell_num+1;
        end
        % Finding configuration number.
        if mod(i,f) ~= 0
            cnf = mod(i,f);
        elseif mod(i,f) == 0
            cnf = f;
        end
        % index number
        [r,c] = ind2sub(size(map),cell_num);
        % Mark configuration
        plot(r-0.5,c-0.5,'o')
        
        % Draw FoV
        start_angle  = FoVfaces.ang(cnf,1);     % first angle
        sector_angle = FoV;                     % increment in first angle
        [xx, yy] = SecDraw(start_angle, sector_angle, radius);
        h1(1) = plot(xx+r-0.5,yy+c-0.5,'color','b'); hold on;
        
        % Visible cells
        vis = find(V(:,i));
        if vis
            [vis_r,vis_c] = ind2sub(size(map),vis);
            h1(2) = plot(vis_r-0.5,vis_c-0.5,'+r');
        end
        pause
        delete(h1(1));
        if vis
            delete(h1(2));
        end
    end
end



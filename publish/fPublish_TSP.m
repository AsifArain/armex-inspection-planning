function [ tPublishTSP ] = fPublish_TSP( map_env,...
                                         FoV,...
                                         FoVfaces,...
                                         V,...
                                         C,...
                                         SensorRange,...
                                         OBS,...
                                         TRoute,...
                                         ParamPublish )
%fPublish_AGPTSP_R1 publish graphical results of Art Gallery Problem.
% Date: 2014-01-09, Rev 1
% 
% n  := total number of cells in the map.
% nf := number of free cells in the map.
% f  := total number of sensing configurations per cell.
% 
% OUTPUTS:
% ----------------------------------------------------------------------------------------
% Graphical Display.
% TPostSolver: [scalar][sec] Total computation time for publishing the results.
% 
% INPUTS:
% ----------------------------------------------------------------------------------------
% map             : [],m][binary] l times m map. 1 := free cell, 0 := occupied cell.
% FoV             : [scalar] [degree] Field of view of a sensing configuration.
% FoVfaces        : [][] .num := number of faces/orientations of FoV.
%                        .alf := first/start angle for each FoV orientation.
%                        .bay := second/final angle for each FoV orientation.
%                        .ang := set of first/start and second/final angles
%                        of each FoV orientation. 
% V               : [nf,nf*f][binary] Visibility matrix. 1 := visibile, 0
% := not visibile. rows are cells and colums are conf. 
% C               : [nf*f,1][binary] Include/exclude status of all sensing
% conf. 1 := selected, 0 := not selected. 
% SensorRange     : [][] .m   := sensor range in meter,
%                        .cell:= sensor range in number of cells.
% OBS             : [][] .ind Indices of occupied cells.
%                        .sub Subscripts (row,col) of occupied cells.
% ParamPublish    : [][] Publish parameters.
%                        .MarkVis   := mark visibility (yes/no).
%                        .ConfSymbol:= draw sensing conf symbol (yes/no).
%                        .FoV       := draw field of view (yes/no).
%                        .fillFoV   := fill field of view with colors (yes/no).
%                        .Map       := publish map (yes/no).
%                        .CellSize  := cell size.
% 
% NOTES/UPDATES:
% ----------------------------------------------------------------------------------------
% 

%% CONSTANTS/SETTINGS:
% ----------------------------------------------------------------------------------------
tPublishTSPi = tic;
[map_row, map_col] = size(map_env);
% n = numel(map);
f = FoVfaces.num;
% c = n*f;
hold on;

%% MAP:
% ----------------------------------------------------------------------------------------
if ParamPublish.Map
    hold on;
    for iC = 1:map_row
        for i = 1:map_col
            % if obstacle/free
            if ~map_env(iC,i) % occupied
                col = [0 0 0];%[0.4 0.4 0.4];
            elseif map_env(iC,i) % free
                col = [0.9 0.9 0.9];
                %col = [1 1 1];
            end
            plot(iC-0.5,i-0.5,...
                'S',...
                'Color',col,...
                'MarkerFaceColor',col,...
                'MarkerSize',ParamPublish.CellSize);
        end
    end
   axis equal
end

%% TRAVELING ROUTE:
% ----------------------------------------------------------------------------------------
if ParamPublish.TravRoute
    
    for l = 1:numel(TRoute)-1
        iCELL = fix(TRoute(l)/f)+(~~mod(TRoute(l),f));
        iCONF = mod(TRoute(l),f)+(~mod(TRoute(l),f)*f);
        [iCELLr,iCELLc] = ind2sub(size(map_env),iCELL);

        if iCONF == 1 && f>1
            iCELLc = iCELLc+0.25;
        elseif iCONF == 2
            iCELLr = iCELLr-0.25;
        elseif iCONF == 3
            iCELLc = iCELLc-0.25;
        elseif iCONF == 4
            iCELLr = iCELLr+0.25;
        end

        jCELL = fix(TRoute(l+1)/f)+(~~mod(TRoute(l+1),f));
        jCONF = mod(TRoute(l+1),f)+(~mod(TRoute(l+1),f)*f);
        [jCELLr,jCELLc] = ind2sub(size(map_env),jCELL);

        if jCONF == 1 && f>1
            jCELLc = jCELLc+0.25;
        elseif jCONF == 2
            jCELLr = jCELLr-0.25;
        elseif jCONF == 3
            jCELLc = jCELLc-0.25;
        elseif jCONF == 4
            jCELLr = jCELLr+0.25;
        end

        x = [iCELLr jCELLr];
        y = [iCELLc jCELLc];
        plot(x-0.5,y-0.5,'-','color',[0.4 0.4 0.4],'LineWidth',1);
        %pause
    end
end

%% SENSING CONFIGURATIONS:
% ----------------------------------------------------------------------------------------

% Rearranging conf with obstacles.
for i = 1:numel(OBS.ind)    
    C = [C(1:(OBS.ind(i)-1)*f,:); zeros(f,length(C(1,:))); C((OBS.ind(i)-1)*f+1:end,:)];
end

iC  = find(C);                      % IDs of selected conf.
col = colormap(lines(numel(iC)));    % colors for the conf.


nCELL = fix(iC/f)+(~~mod(iC,f));   % cell num for each conf num.
nCONF = mod(iC,f)+(~mod(iC,f)*f);  % conf num within a cell.
[nCELLr,nCELLc] = ind2sub(size(map_env),nCELL); % row and col of cell num.

nCell_r = nCELLr; % keep another variable (row num).
nCell_c = nCELLc; % keep another variable (col num).

% To mark the conf according to its ID within a cell.
nCELLc = nCELLc+(0.25*(nCONF==1 & f>1));
nCELLr = nCELLr-(0.25*(nCONF==2));
nCELLc = nCELLc-(0.25*(nCONF == 3));
nCELLr = nCELLr+(0.25*(nCONF == 4));

dirXY = 4*[0 1; -1 0; 0 -1; 1 0]; % dX and dY for vectore (in accordance with FoVfaces.lgt)
if ParamPublish.ConfSymbol
    % mark the conf (accroding to ID with a cell).
    plot(nCELLr-0.5,nCELLc-0.5,'ob','MarkerFaceColor','b','MarkerSize',5); 
    for i = 1:length(iC)
        % Draw FoV
        %start_angle  = FoVfaces.ang(nCONF(i),1);     % first/start angle.
        %sector_angle = FoV;                          % increment in the first/start angle.
        %MaxPts       = 2;                            % maximum points to plot the FoV.
        %[xx,yy]      = SecDraw(start_angle,sector_angle,0.5,MaxPts); % x,y points.
        %plot(xx+nCELLr(i)-0.5,yy+nCELLc(i)-0.5,'color','k','LineWidth',1); %2.5
%         h = quiver(nCELLr(i)-0.5,nCELLc(i)-0.5,...
%             dirXY(nCONF(i),1),dirXY(nCONF(i),2),'LineWidth',1,'color','b');
        %adjust_quiver_arrowhead_size(h,5.5);
%         adjust_quiver_arrowhead_size(h,1);
        
    end
end

if ParamPublish.FoV
    % mark the conf (at center of the cell).
    plot(nCell_r-0.5,nCell_c-0.5,'ok','MarkerFaceColor','k','MarkerSize',3); 
    for i = 1:length(iC)
        start_angle  = FoVfaces.ang(nCONF(i),1);     % first/start angle.
        sector_angle = FoV;                        % increment in the first/start angle.
        MaxPts       = 100;                          % maximum points to plot the FoV.
        [xx,yy]      = SecDraw(start_angle,sector_angle,SensorRange.cell,MaxPts);
        %plot(xx+nCell_r(i)-0.5,yy+nCell_c(i)-0.5,'color',col(i,:),'LineWidth',1);
        plot(xx+nCell_r(i)-0.5,yy+nCell_c(i)-0.5,'color',[0 0 1],'LineWidth',1);        
    end
end

%%%%%%% To be adjusted %%%%%%
if ParamPublish.FoV
    % mark the conf (at center of the cell).
    plot(nCell_r-0.5,nCell_c-0.5,'ok','MarkerFaceColor','k','MarkerSize',3); 
    for i = 1:4
        start_angle  = FoVfaces.ang(i,1)          % first/start angle.
        sector_angle = 5;                         % increment in the first/start angle.
        MaxPts       = 100;                       % maximum points to plot the FoV.
        [xx,yy]      = SecDraw(start_angle,sector_angle,5,MaxPts); % x,y points.
        %plot(xx+nCell_r(i)-0.5,yy+nCell_c(i)-0.5,'color',col(i,:),'LineWidth',1);
        plot(xx+nCell_r(i)-0.5,yy+nCell_c(i)-0.5,'color',[0 0 1],'LineWidth',1);        
    end
end


if ParamPublish.fillFoV
    for i = 1:length(iC)
        % Draw FoV
        start_angle  = FoVfaces.ang(nCONF(i),1); % first/start angle.
        sector_angle = FoV;                      % increment in the first/start angle.
        MaxPts       = 100;                      % maximum points to plot the FoV.
        % x,y points.
        [xx,yy]      = SecDraw(start_angle,sector_angle,SensorRange.cell,MaxPts);
        % fill the FoV.
        fill(xx+nCell_r(i)-0.5,yy+nCell_c(i)-0.5,col(i,:)); alpha(0.2);          
    end
end

if ParamPublish.MarkVis
    
    % Rearranging visibility matrix with obstacles.
    for iC = 1:numel(OBS.ind)
        V = [V(:,1:(OBS.ind(iC)-1)*f) zeros(length(V(:,1)),f) V(:,(OBS.ind(iC)-1)*f+1:end)];    
        V = [V(1:(OBS.ind(iC)-1),:); zeros(1,length(V(1,:))); V((OBS.ind(iC)):end,:)];
    end
    
    % Visible cells
    numPosC = (ceil(sqrt(numel(find(C))))); % size of a square matrix according to num of conf.
    PosC_xy = 1/(numPosC+1):1/(numPosC+1):1-(1/(numPosC+1)); % conf position num.
    [p,q] = meshgrid(PosC_xy,PosC_xy); % 
    PosC = [p(:) q(:)]; % 

    a = 0;
    for iC = (find(C))'
        vis = find(V(:,iC));
        a = a+1;
        [vis_r,vis_c] = ind2sub(size(map_env),vis);
        %plot(vis_r-0.5,vis_c-0.5,'+','Color',col(a,:));
        plot(vis_r-PosC(a,1),vis_c-PosC(a,2),'+','Color',col(a,:));
    end
    
end


%% PUBLISH SETTING
% ----------------------------------------------------------------------------------------

xlim([-1 length(map_env(:,1))+1])
ylim([-1 length(map_env(1,:))+1])
set(gca, ...
  'GridLineStyle','-',....
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'on'      , ...
  'XGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'XTick'       , [], ...  
  'YTick'       , [], ...
  'LineWidth'   , 1         );

% axis off

%% PRINT
% ----------------------------------------------------------------------------------------

set(gcf,'PaperPositionMode','auto','PaperSize',[11.7 8.3]);

if ParamPublish.PrintPDF
    %print -painters -dpdf -r1049 Name.pdf
    eval(['print -painters -dpdf -r1049 ',ParamPublish.PrintName,'.pdf']);
end

if ParamPublish.PrintEPS
    %print -depsc Name.eps
    eval(['print -depsc ',ParamPublish.PrintName,'.eps']);
end

if ParamPublish.PrintTIKZ
    %
end


%% 
% Total computation time taken by the fResults.
tPublishTSP =  1e-4*(round(toc(tPublishTSPi)*1e4));


end

function [xx, yy] = SecDraw( start_angle,sector_angle,radius,MaxPts )
%  Input:
%            start_angle  - starting angle in degrees, [0 is default]
%                           NOTE: this can be positive or negative
%            sector_angle - central angle of sector in degrees
%                           NOTE: if negative, plotting will be clockwise
%                                 direction.                                 
%            radius       - radius of sector (default 1)
%
% Output:                           
%            returns the corner points (xx, yy) of a patch for a sector,
%            which can plotted using the patch(...) function
%            NOTE:
%              1. if less than 2 or NO output variables are given, 
%                 secdraw(...) draws the sector on the current (or a new)
%                 figure
%              2. if one argument is given, 
%                 the argument is taken as the central angle 

theta0 = pi*0/180;              % offset angle in radians
theta1 = pi*start_angle/180;    % starting angle in radians
theta  = pi*sector_angle/180;   % central angle in radians

angle_ratio = theta/(2*pi);     % ratio of the sector to the complete circle
rho = abs(radius);              % abs(...) can be removed
% MaxPts = 100;                   % maximum number of points in a whole circle.
                                % if set to small # (e.g 10) the resolution
                                % of the curve will be poor. 
                                % generally values greater than 50
                                % will give (very) good resolution

n = abs( ceil( MaxPts*angle_ratio ) );
r = [ 0; rho*ones(n+1,1); 0 ];
theta = theta0 + theta1 + [ 0; angle_ratio*(0:n)'/n; 0 ]*2*pi;

% output
[xx,yy] = pol2cart(theta,r);

% plot if not enough output variable are given
if nargout < 2
    hh=patch(xx, yy, 'b','EdgeColor','b','EraseMode','none');
end

end

% ----------------------- End of the document ----------------------------------


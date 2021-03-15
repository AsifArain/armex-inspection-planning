function [ conf_crossAngles,...
           conf_crossAngles_G ] = fXReEstimatedHotspot( map_env,...
                                                        confKept,...
                                                        hotspot_current,...
                                                        NumOfAllowedConf )
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
% ----------------------------------------------------------------------------------------
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
% visualize    = 0;




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
% CROSS ANGLES
%=======================================================================


fprintf('Cross angles...\n')

[ Hot2Conf_Angles,...
  conf_crossAngles,...
  conf_crossAngles_G ] = fX( confKept,...
                             hotspots,...
                             NumOfAllowedConf,...
                             map_env );



%% Total computation time.
% ----------------------------------------------------------------------------------------
tPreprocess = 1e-4*(round(toc(tPreprocessi)*1e4));

disp(['tPreprocess: ',num2str(tPreprocess),' sec']);


end





function [ Hot2Conf_Angles,...
           conf_crossAngles,...
           conf_crossAngles_G] = fX( confKept,...
                                     hotspots,...
                                     NumOfAllowedConf,...
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
        conf_crossAngles(i,j) =...
            rad2deg(angleAbsDiff(Hot2Conf_Angles(1,i),Hot2Conf_Angles(1,j)));
    end
end


%-- cross angles gain 
disp('----- cross gain')

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
    
end


% guss_gain = cell(1,numel(gain_meu));
guss_gain = zeros(size(conf_crossAngles));
for i = 1:numel(gain_meu)
    %guss_gain = guss_gain +...
    %    gaussmf(angles,[gain_sigma,gain_meu(i)]);
    %guss_gain{i} = gaussmf(conf_crossAngles,[gain_sigma,gain_meu(i)]);
    guss_gain = max(guss_gain,gaussmf(conf_crossAngles,[gain_sigma,gain_meu(i)]));
    %plot(angles,guss_gain{i},'-g','LineWidth',3)
end
% conf_crossAngles_G = max(vec2mat(cell2mat(guss_gain),numel(conf_crossAngles)),[],1);
conf_crossAngles_G = guss_gain;


%{
% ------ cross angles gain 1
disp('----- gain')
% gain_sigma = 10; % standard deviation

if     NumOfAllowedConf == 2
    gain_angles = 90;
    gain_sigma = 10; % standard deviation
    disp('Cross angles for two conf.')
elseif NumOfAllowedConf == 3
    gain_angles = [60,120];
    gain_sigma = 10; % standard deviation
    disp('Cross angles for three conf.')
elseif NumOfAllowedConf == 4
    gain_angles = [45,90,135];
    gain_sigma = 5; % standard deviation
    disp('Cross angles for four conf.')
elseif NumOfAllowedConf == 5
    gain_angles = [36,72,108,144];
    disp('Cross angles for five conf.')
    gain_sigma = 5; % standard deviation
else
    warning('Gain for cross angles is unknown')
end

gain_angles
conf_crossAngles_G = zeros(size(Hot2Conf_Angles,2));
for gain_meu = gain_angles % expected value or center
    gain_meu
    conf_crossAngles_G = conf_crossAngles_G +...
        gaussmf(conf_crossAngles,[gain_sigma,gain_meu]);
end
%}

%FIXME
%DEBUG
%{
xAngles_array = conf_crossAngles(:);
gAngles_array = conf_crossAngles_G(:);
[xAngles_val,xAngles_ind] = sort(xAngles_array);

xAngles = xAngles_array(xAngles_ind);
gAngles = gAngles_array(xAngles_ind);

figure; hold on;
plot(xAngles,'-r');
plot(180*gAngles,'-g');
pause
%}



%----- OLD ------

% % ------ cross angles gain 1
% disp('----- gain')
% % gain_sigma = deg2rad(10); % standard deviation
% % gain_meu   = deg2rad(60); % expected value or center
% gain_sigma = 10; % standard deviation
% gain_meu   = 60; % expected value or center
% conf_crossAngles_G1 = gaussmf(conf_crossAngles,[gain_sigma,gain_meu]);
% 
% 
% % gain_sigma = deg2rad(10); % standard deviation
% % gain_meu   = deg2rad(120); % expected value or center
% gain_sigma = 10; % standard deviation
% gain_meu   = 120; % expected value or center
% conf_crossAngles_G2 = gaussmf(conf_crossAngles,[gain_sigma,gain_meu]);
% 
% conf_crossAngles_G = conf_crossAngles_G1+conf_crossAngles_G2;



%{ 
NumOfAllowedConf = 5;
angles = 0:1:180; % in degrees
guss_gain = zeros(1,size(angles,2));
if     NumOfAllowedConf == 2
    gain_angles = 90;
    gain_sigma = 10; % standard deviation
    disp('Cross angles for two conf.')
elseif NumOfAllowedConf == 3
    gain_angles = [60,120];
    gain_sigma = 10; % standard deviation
    disp('Cross angles for three conf.')
elseif NumOfAllowedConf == 4
    gain_angles = [45,90,135];
    gain_sigma = 5; % standard deviation
    disp('Cross angles for four conf.')
elseif NumOfAllowedConf == 5
    gain_angles = [36,72,108,144];
    gain_sigma = 5; % standard deviation
    disp('Cross angles for five conf.')
else
    warning('Gain for cross angles is unknown')
end

for gain_meu = gain_angles % expected value or center
    gain_meu
    guss_gain = guss_gain +...
        gaussmf(angles,[gain_sigma,gain_meu]);
end

figure; plot(guss_gain)




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

% ----------------------------- End of the document --------------------------------------

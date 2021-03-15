function [ conf_crossAngles,...
           conf_crossAngles_G,...
           Angles_FromHc2Conf] = fXReEstimatedHc( hc,...
                                                  map_env,...
                                                  confKept,...
                                                  oConfOrientations,...
                                                  NumOfAllowedConf,...
                                                  para_ )
% fPreProcess_R4 genereates required parameters for solving a global coverage problem.
% Date: 2014-01-19, Rev 4
% 
% n  := total number of cells in the map.
% nf := number of free cells in the map.
% f  := total number of sensing configurations per cell.
% 
% OUTPUTS:
% 
% INPUTS:
% ----------------------------------------------------------------------------------------

%
% NOTES/UPDATES:
% ----------------------------------------------------------------------------------------
% 
% 

tPreprocessi = tic; % initialize pre computational time.
% visualize    = 0;

%=======================================================================
% CONF ANGLES TOWARDS HC
%=======================================================================

fprintf('Cross angles...\n')



% Hot2Conf_Angles = zeros(1,numel(confKept));

confKeptCell_ind = confKept;
[confKeptCell_r,confKeptCell_c] = ind2sub(size(map_env),confKeptCell_ind);
sensors = [confKeptCell_r',confKeptCell_c'];

% [hotspots_r,hotspots_c] = ind2sub(size(map_env),hotspots);
    
%Hot2Conf_Angles(i,:) = rad2deg(angle2Points([hotspots_r(i),hotspots_c(i)],sensors));
Angles_FromHc2Conf = angle2Points(hc,sensors);


%-------
% DO NOT USE oConfOrientations
% it will create problem for fixed confs.
%-------
% % -- orientations of confs kept
% ConfOrientations = oConfOrientations(confKept);
% % -- angles from hotspot to confs kept
% Angles_FromHc2Conf2 = deg2rad(wrapTo360(ConfOrientations+180));
% [rad2deg(Angles_FromHc2Conf(:)),rad2deg(Angles_FromHc2Conf2(:)),...
%     rad2deg(angleAbsDiff(Angles_FromHc2Conf2(:),Angles_FromHc2Conf(:)))]
% pause


%=======================================================================
% CROSS ANGLES
%=======================================================================

% ----- cross angles 1
disp('----- cross angles')
% conf_crossAngles = zeros(size(Angles_FromHc2Conf,2));
conf_crossAngles = zeros(numel(confKept));

for i = 1:numel(confKept)
    for j = 1:numel(confKept)
        conf_crossAngles(i,j) = ...
            rad2deg(angleAbsDiff(Angles_FromHc2Conf(i),Angles_FromHc2Conf(j)));
    end
end


 
%=======================================================================
% CROSS ANGLES GAIN
%=======================================================================


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
else
    error('Unknown cross angle gains')
    
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
conf_crossAngles_G = zeros(size(conf_crossAngles));
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



%% Total computation time.
% ----------------------------------------------------------------------------------------
tPreprocess = 1e-4*(round(toc(tPreprocessi)*1e4));

disp(['tPreprocess: ',num2str(tPreprocess),' sec']);



end

% ----------------------------- End of the document --------------------------------------

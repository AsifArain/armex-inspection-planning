function [ conf_crossAngles,...
           conf_crossAngles_G,...
           Angles_FromHc2Conf ] = fXHc1t( hc_this_sub,...
                                          map_env,...
                                          confKept,...
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

%***********************************************************************
%
% CONF ANGLES TOWARDS HC
%
%***********************************************************************
fprintf('---- angles towards Hc...\n')

confKeptCell_ind = confKept;
[confKeptCell_r,confKeptCell_c] = ind2sub(size(map_env),confKeptCell_ind);
sensor_placements = [confKeptCell_r',confKeptCell_c'];
Angles_FromHc2Conf = angle2Points(hc_this_sub,sensor_placements);



%***********************************************************************
%
% CROSS ANGLES
%
%***********************************************************************
disp('---- cross angles')

conf_crossAngles = zeros(numel(confKept));
for i = 1:numel(confKept)
    for j = 1:numel(confKept)
        conf_crossAngles(i,j) = ...
            rad2deg(angleAbsDiff(Angles_FromHc2Conf(i),Angles_FromHc2Conf(j)));
    end
end

 
%***********************************************************************
%
% CROSS ANGLES GAIN
%
%***********************************************************************
disp('---- cross angles gain')

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


guss_gain = zeros(size(conf_crossAngles));
for i = 1:numel(gain_meu)
    guss_gain = ...
        max(guss_gain,gaussmf(conf_crossAngles,[gain_sigma,gain_meu(i)]));
end
conf_crossAngles_G = guss_gain;


%% Total computation time.
% ----------------------------------------------------------------------------------------
tPreprocess = 1e-4*(round(toc(tPreprocessi)*1e4));

% disp(['---- tPreprocess: ',num2str(tPreprocess),' sec']);
fprintf('---- tPreprocess: %.4f sec\n\n',tPreprocess);

end

% ----------------------------- End of the document --------------------------------------

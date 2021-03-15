function fPublish_RWL1Iterations( dataPublish,para_,dir_ )
%fPublishRWL1Iterations publishs graphical results of Art Gallery Problem.
% Date: 2019-03-04
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
% V               : [nf,nf*f][binary] Visibility matrix. 1 := visibile,
%                                                        0 := not visibile.  
%                             rows are cells and colums are conf.  
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
% ----------------------------------
%FoVfaces = dataPublish.FoVfaces;
%f = FoVfaces.num;
% ----------------------------------


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                         PUBLISHING PARAMETERS
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
paraPub.zoom_in_environment = 5750;
paraPub.cell_marker_size = 10; %4; %10.25 %4;%4.85;
paraPub.color_scale = 'linear'; %
%beam_color = [255,204,153]/256; % [204,229,255]/256
paraPub.marker_size_cand_conf = 2.5;
paraPub.marker_size_hotspot = 5; %7;
%marker_size_weak_conf = 4;
%marker_size_cond_conf = 5;
paraPub.marker_size_conf_position = 3; %5; % icra-submission size: 3
paraPub.line_width_fov = 1; %2.5
paraPub.symbolic_fov_length = 2.0; %3;
paraPub.line_width_symbolic_fov = 0.5; %1.25;

paraPub.arrow_length = 0.5;
paraPub.line_width_arrow = 2; %0.25
paraPub.arrow_head_size = 20%20;
paraPub.marker_size_arrow = 1; %2

paraPub.font_size_conf_num = 6;%16
paraPub.font_size_conc_scale = 10;
paraPub.conf_num_dist = 0.5; % 0.3
    



load('CIterativeValues.mat');

dir_.logs = dir_.gdplan_logs;    
load([dir_.gdplan_spp,'preprocess_detection.mat'],...
        'oV',...
        'SensorRange',...
        'FoVfaces',...
        'confKept');
%load([dir_.DetectionResults,'spp_detection_solution.mat'],'C')

% Conf = zeros(size(oV,2),1);
% Conf(confKept) = C;
map_env = dataPublish.map_env;
% ----------------------------------
FoVfaces = dataPublish.FoVfaces;
f = FoVfaces.num;



%-----------------------------------------------------------------------------------
%
%           ENVIRONMENT MAP WITH CANDIDATE CONFIGURATIONS
%
%-----------------------------------------------------------------------------------

%-- environment map
fPublishEnv(map_env,para_,paraPub)

%-- confs.
C = ones(size(CIterativeValues,1),1);
Conf = zeros(size(oV,2),1);
Conf(confKept) = C;
C  = Conf;
iC = 1:numel(C);
nCELL = fix(iC/f)+(~~mod(iC,f));   % cell num for each conf num.
nCONF = mod(iC,f)+(~mod(iC,f)*f);  % conf num within a cell.
[nCELLr,nCELLc] = ind2sub(size(map_env),nCELL); % row and col of cell num.

%-- plot poses
for i = 1:length(iC)
    this_color = [1.0,1.0,1.0];    
    if map_env(nCELLr(i),nCELLc(i))
        quiver(nCELLr(i),nCELLc(i),...
               paraPub.arrow_length*cosd(FoVfaces.lgt(nCONF(i))),...
               paraPub.arrow_length*sind(FoVfaces.lgt(nCONF(i))),...
               'LineWidth',paraPub.line_width_arrow,...
               'color',this_color,...
               'MaxHeadSize',paraPub.arrow_head_size,...
               'Marker','o',...
               'MarkerSize',paraPub.marker_size_arrow);
    end
end

%-- save graphics
file_name = sprintf('%s-rwl1-candidate-map',para_.ExperimentTitle);
export_fig([dir_.gdplan_spp,file_name], '-pdf','-r500')
export_fig([dir_.gdplan_spp,file_name], '-png','-r500')



%-----------------------------------------------------------------------------------
%
%           RWL1 ITERATIONS
%
%-----------------------------------------------------------------------------------

for itr = 1:size(CIterativeValues,2)
    
    %-- environment map
    fPublishEnv(map_env,para_,paraPub)

    %-- confs.
    C = CIterativeValues(:,itr);
    Conf = zeros(size(oV,2),1);
    Conf(confKept) = C;
    C  = Conf;
    iC = 1:numel(C);
    nCELL = fix(iC/f)+(~~mod(iC,f));   % cell num for each conf num.
    nCONF = mod(iC,f)+(~mod(iC,f)*f);  % conf num within a cell.
    [nCELLr,nCELLc] = ind2sub(size(map_env),nCELL); % row and col of cell num.

    %-- plot poses
    conf_colormap = flipud(hot(101)); 
    for i = 1:length(iC)
        this_color = conf_colormap(round(100*C(i))+1,:);
        if map_env(nCELLr(i),nCELLc(i))
            quiver(nCELLr(i),nCELLc(i),...
                   paraPub.arrow_length*cosd(FoVfaces.lgt(nCONF(i))),...
                   paraPub.arrow_length*sind(FoVfaces.lgt(nCONF(i))),...
                   'LineWidth',paraPub.line_width_arrow,...
                   'color',this_color,...
                   'MaxHeadSize',paraPub.arrow_head_size,...
                   'Marker','o',...
                   'MarkerSize',paraPub.marker_size_arrow);
        end
    end

    %-- save graphics
    if itr < size(CIterativeValues,2)
        file_name = sprintf('%s-rwl1-iteration%02d',para_.ExperimentTitle,itr);
    else
        file_name = sprintf('%s-rwl1-comb',para_.ExperimentTitle);
    end
    export_fig([dir_.gdplan_spp,file_name], '-pdf','-r500')
    export_fig([dir_.gdplan_spp,file_name], '-png','-r500')

end





end

function fPublishEnv(map_env,para_,paraPub)

h = figure('name',para_.PublishFigureTitle);
imshow(map_env','InitialMagnification',paraPub.zoom_in_environment); 
hold on;
set(gca,'YDir','normal')
colormap('gray'); % not go back to gray scale.
%caxis([-5 2]) % for resulting gdm
%caxis([-12 4]) % for sensor placement solution
% caxis([-0.25 1]) % for sensor placement solution
caxis([-1 1.25])
% caxis([-2 2])
%caxis([-10 1.5]) % for gdm insets

end

% ----------------------- End of the document ----------------------------------

function fPublish( dataPublish,para_,dir_ )
%fPublish_AGPTSP_R1 publish graphical results of Art Gallery Problem.
% Date: 2014-01-09, Rev 1
% 
% n  := total number of cells in the map.
% nf := number of free cells in the map.
% f  := total number of sensing configurations per cell.
% 
% OUTPUTS:
% ---------------------------------------------------------------------------------------------------------------------------
% Graphical Display.
% TPostSolver: [scalar][sec] Total computation time for publishing the results.
% 
% INPUTS:
% ---------------------------------------------------------------------------------------------------------------------------
% map             : [],m][binary] l times m map. 1 := free cell, 0 := occupied cell.
% FoV             : [scalar] [degree] Field of view of a sensing configuration.
% FoVfaces        : [][] .num := number of faces/orientations of FoV.
%                        .alf := first/start angle for each FoV orientation.
%                        .bay := second/final angle for each FoV orientation.
%                        .ang := set of first/start and second/final angles of each FoV orientation.
% V               : [nf,nf*f][binary] Visibility matrix. 1 := visibile, 0 := not visibile. 
%                             rows are cells and colums are conf.  
% C               : [nf*f,1][binary] Include/exclude status of all sensing conf. 1 := selected, 0 := not selected.
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
% ---------------------------------------------------------------------------------------------------------------------------
% 
% ----------------------------------
FoVfaces = dataPublish.FoVfaces;
f = FoVfaces.num;
% ----------------------------------


% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                         PUBLISHING PARAMETERS
% 
% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
zoom_in_environment = 350; %800
cell_marker_size = 2.75; %4.85 
color_scale = 'linear'; %
%color_scale = 'nonlinear'; % main fig.
beam_color = [255,204,153]/256; % [204,229,255]/256
marker_size_cand_conf = 2.5;
marker_size_hotspot = 7;
marker_size_weak_conf = 4;
marker_size_cond_conf = 5;
marker_size_conf_position = 2; %5; % icra-submission size: 3
line_width_fov = 1; %2.5
symbolic_fov_length = 3;
line_width_symbolic_fov = 0.5; %1.25;

arrow_length = 3;
line_width_arrow = 0.5; %0.25
arrow_head_size = 20;
marker_size_arrow = 1; %2

font_size_conf_num = 6;%16
conf_num_dist = 0.5; % 0.3

% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             ENVIRONMENT MAP:
% 
% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
map_env = dlmread([dir_.env,para_.EnvironmentFigFile]);
map_env = 1-map_env;
[map_row,map_col] = size(map_env);

figure; imshow(map_env','InitialMagnification',zoom_in_environment); hold on;
set(gca,'YDir','normal')
colormap('gray'); % not go back to gray scale.
%caxis([-5 2]) % for resulting gdm
%caxis([-12 4]) % for sensor placement solution
caxis([-0.25 1]) % for sensor placement solution
% caxis([-2 2])
%caxis([-10 1.5]) % for gdm insets

% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             GAS CONCENTRATION MAP
% 
% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishGDM
    
    % -- reconstruction        
    load([dir_.results,para_.ReconstructionFile],...
        'mean_map',...
        'cell_coords_x',...
        'cell_coords_y');
    map_gd = mean_map;
    
    
    % ------------- select color map for gdm --------------
    % colormap hot; 
    % colormap gray; 
    % colormap(flipud(colormap('hot')));
    gdm_colormap = colormap(flipud(colormap('hot'))); % for main fig. 
    %gdm_colormap = colormap(flipud(colormap('bone')));
    %gdm_colormap = colormap(flipud(colormap('copper')));
    %gdm_colormap = colormap(flipud(colormap('pink')));
    %gdm_colormap = colormap(flipud(colormap('cool')));
    %gdm_colormap = colormap(flipud(colormap('spring')));
    %gdm_colormap = colormap(flipud(colormap('summer'))); % for insets
    %gdm_colormap = colormap(flipud(colormap('autumn')));
    %gdm_colormap = colormap('summer');
    %gdm_colormap = colormap('autumn');
    %colormap('gray'); % not go back to gray scale.
    %caxis([-5 2]) % for resulting gdm
    %caxis([-12 4]) % for sensor placement solution
    %caxis([-2 2]) % for sensor placement solution
    %caxis([-10 1.5]) % for gdm insets
    % load gdm_colormap.mat
    % delta = max(map_gd(:))/(size(gdm_colormap,1)-1);
    %load publish/gdm_colormap.mat
    % -----------------------------------------------------
    
    
    % ------ color step size --------------------------------
    %max_concentration = 8000; % inset gdm
    max_concentration = max(map_gd(:)) % main fig.
    max_concentration = 2.5e+04
    %max_concentration = 4.7135e+10
    %pause
    delta = max_concentration/(size(gdm_colormap,1)-1);    
    %delta = max(map_gd(:))/(size(gdm_colormap,1)-1);
    %delta = 500/(size(gdm_colormap,1)-1); % lets limit to 500ppm
    
    % lets plot for the scale of 500 ppm.
    %delta = 3000/(size(gdm_colormap,1)-1);
    %max(map_gd(:))
    %pause
    %delta = 1.78e+04/(size(gdm_colormap,1)-1)
    
    
    % ------------------- nonlinear color code ---------------------  
    if strcmp(color_scale,'nonlinear')
        %b = (exp(256/1))/255;
        %a = 256 /(log(b*256));

        x = 1:1:256;
        %epsilon = 1/(exp(1)-1); % standard
        %epsilon = 1/(exp(0.05)-1); % slower
        epsilon = 1/(exp(0.10)-1); % slower
        %epsilon = 1/(exp(0.15)-1); % slower
        %epsilon = 1/(exp(0.20)-1); % slower
        %epsilon = 1/(exp(0.90)-1); % slower
        %epsilon = 1/(exp(1.5)-1); % standard
        normlogsum_penalty = log(1+(abs(x)/epsilon));
        %normlogsum_penalty = 1.05.^(1+(abs(x)/epsilon));
        normlogsum_penalty = normlogsum_penalty/(normlogsum_penalty(x==max(x))); % normalized
        normlogsum_penalty = round(256*normlogsum_penalty);
        normlogsum_penalty(1)
        normlogsum_penalty(end)
        %figure; plot(x,normlogsum_penalty); % to test

        % publish nonlinear steps
        linear_percentage    = 10:10:100;
        linear_percentage_ind  = round((linear_percentage/100)*numel(x));
        nonlinear_percentage_ind = normlogsum_penalty(linear_percentage_ind);
        nonlinear_percentage = (nonlinear_percentage_ind/numel(x))*100;
        nonlinear_percentage = round(nonlinear_percentage);
        %disp('nonlinear_percentage = ',num2str(nonlinear_percentage))
        disp_stng = ['nonlinear color code scale is: ',num2str(nonlinear_percentage)];
        %disp(disp_stng)
        disp(disp_stng)
    end
    % ----------------------------------------------------------------
    
    
    for i = 1:numel(cell_coords_y)%map_row
        for j = 1:numel(cell_coords_x)%map_col
            % if positive concentration
            
            if map_gd(i,j)>=0
                
                % ---------------------- useful trash -------------------
                %col = gdm_colormap(round(map_gd(i,j)/delta)+1,:);
                
                % --- if concentration is greater than 1000 ppm, its still 1000 ppm
                %col = gdm_colormap( round(min(map_gd(i,j),3000)/delta)+1, :);
                %col = gdm_colormap( round(min(map_gd(i,j),1.78e+04)/delta)+1, :);
                %col = gdm_colormap( round(min(map_gd(i,j),inf)/delta)+1, :);
                % -------------------------------------------------------
                
                % ----------------------------------------------------
                % linear color code
                % ----------------------------------------------------
                if strcmp(color_scale,'linear')
                % {
                linear_color_number = round( min(map_gd(i,j),max_concentration)/delta)+1;
                col = gdm_colormap( linear_color_number,:);                
                %}
                end
                % ----------------------------------------------------
                
                % ----------------------------------------------------
                % nonlinear color code
                % ----------------------------------------------------
                if strcmp(color_scale,'nonlinear')
                % {
                linear_color_number = round( min(map_gd(i,j),max_concentration)/delta)+1;
                %nonlinear_color_number = round( a*log(b*linear_color_number) );
                %nonlinear_color_number = linear_color_number;
                nonlinear_color_number = normlogsum_penalty(linear_color_number);
                %if nonlinear_color_number <1
                %    nonlinear_color_number = 1;
                %end
                col = gdm_colormap( nonlinear_color_number,:);                
                %}
                end
                % ----------------------------------------------------
                
                %----- to avoid black borders of earch cell -----
                % {
                if sum(col) == 3
                    col = col-(1e-10);
                end
                %}
                % ----- plot ------------
                %plot(j,i,'s','Color',col,'MarkerFaceColor',col,'MarkerSize',4);
                %plot(i,j,'s','Color',col,'MarkerFaceColor',col,'MarkerSize',4);
                %plot(cell_coords_x(i),cell_coords_y(j),'s','Color',col,'MarkerFaceColor',col,'MarkerSize',4);
                
                %plot(cell_coords_y(i)+0.5,cell_coords_x(j)+0.5,...
                %    's','Color',col,'MarkerFaceColor',col,'MarkerSize',5);
                
                %plot(cell_coords_x(j)+0.5,cell_coords_y(i)+0.5,...
                %    's','Color',col,'MarkerFaceColor',col,'MarkerSize',4.75);
                plot(cell_coords_x(j)+0.5,cell_coords_y(i)+0.5,...
                    's',...
                    'MarkerEdgeColor',col,...
                    'MarkerFaceColor',col,...
                    'MarkerSize',cell_marker_size);
            end
        end
    end   
end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             GROUND TRUTH
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishGroundTruth
    
    % -- ground truth
    gt = dlmread([dir_.logs,para_.GroundTruthFile], '', 1,0);
    gt_z = gt(gt(:,3)==para_.GroundTruthZSlice,[1,2,4]);
        
    
    % ------------- select color map for gdm --------------
    % colormap hot; 
    % colormap gray; 
    % colormap(flipud(colormap('hot')));
    gdm_colormap = colormap(flipud(colormap('hot'))); % for main fig. 
    %gdm_colormap = colormap(flipud(colormap('bone')));
    %gdm_colormap = colormap(flipud(colormap('copper')));
    %gdm_colormap = colormap(flipud(colormap('pink')));
    %gdm_colormap = colormap(flipud(colormap('cool')));
    %gdm_colormap = colormap(flipud(colormap('spring')));
    %gdm_colormap = colormap(flipud(colormap('summer'))); % for insets
    %gdm_colormap = colormap(flipud(colormap('autumn')));
    %gdm_colormap = colormap('summer');
    %gdm_colormap = colormap('autumn');
    colormap('gray'); % not go back to gray scale.
    %caxis([-5 2]) % for resulting gdm
    %caxis([-12 4]) % for sensor placement solution
    caxis([-2 2]) % for sensor placement solution
    %caxis([-10 1.5]) % for gdm insets
    % load gdm_colormap.mat
    % delta = max(map_gd(:))/(size(gdm_colormap,1)-1);
    %load publish/gdm_colormap.mat
    % -----------------------------------------------------
    
    
    % ------ color step size --------------------------------
    %max_concentration = 8000; % inset gdm
    
    %max_concentration = max(map_gd(:)); % main fig.
    %max_concentration = max(gt_z(:,3))
    
    %gt_z(:,3) = gt_z(:,3)/100000;    
    max_concentration = max(gt_z(:,3))
    %max_concentration = 632.8936
    
    %pause
    delta = max_concentration/(size(gdm_colormap,1)-1);    
    %delta = max(map_gd(:))/(size(gdm_colormap,1)-1);
    %delta = 500/(size(gdm_colormap,1)-1); % lets limit to 500ppm
    
    % lets plot for the scale of 500 ppm.
    %delta = 3000/(size(gdm_colormap,1)-1);
    %max(map_gd(:))
    %pause
    %delta = 1.78e+04/(size(gdm_colormap,1)-1)
    
    
    % ------------------- nonlinear color code ---------------------  
    if strcmp(color_scale,'nonlinear')
        %b = (exp(256/1))/255;
        %a = 256 /(log(b*256));

        x = 1:1:256;
        %epsilon = 1/(exp(1)-1); % standard
        %epsilon = 1/(exp(0.05)-1); % slower
        epsilon = 1/(exp(0.10)-1); % slower
        %epsilon = 1/(exp(0.15)-1); % slower
        %epsilon = 1/(exp(0.20)-1); % slower
        %epsilon = 1/(exp(0.90)-1); % slower
        %epsilon = 1/(exp(1.5)-1); % standard
        normlogsum_penalty = log(1+(abs(x)/epsilon));
        %normlogsum_penalty = 1.05.^(1+(abs(x)/epsilon));
        normlogsum_penalty = normlogsum_penalty/(normlogsum_penalty(x==max(x))); % normalized
        normlogsum_penalty = round(256*normlogsum_penalty);
        normlogsum_penalty(1)
        normlogsum_penalty(end)
        %figure; plot(x,normlogsum_penalty); % to test

        % publish nonlinear steps
        linear_percentage    = 10:10:100;
        linear_percentage_ind  = round((linear_percentage/100)*numel(x));
        nonlinear_percentage_ind = normlogsum_penalty(linear_percentage_ind);
        nonlinear_percentage = (nonlinear_percentage_ind/numel(x))*100;
        nonlinear_percentage = round(nonlinear_percentage);
        %disp('nonlinear_percentage = ',num2str(nonlinear_percentage))
        disp_stng = ['nonlinear color code scale is: ',num2str(nonlinear_percentage)];
        %disp(disp_stng)
        disp(disp_stng)
    end
    % ----------------------------------------------------------------
        
    for i = 1:size(gt_z,1)
        %i
        if gt_z(i,3) >= 0
                        
            % ----------------------------------------------------
            % linear color code
            % ----------------------------------------------------
            if strcmp(color_scale,'linear')
            % {
            linear_color_number = round( min(gt_z(i,3),max_concentration)/delta)+1;
            col = gdm_colormap( linear_color_number, :);
            %}
            end
            % ----------------------------------------------------

            % ----------------------------------------------------
            % nonlinear color code
            % ----------------------------------------------------
            if strcmp(color_scale,'nonlinear')
            % {
            linear_color_number = round( min(gt_z(i,3),max_concentration)/delta)+1;
            %nonlinear_color_number = round( a*log(b*linear_color_number) );
            %nonlinear_color_number = linear_color_number;
            nonlinear_color_number = normlogsum_penalty(linear_color_number);
            %if nonlinear_color_number <1
            %    nonlinear_color_number = 1;
            %end
            col = gdm_colormap( nonlinear_color_number, :);                
            %}
            end
            % ----------------------------------------------------
            
            %----- to avoid black borders of earch cell -----
            % {
            if sum(col) == 3
                col = col-(1e-10);
            end
            %}
            
            %plot(gt_z(i,1),gt_z(i,2),...
            %    's',...
            %    'Color',col,...
            %    'MarkerFaceColor',col,...
            %    'MarkerSize',4.5);
            plot(gt_z(i,1)+0.5,gt_z(i,2)+0.5,...
                's',...
                'MarkerEdgeColor',col,... %[0.2,0.2,0.2]
                'MarkerFaceColor',col,...
                'MarkerSize',cell_marker_size);
            
        end
    end
end


%_________________________________________________________________________________________
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             SENSING COVERAGE
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishSensingCoverage
    
    
    files = dir([dir_.logs,'sensor*'])    
    M = [];
    for i = 1:size(files,1)%6

        %sensor01_frame
        %file_name = sprintf('sensor%02d_frame.dat',i)
        file_name = files(i).name
        current_M = dlmread([dir_.logs,file_name]);
        %current_M(:,5) = current_M(:,5)/100000;
        %current_M(:,5) = current_M(:,5)/100000;

        M = [M;current_M];

        % publish (simulator results)
        M_publish_x = current_M(:,[1,3]);
        M_publish_y = current_M(:,[2,4]);

        % matlab resutls
        %M_publish_x = current_M(:,[3,5]);
        %M_publish_y = current_M(:,[4,6]);
        for j = 1:size(current_M,1)    
            plot(M_publish_x(j,:)+0.5,M_publish_y(j,:)+0.5,'-','Color',beam_color); 
        end
    end
    
    M_high = M(M(:,end)>=(para_.PublishHighConcentrationTH),:);
    for j = 1:size(M_high,1)
        plot(M_high(j,[1,3])+0.5,M_high(j,[2,4])+0.5,'-r');
    end
  

end


% --
%{
figure; hold on;
[env_row,env_col] = size(map_env)

for i = 1:env_row
    for j = 1:env_col

        if map_env(i,j)==0
            col = [0,0,0];
        elseif map_env(i,j)==1
            col = [1,1,1];
        end
        % ----------------------------------------------------
        plot(i,j,'s','Color',col,'MarkerFaceColor',col,'MarkerSize',5);
    end
end
axis equal
%}

% close all
% figure; imshow(map_env','InitialMagnification',800); hold on; % indoor/outdoor
% set(gca,'YDir','normal')


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             CANDIDATE CONFIGURATIONS MAP
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishCandConf
    
    map_candConf = dataPublish.map_candConf;
    
    for i = 1:map_row
        for j = 1:map_col
            % if obstacle/free            
            if map_candConf(i,j) % free                
                plot(i,j,'xb','MarkerSize',marker_size_cand_conf);
            end
        end
    end
   %axis equal
end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             HOTSPOTS/CENTER OF MASSES:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishHotspots
    
    hotspots = dataPublish.hotspots;
    
    [hot_r,hot_c] = ind2sub(size(map_env),hotspots);
    plot(hot_r,hot_c,'or','MarkerSize',marker_size_hotspot); %'MarkerFaceColor',[1,0,0]
    
end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             MARK SUBJECT HOTSPOTS:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishSubH
    
    subHotspots = dataPublish.subHotspots;    
    plot(subHotspots.sub(:,1),subHotspots.sub(:,2),...
        'og','MarkerSize',marker_size_hotspot);
end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             ROBOT START POSITION:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishStartPosition
    
    robotStartCellInMap = dataPublish.robotStartCellInMap;
    plot(robotStartCellInMap(1),robotStartCellInMap(2),...
        'ob','MarkerFaceColor',[0,0,1],'MarkerSize',marker_size_conf_position);
    
    text(robotStartCellInMap(1)+0.3,robotStartCellInMap(2),...
        '0','HorizontalAlignment','left','VerticalAlignment','middle',...
        'FontWeight','demi',...
        'FontSize',font_size_conf_num,...
        'Color',[0,0,0.75],...
        'Interpreter','latex',...
        'FontWeight','bold');
end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                        MARK WEAK/REDUNDANT CONFIGURATIONS:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishWeakConf
    
    weakConfGlobal_ind = dataPublish.weakConfGlobal_ind;

    
    % conf cell ind
    weakConfCell = fix(weakConfGlobal_ind/f)+(~~mod(weakConfGlobal_ind,f));   % cell num for each conf num.
    % conf num within a cell
    weakConfNum = mod(weakConfGlobal_ind,f)+(~mod(weakConfGlobal_ind,f)*f);  % conf num within a cell.
    % row and col of conf cell
    [weakConfCell_r,weakConfCell_c] = ind2sub(size(map_env),weakConfCell); % row and col of cell num
    
    % mark positions
    plot(weakConfCell_r,weakConfCell_c,...
        'or','MarkerFaceColor','r','MarkerSize',marker_size_weak_conf);
    % Draw arrow    
    for i = 1:numel(weakConfGlobal_ind)
        
%         h = quiver(weakConfCell_c(i),weakConfCell_r(i),...
%             2*cosd(360-FoVfaces.ang(weakConfNum(i),1)),...
%             2*sind(360-FoVfaces.ang(weakConfNum(i),1)),...
%             'LineWidth',1,'color','r');    
%         adjust_quiver_arrowhead_size(h,9.5);
%         %adjust_quiver_arrowhead_size(h,1);
        
        
        % Draw arrow
        quiver(weakConfCell_r(i),weakConfCell_c(i),...
           arrow_length*cosd(FoVfaces.lgt(weakConfNum(i))),...
           arrow_length*sind(FoVfaces.lgt(weakConfNum(i))),...
           'LineWidth',line_width_arrow,...
           'color','b',...
           'MaxHeadSize',arrow_head_size,...
           'Marker','o',...
           'MarkerSize',marker_size_arrow);

    end    
end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                       MARK CONDENSED/STRONG CONFIGURATIONS:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishCondConf
    
    condensedConf_ind = dataPublish.condensedConf_ind;

    
    % conf cell ind
    condConfCell = fix(condensedConf_ind/f)+(~~mod(condensedConf_ind,f));   % cell num for each conf num.
    % conf num within a cell
    condConfNum = mod(condensedConf_ind,f)+(~mod(condensedConf_ind,f)*f);  % conf num within a cell.
    % row and col of conf cell
    [condConfCell_r,condConfCell_c] = ind2sub(size(map_env),condConfCell); % row and col of cell num
    
    % -- mark positions   
    plot(condConfCell_r,condConfCell_c,...
        'o','color',[0.5,0.5,0.5],...
        'MarkerFaceColor',[0.5,0.5,0.5],...
        'MarkerSize',marker_size_cond_conf);
    
    % Draw arrow
    %{
    for i = 1:numel(condensedConf_ind)        
        h = quiver(condConfCell_c(i),condConfCell_r(i),...
            2*cosd(360-FoVfaces.ang(condConfNum(i),1)),2*sind(360-FoVfaces.ang(condConfNum(i),1)),...
            'LineWidth',1,'color',[0.5,0.5,0.5]);    
        adjust_quiver_arrowhead_size(h,9.5);
        %adjust_quiver_arrowhead_size(h,1);
    end
    %}
end



%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             MARK CONFIGURATIONS:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishConfPosition
    
    Conf = dataPublish.Conf;
    C = Conf;
    iC  = find(C);                      % IDs of selected conf.
    nCELL = iC; %fix(iC/f)+(~~mod(iC,f));   % cell num for each conf num.
    nCONF = 1*ones(1,numel(iC)) %mod(iC,f)+(~mod(iC,f)*f);  % conf num within a cell.
    [nCELLr,nCELLc] = ind2sub(size(map_env),nCELL); % row and col of cell num.
    
    plot(nCELLr,nCELLc,...
        'ob','MarkerFaceColor','b','MarkerSize',marker_size_conf_position);
end



%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             DRAW FIELD OF VIEW:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishFoV
    
    FoV = dataPublish.FoV;
    SensorRange = dataPublish.SensorRange;
    
    for i = 1:length(iC)
                
        % Draw FoV
        start_angle  = FoVfaces.ang(nCONF(i),1); % first/start angle.
        sector_angle = FoV;                      % increment in the first/start angle.
        MaxPts       = 1000;                     % maximum points to plot the FoV.
        [xx,yy]      = SecDraw(start_angle,sector_angle,SensorRange.cell,MaxPts); % x,y points.
        plot(xx+nCELLr(i),yy+nCELLc(i),'color','b','LineWidth',line_width_fov);
    end
end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             DRAW SYMBOLIC FIELD OF VIEW:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishSymbolicFoV
    
    FoV = dataPublish.FoV;
    
    %arrow_length = 2;
    %line_length  = 2;
    
    ConfOrientations = dataPublish.Orientations;
    
    for i = 1:length(iC)
        
        % -- drat sumbolic field of view
        
        
        start_angle  = ConfOrientations-FoV/2; %FoVfaces.ang(nCONF(i),1);
        end___angle  = ConfOrientations+FoV/2; %start_angle+FoV;
        r            = symbolic_fov_length;
        x1           = r*cosd(start_angle);
        y1           = r*sind(start_angle);
        x2           = r*cosd(end___angle);
        y2           = r*sind(end___angle);
        plot([nCELLr(i),x1+nCELLr(i)],...
             [nCELLc(i),y1+nCELLc(i)],...
              'color',[0.0 0.0 0.6],'LineWidth',line_width_symbolic_fov); 
        plot([nCELLr(i),x2+nCELLr(i)],...
             [nCELLc(i),y2+nCELLc(i)],...
              'color',[0.0 0.0 0.6],'LineWidth',line_width_symbolic_fov); 
    end
end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             DRAW CONF ARROWS:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishArrow
    
    ConfNumbers      = dataPublish.ConfNum;
    ConfOrientations = dataPublish.ConfOrn;
    
    [ConfCell_r,ConfCell_c] = ind2sub(size(map_env),ConfNumbers);
    
    for i = 1:length(ConfNumbers)
        
        conf_to_plot = [ConfNumbers(i),ConfCell_r(i),ConfCell_c(i),ConfOrientations(i)]
        
        %FoVfaces.lgt
        %nCONF(i)
        %FoVfaces.lgt(nCONF(i))
        
        % Draw arrow
        quiver(ConfCell_r(i),ConfCell_c(i),...
               arrow_length*cosd(ConfOrientations(i)),...
               arrow_length*sind(ConfOrientations(i)),...
               'LineWidth',line_width_arrow,...
               'color','b',...
               'MaxHeadSize',arrow_head_size,...
               'Marker','o',...
               'MarkerSize',marker_size_arrow);
        pause
    end
end

%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             PRINT CONF SEQUENCE NUMBERS:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishConfNum
    
    
    confSequenceNum = dataPublish.confSequenceNum;
    all_conf = zeros(length(iC),4); 
    
    for i = 1:length(iC)

        % conf number.       
        x3 = nCELLr(i)-(conf_num_dist*cosd(FoVfaces.lgt(nCONF(i))));
        y3 = nCELLc(i)-(conf_num_dist*sind(FoVfaces.lgt(nCONF(i))));
        %{
                            horizontal
                        (right,center,left)        
                <-------------------------------|
                                                |
                                                | vertical
                                                | (bottom,
                                                |  middle,
                                                |  top)
                                                µ 
        %}
        
        if     (sind(FoVfaces.lgt(nCONF(i)))==0 && cosd(FoVfaces.lgt(nCONF(i)))==1)
            HoriAlign = 'right';
            VertAlign = 'middle';
        elseif (sind(FoVfaces.lgt(nCONF(i)))==1 && cosd(FoVfaces.lgt(nCONF(i)))==0)
            HoriAlign = 'center';
            VertAlign = 'top'; 
        elseif (sind(FoVfaces.lgt(nCONF(i)))==0 && cosd(FoVfaces.lgt(nCONF(i)))==-1)
            HoriAlign = 'left';
            VertAlign = 'middle';
        elseif (sind(FoVfaces.lgt(nCONF(i)))==-1 && cosd(FoVfaces.lgt(nCONF(i)))==0)            
            HoriAlign = 'center';
            VertAlign = 'bottom';
        else
            error('WTF');
        end
        
        text(x3,y3,num2str(confSequenceNum(i)),...
            'HorizontalAlignment',HoriAlign,...
            'VerticalAlignment',VertAlign,...
            'FontWeight','demi',...
            'FontSize',font_size_conf_num,...
            'Color',[0,0,0.75],...
            'Interpreter','latex',...
            'EdgeColor',[0.25,0.25,1.0],...
            'Margin',1,...
            'FontWeight','bold');
        
                
        this_conf = [confSequenceNum(i),nCELLr(i), nCELLc(i), FoVfaces.lgt(nCONF(i))];
        all_conf(i,:) = this_conf;
    end
    
    % all_conf
    all_conf = sortrows(all_conf,1)

    % -- save
    % {
    file_name = 'robot_plan_detection.dat';
    %file_name = 'robot_plan_tomography.dat';
    fileID = fopen([dir_.results,file_name],'wt'); fclose(fileID); % writing to follow
    dlmwrite([dir_.results,file_name],all_conf,'-append','delimiter',' ');
    fileID = fopen([dir_.results,file_name],'a'); fclose(fileID);
    %}
    
end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                            IN-OUT TSP ANGLES FOR HOTSPOTS:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishHotspotsInOutTSPAngles
    
    
    
    InOutTSPAnglesMean  = dataPublish.InOutTSPAnglesMean;
    
    
    
    for i = 1:size(InOutTSPAnglesMean,1)
        
        %hotspotCell_ind = hotspotSEQ(i)
        hotspotCell_ind = hotspots(i);
        [hotspotCell_r,hotspotCell_c] = ind2sub(size(map_env),hotspotCell_ind);
        
        % entrance mark
        entranceMark_r = hotspotCell_r+(4*cosd(InOutTSPAnglesMean(i,1)));
        entranceMark_c = hotspotCell_c+(4*sind(InOutTSPAnglesMean(i,1)));
        
        % angle from enterance mark to hotspot
        angleMarkToHotspot = rad2deg(atan2(hotspotCell_c-entranceMark_c,hotspotCell_r-entranceMark_r));
        
        % exit mark
        exitMark_r = hotspotCell_r+(4*cosd(InOutTSPAnglesMean(i,2)));
        exitMark_c = hotspotCell_c+(4*sind(InOutTSPAnglesMean(i,2)));
        
       
        % plot entrance arrow
        enterH = quiver(entranceMark_c,entranceMark_r,...
            3.5*sind(angleMarkToHotspot),3.5*cosd(angleMarkToHotspot),...
            'LineWidth',1,'color','b');
        adjust_quiver_arrowhead_size(enterH,8.5);
        
        % plot exit arrow
        exitH = quiver(hotspotCell_c+(0.5*sind(InOutTSPAnglesMean(i,2))),hotspotCell_r+(0.5*cosd(InOutTSPAnglesMean(i,2))),...
            4*sind(InOutTSPAnglesMean(i,2)),4*cosd(InOutTSPAnglesMean(i,2)),...
            'LineWidth',1,'color','b');
        adjust_quiver_arrowhead_size(exitH,8.5);
        
        
    end
end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                       REFERENCE POINT FOR DISTANCE MATRIX:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% ---- mark reference point to calculate distance
if para_.PublishReferenceDistPoint    
    
    refDistPoint = dataPublish.refDistPoint;
    plot(refDistPoint(1),refDistPoint(2),'ob');
    
end


% hold off
% alpha(.1)

end

% ----------------------- End of the document ----------------------------------
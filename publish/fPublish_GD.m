function fPublish_GD( dataPublish,para_,dir_ )
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


% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                         PUBLISHING PARAMETERS
% 
% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% zoom_in_environment = 350; %800
% cell_marker_size = 2.75; %4.85 
% color_scale = 'linear'; %
% %color_scale = 'nonlinear'; % main fig.
% beam_color = [255,204,153]/256; % [204,229,255]/256
% marker_size_cand_conf = 2.5;
% marker_size_hotspot = 7;
% marker_size_weak_conf = 4;
% marker_size_cond_conf = 5;
% marker_size_conf_position = 2; %5; % icra-submission size: 3
% line_width_fov = 1; %2.5
% symbolic_fov_length = 3;
% line_width_symbolic_fov = 0.5; %1.25;
% 
% arrow_length = 3;
% line_width_arrow = 0.5; %0.25
% arrow_head_size = 20;
% marker_size_arrow = 1; %2
% 
% font_size_conf_num = 6;%16
% conf_num_dist = 0.5; % 0.3

if strcmp(para_.PublishTextWidth,'single-column') % -- single column.
    
    zoom_in_environment = 350; %750;
    cell_marker_size = 2.95; %10; %4; %10.25 %4;%4.85;
    color_scale = 'linear'; %
    %beam_color = [255,204,153]/256; % [204,229,255]/256
    marker_size_cand_conf = 2.5;
    marker_size_hotspot = 5; %7;
    %marker_size_weak_conf = 4;
    %marker_size_cond_conf = 5;
    marker_size_conf_position = 2; %3; %5; % icra-submission size: 3
    line_width_fov = 1; %2.5
    symbolic_fov_length = 2.0; %3;
    line_width_symbolic_fov = 1; %0.5; %1.25;

    arrow_length = 3;
    line_width_arrow = 1; %0.5; %0.25
    arrow_head_size = 30; %20;
    marker_size_arrow = 1; %2

    font_size_conf_num = 10; %6;%16
    font_size_conc_scale = 10;
    conf_num_dist = 1.5; %0.5; % 0.3
    
elseif strcmp(para_.PublishTextWidth,'two-columns') % -- two columns.
    
    zoom_in_environment = 350; %800
    cell_marker_size = 2.75; %4.85 
    color_scale = 'linear'; %
    %color_scale = 'nonlinear'; % main fig.
    %beam_color = [255,204,153]/256; % [204,229,255]/256
    marker_size_cand_conf = 2.5;
    marker_size_hotspot = 3; %7;
    %marker_size_weak_conf = 4;
    %marker_size_cond_conf = 5;
    marker_size_conf_position = 2; %5; % icra-submission size: 3
    line_width_fov = 1; %2.5
    symbolic_fov_length = 2;
    line_width_symbolic_fov = 0.5; %1.25;

    arrow_length = 3;
    line_width_arrow = 0.5; %0.25
    arrow_head_size = 20;
    marker_size_arrow = 1; %2

    font_size_conf_num = 5;%16
    font_size_conc_scale = 8;
    conf_num_dist = 0.5; % 0.3

end



% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             ENVIRONMENT MAP:
% 
% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% map_env = dlmread([dir_.env,para_.EnvironmentFigFile]);
% map_env = 1-map_env;
map_env = dataPublish.map_env;
[map_row,map_col] = size(map_env);

h = figure('name',para_.PublishFigureTitle);
imshow(map_env','InitialMagnification',zoom_in_environment); 
hold on;
set(gca,'YDir','normal')
colormap('gray'); % not go back to gray scale.
%caxis([-5 2]) % for resulting gdm
%caxis([-12 4]) % for sensor placement solution
% caxis([-0.25 1]) % for sensor placement solution
% caxis([-1 1.25])
caxis([-0.5 1.0])
% caxis([-2 2])
%caxis([-10 1.5]) % for gdm insets

%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             GAS CONCENTRATION MAP
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishGDM
    
    % -- reconstruction        
    load([dir_.gdplan_recon,para_.ReconstructionFile],...
        'mean_map',...
        'cell_coords_x',...
        'cell_coords_y');
    map_gd = mean_map';
    
    
    % ------------- select color map for gdm --------------
    gdm_colormap = flipud(colormap('hot'));
    % -----------------------------------------------------
    
    
    % ------ color step size --------------------------------    
    max_concentration = max(map_gd(:));
    %max_concentration = 500;
    delta = max_concentration/(size(gdm_colormap,1)-1);
    
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
        
    for i = 1:numel(cell_coords_x)%1:numel(cell_coords_y)%map_row
        for j = 1:numel(cell_coords_y)%1:numel(cell_coords_x)%map_col
            % if positive concentration
            
            if map_gd(i,j) >= 0 %>=0
                
                if map_env(i,j) > 0
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
                plot(cell_coords_x(i)+0.0,cell_coords_y(j)+0.0,...
                    's',...
                    'MarkerEdgeColor',col,...
                    'MarkerFaceColor',col,...
                    'MarkerSize',cell_marker_size);
                end
                
            end
        end
    end
    
    %{
    [obs_r,obs_c] = find(map_env==0);
    plot(obs_r,obs_c,...
            's',...
            'MarkerEdgeColor',[0.5,0.5,0.5],...
            'MarkerFaceColor',[0.5,0.5,0.5],...
            'MarkerSize',cell_marker_size+0.1);
    pause
    %}
    
    % --- concentration scale.
    %{
    conc_text = sprintf('Concentration scale: %d -- %d ppm',0,max_concentration);
    text(map_row-2,2,conc_text,...
            'HorizontalAlignment','right',...
            'VerticalAlignment','bottom',...
            'FontWeight','demi',...
            'FontSize',font_size_conc_scale,...
            'Color',[0,0,0],...
            'BackgroundColor',[1,1,1],...
            'Interpreter','latex',...
            'EdgeColor',[0,0,0],...
            'Margin',1,...
            'FontWeight','bold');
    %}
end



%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             GROUND TRUTH
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishGroundTruth
    
    
    % -- ground truth
    gt = dlmread([dir_.gt,'map_con.dat']);
    
    % ------------- select color map for gdm --------------
    gdm_colormap = flipud(hot(512));
    % -----------------------------------------------------
    
    
    % ------ color step size --------------------------------    
    max_concentration = max(gt(:));
    delta = max_concentration/(size(gdm_colormap,1)-1);
    
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
    
    [gt_row,gt_col] = size(gt);
    
    for i = 1:gt_row
        for j = 1:gt_col
            % if positive concentration
            
            if gt(i,j)>=0
                
                if map_env(i,j)>0
                
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
                linear_color_number = round( min(gt(i,j),max_concentration)/delta)+1;
                col = gdm_colormap( linear_color_number,:);                
                %}
                end
                % ----------------------------------------------------
                
                % ----------------------------------------------------
                % nonlinear color code
                % ----------------------------------------------------
                if strcmp(color_scale,'nonlinear')
                % {
                linear_color_number = round( min(gt(i,j),max_concentration)/delta)+1;
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
                plot(i+0.0,j+0.0,...
                    's',...
                    'MarkerEdgeColor',col,...
                    'MarkerFaceColor',col,...
                    'MarkerSize',cell_marker_size);
                end
            end
        end
    end
    cell_marker_size
    % --- concentration scale.
    %{
    conc_text = sprintf('Concentration scale: %d -- %.2d ppm',0,max_concentration);
    text(map_row-2,2,conc_text,...
            'HorizontalAlignment','right',...
            'VerticalAlignment','bottom',...
            'FontWeight','demi',...
            'FontSize',font_size_conc_scale,...
            'Color',[0,0,0],...
            'BackgroundColor',[1,1,1],...
            'Interpreter','latex',...
            'EdgeColor',[0,0,0],...
            'Margin',1,...
            'FontWeight','bold');
    %}
end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             SENSING COVERAGE
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishSensingCoverage
    
    
    files = dir([dir_.logs,para_.MeasurementFiles]);
    M = [];
    for i = 1:size(files,1)%6

        %sensor01_frame
        %file_name = sprintf('sensor%02d_frame.dat',i)
        file_name = files(i).name;
        current_M = dlmread([dir_.logs,file_name]);
        %current_M(:,5) = current_M(:,5)/100000;
        %current_M(:,5) = current_M(:,5)/100000;
        
        %max(current_M(:,8))
        
        M = [M;current_M];
        
        %{
        % publish (simulator results)
        M_publish_x = current_M(:,[1,3]);
        M_publish_y = current_M(:,[2,4]);

        % matlab resutls
        %M_publish_x = current_M(:,[3,5]);
        %M_publish_y = current_M(:,[4,6]);
        for j = 1:size(current_M,1)    
            plot(M_publish_x(j,:)+0.5,M_publish_y(j,:)+0.5,'-','Color',beam_color); 
        end
        %}
    end
    
    %{
    M_high = M(M(:,end)>=(para_.PublishHighConcentrationTH),:);
    for j = 1:size(M_high,1)
        plot(M_high(j,[1,3])+0.5,M_high(j,[2,4])+0.5,'-r');
    end
    %}
    
    % {    
    start_x__all = M(:,3);
    start_y__all = M(:,4);
    end___x__all = M(:,5);
    end___y__all = M(:,6);
    weight___all = M(:,8)/max(M(:,8));
    
    %max_concentration_optical_beams= max(M(:,8));

    beam_colormap = flipud(bone(512));
    for i = 1:size(start_x__all,1)
        plot([start_x__all(i),end___x__all(i)]-0.5,...
             [start_y__all(i),end___y__all(i)]-0.5,...
             '-','Color',beam_colormap(round(255*weight___all(i))+100,:));
    end
    %}
    %{
    % --- concentration scale.
    conc_text = sprintf('Concentration scale: %d -- %.2d ppm',0,max(M(:,8)));
    text(map_row-2,2,conc_text,...
            'HorizontalAlignment','right',...
            'VerticalAlignment','bottom',...
            'FontWeight','demi',...
            'FontSize',font_size_conc_scale,...
            'Color',[0,0,0],...
            'BackgroundColor',[1,1,1],...
            'Interpreter','latex',...
            'EdgeColor',[0,0,0],...
            'Margin',1,...
            'FontWeight','bold');
    %}
end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             DETECTION EVENTS
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishDetectionEvents
    
    
    files = dir([dir_.logs,para_.MeasurementFiles]);
    M = [];
    
    conf_detection_flag = zeros(size(files,1),1);
    
    for i = 1:size(files,1)%6

        %sensor01_frame
        %file_name = sprintf('sensor%02d_frame.dat',i)
        file_name = files(i).name;
        current_M = dlmread([dir_.logs,file_name]);
        %current_M(:,5) = current_M(:,5)/100000;
        %current_M(:,5) = current_M(:,5)/100000;
        
        if max(current_M(:,8)) >= 1000
            conf_detection_flag(i) = 1;
        end
        
        M = [M;current_M];        
    end
    
    % {    
    start_x__all = M(:,3);
    start_y__all = M(:,4);
    end___x__all = M(:,5);
    end___y__all = M(:,6);
    ppm______all = M(:,8);
    
    %max_concentration_optical_beams= max(M(:,8))
    %beam_colormap = flipud(bone(512));
    
    for i = 1:size(start_x__all,1)
        
        if ppm______all(i) >= 500
        
        %plot([start_x__all(i),end___x__all(i)]-0.5,...
        %     [start_y__all(i),end___y__all(i)]-0.5,...
        %     '-','Color',[255,69,0]/255);
         
        end
    end
    %}
    %{
    % --- concentration scale.
    conc_text = sprintf('Concentration scale: %d -- %.2d ppm',0,max(M(:,8)));
    text(map_row-2,2,conc_text,...
            'HorizontalAlignment','right',...
            'VerticalAlignment','bottom',...
            'FontWeight','demi',...
            'FontSize',font_size_conc_scale,...
            'Color',[0,0,0],...
            'BackgroundColor',[1,1,1],...
            'Interpreter','latex',...
            'EdgeColor',[0,0,0],...
            'Margin',1,...
            'FontWeight','bold');
    %}
end



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
%                             MARK CONFIGURATIONS:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishConfPosition
    
    % ----------------------------------
    FoVfaces = dataPublish.FoVfaces;
    f = FoVfaces.num;
    % ----------------------------------

    Conf  = dataPublish.Conf;
    C     = Conf;
    iC    = find(C);                   % IDs of selected conf.
    nCELL = fix(iC/f)+(~~mod(iC,f));   % cell num for each conf num.
    nCONF = mod(iC,f)+(~mod(iC,f)*f);  % conf num within a cell.
    [nCELLr,nCELLc] = ind2sub(size(map_env),nCELL); % row and col of cell num.
    
    plot(nCELLr,nCELLc,...
        'ob','MarkerFaceColor','b','MarkerSize',marker_size_conf_position);
end



%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             DRAW FIELD OF VIEW:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% if para_.PublishFoV
%     
%     FoV = dataPublish.FoV;
%     SensorRange = dataPublish.SensorRange;
%     
%     for i = 1:length(iC)
%                 
%         % Draw FoV
%         start_angle  = FoVfaces.ang(nCONF(i),1); % first/start angle.
%         sector_angle = FoV;                      % increment in the first/start angle.
%         MaxPts       = 1000;                     % maximum points to plot the FoV.
%         [xx,yy]      = SecDraw(start_angle,sector_angle,SensorRange.cell,MaxPts); % x,y points.
%         %[xx,yy]      = SecDraw(start_angle,sector_angle,SensorRange.cell); % x,y points.
%         plot(xx+nCELLr(i),yy+nCELLc(i),'color','b','LineWidth',line_width_fov);
%     end
% end

if para_.PublishFoV
    
    FoV = dataPublish.FoV;
    SensorRange = dataPublish.SensorRange;
    
    
    %------------------------------------------------
    % AREA INSIDE FIELD OF VIEW WITHOUT OCCLUSION
    %------------------------------------------------
    %{
    for i = 1:length(iC)                
        % Draw FoV
        start_angle  = FoVfaces.ang(nCONF(i),1); % first/start angle.
        sector_angle = FoV;                      % increment in the first/start angle.
        MaxPts       = 1000;                     % maximum points to plot the FoV.
        [xx,yy]      = SecDraw(start_angle,sector_angle,SensorRange.cell,MaxPts); % x,y points.
        %hfill = fill(xx+nCELLr(i),yy+nCELLc(i),[.4,.4,1]);
        hfill = fill(xx+nCELLr(i),yy+nCELLc(i),[1,.4,0.4]);
        set(hfill,'facealpha',0.25)
        pause
    end
    %}
    
    %------------------------------------------------
    % AREA INSIDE FIELD OF VIEW WITH OCCLUSION
    % (TEMPORARY SOLUTION BASED ON COLLECTED MEASUREMENTS)
    %------------------------------------------------
    %files = dir([dir_.logs,para_.MeasurementFiles]);
    
    meas_files = 'measurement_conf*.dat*';
    files = dir([dir_.logs,meas_files]);    
    
    %fill_color = [1.0,0.4,0.4];
    fill_color = [0.4,0.4,1.0];
    
    for i = 1:size(files,1)

        file_name = files(i).name;
        current_M = dlmread([dir_.logs,file_name]);
        
        start_x__all = current_M(:,3);
        start_y__all = current_M(:,4);
        end___x__all = current_M(:,5);
        end___y__all = current_M(:,6);
        
        hfill = fill([start_x__all(1);end___x__all]-0.5,...
                     [start_y__all(1);end___y__all]-0.5,...
                     fill_color,'EdgeColor','none');
        set(hfill,'facealpha',0.25)
        %pause
    end
    
    
    
    %------------------------------------------------
    % CIRCULAR SECTOR
    %------------------------------------------------
    for i = 1:length(iC)
                
        % Draw FoV
        start_angle  = FoVfaces.ang(nCONF(i),1); % first/start angle.
        sector_angle = FoV;                      % increment in the first/start angle.
        MaxPts       = 1000;                     % maximum points to plot the FoV.
        [xx,yy]      = SecDraw(start_angle,sector_angle,SensorRange.cell,MaxPts); % x,y points.
        %[xx,yy]      = SecDraw(start_angle,sector_angle,SensorRange.cell); % x,y points.
        %plot(xx+nCELLr(i),yy+nCELLc(i),'color','b','LineWidth',line_width_fov);
        plot(xx+nCELLr(i),yy+nCELLc(i),'color',[0,0,.6],'LineWidth',line_width_fov);
        %hfill = fill(xx+nCELLr(i),yy+nCELLc(i),'r')
        %set(hfill,'facealpha',0.5)
        %pause
    end
end

%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             DRAW SYMBOLIC FIELD OF VIEW:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishSymbolicFoV
    
    FoV = dataPublish.FoV;
    
    
    for i = 1:length(iC)
        
        % -- drat sumbolic field of view
        start_angle  = FoVfaces.ang(nCONF(i),1);
        end___angle  = start_angle+FoV;
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
    
    for i = 1:length(iC)
        
        % Draw arrow
        quiver(nCELLr(i),nCELLc(i),...
               arrow_length*cosd(FoVfaces.lgt(nCONF(i))),...
               arrow_length*sind(FoVfaces.lgt(nCONF(i))),...
               'LineWidth',line_width_arrow,...
               'color','b',...
               'MaxHeadSize',arrow_head_size,...
               'Marker','o',...
               'MarkerSize',marker_size_arrow);
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
        
        %{
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
        %}
        
        sin_test = round(sind(FoVfaces.lgt(nCONF(i)))*1e4)/1e4;
        cos_test = round(cosd(FoVfaces.lgt(nCONF(i)))*1e4)/1e4;
        
        if     sin_test ==  0.0000 && cos_test ==  1.0000 % 000 deg
            HoriAlign = 'right';
            VertAlign = 'middle';
        elseif sin_test ==  0.7071 && cos_test ==  0.7071 % 045 deg
            HoriAlign = 'right';
            VertAlign = 'top'; 
        elseif sin_test ==  1.0000 && cos_test ==  0.0000 % 090 deg
            HoriAlign = 'center';
            VertAlign = 'top';
        elseif sin_test ==  0.7071 && cos_test == -0.7071 % 135 deg
            HoriAlign = 'left';
            VertAlign = 'top'; 
        elseif sin_test ==  0.0000 && cos_test == -1.0000 % 180 deg
            HoriAlign = 'left';
            VertAlign = 'middle';
        elseif sin_test == -0.7071 && cos_test == -0.7071 % 225 deg
            HoriAlign = 'left';
            VertAlign = 'bottom'; 
        elseif sin_test == -1.0000 && cos_test ==  0.0000 % 270 deg
            HoriAlign = 'center';
            VertAlign = 'bottom';
        elseif sin_test == -0.7071 && cos_test ==  0.7071 % 315 deg
            HoriAlign = 'right';
            VertAlign = 'middle';
        else
            error('WTF');
        end
        
        
        
        text(x3,y3,num2str(confSequenceNum(i)),...
            'HorizontalAlignment',HoriAlign,...
            'VerticalAlignment',VertAlign,...
            'FontWeight','demi',...
            'FontSize',font_size_conf_num,...
            'Color',[0,0,0.75],...
            'Interpreter','latex',...           %'EdgeColor',[0.25,0.25,1.0],...
            'Margin',1,...
            'FontWeight','bold');
        
                
        this_conf = [confSequenceNum(i),nCELLr(i), nCELLc(i), FoVfaces.lgt(nCONF(i))];
        all_conf(i,:) = this_conf;
        
        
        
    end
    
    % all_conf
    all_conf = sortrows(all_conf,1);

    % -- save
    %{
    file_name = 'robot_plan_detection.dat';
    %file_name = 'robot_plan_tomography.dat';
    fileID = fopen([dir_.results,file_name],'wt'); fclose(fileID); % writing to follow
    dlmwrite([dir_.results,file_name],all_conf,'-append','delimiter',' ');
    fileID = fopen([dir_.results,file_name],'a'); fclose(fileID);
    %}
    
end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             DETECTION FIELD OF VIEW:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishDetectionEvents
    
    
    % ----------------------------------
    FoVfaces = dataPublish.FoVfaces;
    f = FoVfaces.num;
    % ----------------------------------
    
    FoV = dataPublish.FoV;
    SensorRange = dataPublish.SensorRange;
    
    confSequenceNum = dataPublish.confSequenceNum;
    
    
    Conf  = dataPublish.Conf;
    C     = Conf;
    iC    = find(C);                   % IDs of selected conf.
    nCELL = fix(iC/f)+(~~mod(iC,f));   % cell num for each conf num.
    nCONF = mod(iC,f)+(~mod(iC,f)*f);  % conf num within a cell.
    [nCELLr,nCELLc] = ind2sub(size(map_env),nCELL); % row and col of cell num.
    
    for i = 1:length(iC)
        
        if conf_detection_flag(confSequenceNum(i)) > 0
                
        % Draw FoV
        start_angle  = FoVfaces.ang(nCONF(i),1); % first/start angle.
        sector_angle = FoV;                      % increment in the first/start angle.
        MaxPts       = 1000;                     % maximum points to plot the FoV.
        [xx,yy]      = SecDraw(start_angle,sector_angle,SensorRange.cell,MaxPts); % x,y points.
        %[xx,yy]      = SecDraw(start_angle,sector_angle,SensorRange.cell); % x,y points.
        plot(xx+nCELLr(i),yy+nCELLc(i),'color','r','LineWidth',line_width_fov);
        
        plot(nCELLr(i),nCELLc(i),...
            'or','MarkerFaceColor','r','MarkerSize',marker_size_conf_position);
    
        end
    end
end



end

% ----------------------- End of the document ----------------------------------

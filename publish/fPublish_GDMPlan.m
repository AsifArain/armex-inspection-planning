function fPublish_GDMPlan( dataPublish,para_,dir_ )
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
% FoVfaces = dataPublish.FoVfaces;
% f = FoVfaces.num;
% ----------------------------------


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                         PUBLISHING PARAMETERS
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤



if strcmp(para_.PublishTextWidth,'single-column') % -- single column.
    
       
    %{
    zoom_in_environment = 750;
    cell_marker_size = 10 %4.85;
    color_scale = 'linear'; %
    %beam_color = [255,204,153]/256; % [204,229,255]/256
    marker_size_cand_conf = 2.5;
    marker_size_hotspot = 10 %5; %7;
    line_width_hotspot = 2;
    %marker_size_weak_conf = 4;
    %marker_size_cond_conf = 5;
    marker_size_conf_position = 5 %2; %5; % icra-submission size: 3
    line_width_fov = 1; %2.5
    symbolic_fov_length = 2.0; %3;
    line_width_symbolic_fov = 1.5 % 0.5; %1.25;

    arrow_length = 3;
    line_width_arrow = 1.5 %0.5; %0.25
    arrow_head_size = 40 %30 %20;
    marker_size_arrow = 10 %1; %2

    font_size_conf_num = 15 %10%6;%16
    font_size_conc_scale = 10;
    %conf_num_dist = 0.5; % 0.3
    %conf_num_dist_adaptive_planned = 0.75;
    %}
    
    zoom_in_environment = 750;
    cell_marker_size = 10; %4.85;
    color_scale = 'linear'; %
    %beam_color = [255,204,153]/256; % [204,229,255]/256
    marker_size_cand_conf = 2.5;
    marker_size_hotspot = 5; %5; %7;
    line_width_hotspot = 1;
    %marker_size_weak_conf = 4;
    %marker_size_cond_conf = 5;
    marker_size_conf_position = 2; %5 %2; %5; % icra-submission size: 3
    line_width_fov = 1; %2.5
    symbolic_fov_length = 2.0; %3;
    line_width_symbolic_fov = 1.5; % 0.5; %1.25;

    arrow_length = 3;
    line_width_arrow = 1.5; %0.5; %0.25
    arrow_head_size = 40; %30 %20;
    marker_size_arrow = 10; %1; %2

    font_size_conf_num = 18; %15 %10%6;%16
    font_size_conc_scale = 10;
    %conf_num_dist = 0.5; % 0.3
    %conf_num_dist_adaptive_planned = 0.75;
    
    
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
    %conf_num_dist = 0.5; % 0.3
    %conf_num_dist_adaptive_planned = 0.75;
    
end

color_conf                  = [0.0 0.0 0.6];
color_conf_redundatn        = [135,206,235]/255; %[0.7,0.7,0.7];
color_conf_adaptive_planned = [0.4,0.4,0.4];

distance_conf_seq                  = 0.5; % 0.3
distance_conf_seq_adaptive_planned = 0.6;


% initialize legends
%----------------------------
h_hot = [];
h_pos = [];
h_travel = [];
h_beam = [];
h_sfov = [];
h_wsfov = [];

%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             ENVIRONMENT MAP:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% map_env = dlmread([dir_.env,para_.EnvironmentFigFile]);
% map_env = 1-map_env;
map_env = dataPublish.map_env;
[map_row,map_col] = size(map_env);

h = figure('name',para_.PublishFigureTitle); 
imshow(map_env','InitialMagnification',zoom_in_environment);
hold on;
set(gca,'YDir','normal');
colormap('gray'); 

%caxis([-5 2]) % for resulting gdm
%caxis([-12 4]) % for sensor placement solution
% caxis([-0.25 1]) % for sensor placement solution
% caxis([-2 2])
% caxis([-1 1.25])
caxis([-0.5 1.0])
%caxis([-10 1.5]) % for gdm insets


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             GAS CONCENTRATION MAP
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishGDM
    
    %-- reconstruction 
    load([dir_.recon,para_.ReconstructionFile],...
        'mean_map',...
        'cell_coords_x',...
        'cell_coords_y');
    map_recon = mean_map';
        
    %map_recon = load([dir_.recon,'reconstruction_conf01.dat']);
    %cell_coords_x = load([dir_.recon,'x_coord_conf01.dat']);
    %cell_coords_y = load([dir_.recon,'y_coord_conf01.dat']);
    
%     map_recon = load([dir_.recon,sprintf('reconstruction_conf%02d.dat',dataPublish.publish_gdm_for_executed_conf_num)]);
%     cell_coords_x = load([dir_.recon,sprintf('x_coord_conf%02d.dat',dataPublish.publish_gdm_for_executed_conf_num)]);
%     cell_coords_y = load([dir_.recon,sprintf('y_coord_conf%02d.dat',dataPublish.publish_gdm_for_executed_conf_num)]);
    
    %dlmwrite('recon_publish.dat',map_recon,'delimiter',' ');
    
    %--------------- select color map for gdm --------------
    %gdm_colormap = flipud(colormap('hot')); % for main fig. 
    gdm_colormap = flipud(hot(512)); % for main fig. 
    %gdm_colormap = flipud(autumn(512)); % for main fig. 
    %gdm_colormap = flipud(summer(512)); % for main fig. 
    %gdm_colormap = flipud(colormap('summer')); % for insets
    %-------------------------------------------------------
        
    %------ color step size --------------------------------
    %max_concentration = 10000; % inset gdm
    %gt = dlmread([dir_.GroundTruthLogs,'map_con.dat']);
    %max_concentration = max(gt(:))
    %max_concentration = max(map_recon(:)); % main fig.
    %max_concentration = 500;
    max(map_recon(:))
    para_.PublishUpperPPM
    %pause
    max_concentration = min( [max(map_recon(:)),para_.PublishUpperPPM] ) % main fig.
    %max_concentration = max(map_recon(:))
    %max_concentration = 5000
    delta = max_concentration/(size(gdm_colormap,1)-1);    
    %delta = max(map_gd(:))/(size(gdm_colormap,1)-1);
    %delta = 500/(size(gdm_colormap,1)-1); % lets limit to 500ppm
    
    %------------------- nonlinear color code ---------------------  
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
    %--------------------------------------------------------------
    
    xyzPoints = zeros(numel(cell_coords_x)*numel(cell_coords_y),3);
    colorPoints = zeros(numel(cell_coords_x)*numel(cell_coords_y),3);
    
    
    reconColorR = zeros(numel(cell_coords_x)*numel(cell_coords_y),1);
    reconColorG = zeros(numel(cell_coords_x)*numel(cell_coords_y),1);
    reconColorB = zeros(numel(cell_coords_x)*numel(cell_coords_y),1);
    
    k = 1;
    for i = 1:numel(cell_coords_x)%1:numel(cell_coords_y)%map_row
        for j = 1:numel(cell_coords_y)%1:numel(cell_coords_x)%map_col
            % if positive concentration
            
            if map_recon(i,j)>=0
                
                if map_env(i,j)>0
                %------------------- useful trash -------------------
                %col = gdm_colormap(round(map_gd(i,j)/delta)+1,:);
                
                % --- if concentration is greater than 1000 ppm, its still 1000 ppm
                %col = gdm_colormap( round(min(map_gd(i,j),3000)/delta)+1, :);
                %col = gdm_colormap( round(min(map_gd(i,j),1.78e+04)/delta)+1, :);
                %col = gdm_colormap( round(min(map_gd(i,j),inf)/delta)+1, :);
                %-----------------------------------------------------
                
                %----------------------------------------------------
                % linear color code
                %----------------------------------------------------
                if strcmp(color_scale,'linear')
                % {
                linear_color_number = round( min(map_recon(i,j),max_concentration)/delta)+1;
                col = gdm_colormap( linear_color_number,:);                
                %}
                end
                %----------------------------------------------------
                
                %----------------------------------------------------
                % nonlinear color code
                %----------------------------------------------------
                if strcmp(color_scale,'nonlinear')
                % {
                linear_color_number = round( min(map_recon(i,j),max_concentration)/delta)+1;
                %nonlinear_color_number = round( a*log(b*linear_color_number) );
                %nonlinear_color_number = linear_color_number;
                nonlinear_color_number = normlogsum_penalty(linear_color_number);
                %if nonlinear_color_number <1
                %    nonlinear_color_number = 1;
                %end
                col = gdm_colormap( nonlinear_color_number,:);                
                %}
                end
                %----------------------------------------------------
                
                %----- to avoid black borders of earch cell -----
                % {
                if sum(col) == 3
                    col = col-(1e-10);
                end
                %}
                %----- plot ------------
                xyzPoints(k,:) = [cell_coords_x(i),cell_coords_y(j),0];                
                colorPoints(k,:) = col;                
                reconColorR(k) = round(col(1)*255);
                reconColorG(k) = round(col(2)*255);
                reconColorB(k) = round(col(3)*255);
                
                plot(cell_coords_x(i)+0.0,cell_coords_y(j)+0.0,...
                    's',...
                    'MarkerEdgeColor',col,...
                    'MarkerFaceColor',col,...
                    'MarkerSize',cell_marker_size);
                k = k+1;
                
                end
            end
        end
    end
    
    %{
    map_recon = map_recon';
    linear_color_number = round(map_recon(:)./delta)+1;
    
    colR = round(gdm_colormap(linear_color_number,1)*255);
    colG = round(gdm_colormap(linear_color_number,2)*255);
    colB = round(gdm_colormap(linear_color_number,3)*255);
    
    thisfile = sprintf('reconstructionColorR.dat');
    fileID = fopen(thisfile,'wt'); fclose(fileID);
    dlmwrite(thisfile,colR,'-append','delimiter',' ');
    fileID = fopen(thisfile,'a'); fclose(fileID);
    
    thisfile = sprintf('reconstructionColorG.dat');
    fileID = fopen(thisfile,'wt'); fclose(fileID);
    dlmwrite(thisfile,colG,'-append','delimiter',' ');
    fileID = fopen(thisfile,'a'); fclose(fileID);
    
    thisfile = sprintf('reconstructionColorB.dat');
    fileID = fopen(thisfile,'wt'); fclose(fileID);
    dlmwrite(thisfile,colB,'-append','delimiter',' ');
    fileID = fopen(thisfile,'a'); fclose(fileID);
    
    %colR(find(colR(:)~=255))
    
    %pause
    %ptCloud = pointCloud(xyzPoints,'Color',colorPoints)
    %pcwrite(ptCloud,'testPlyFile','PLYFormat','binary')
    
    
    test_recon = map_recon;%+rand(size(map_recon));
    max(test_recon(:))
    %test_recon = test_recon./(max(test_recon(:)));
    
    %size(test_recon)
    %pause
    
    thisfile = sprintf('reconstructionColorR.dat');
    fileID = fopen(thisfile,'wt'); fclose(fileID);
    dlmwrite(thisfile,reconColorR,'-append','delimiter',' ');
    fileID = fopen(thisfile,'a'); fclose(fileID);
    
    thisfile = sprintf('reconstructionColorG.dat');
    fileID = fopen(thisfile,'wt'); fclose(fileID);
    dlmwrite(thisfile,reconColorG,'-append','delimiter',' ');
    fileID = fopen(thisfile,'a'); fclose(fileID);
    
    thisfile = sprintf('reconstructionColorB.dat');
    fileID = fopen(thisfile,'wt'); fclose(fileID);
    dlmwrite(thisfile,reconColorB,'-append','delimiter',' ');
    fileID = fopen(thisfile,'a'); fclose(fileID);
    
    thisfile = sprintf('reconstruction.dat');
    fileID = fopen(thisfile,'wt'); fclose(fileID);
    dlmwrite(thisfile,test_recon,'-append','delimiter',' ');
    fileID = fopen(thisfile,'a'); fclose(fileID);
    
    thisfile = sprintf('x_coord.dat');
    fileID = fopen(thisfile,'wt'); fclose(fileID);
    dlmwrite(thisfile,cell_coords_x,'-append','delimiter',' ');
    fileID = fopen(thisfile,'a'); fclose(fileID);
    
    thisfile = sprintf('y_coord.dat');
    fileID = fopen(thisfile,'wt'); fclose(fileID);
    dlmwrite(thisfile,cell_coords_y,'-append','delimiter',' ');
    fileID = fopen(thisfile,'a'); fclose(fileID);
    %}
    %pause
    
    % --- concentration scale.
    conc_text = sprintf('Concentration scale: %d -- %.2d ppm',0,max_concentration);
    %{
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


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             GROUND TRUTH
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishGroundTruth
    
    
    % -- ground truth
    gt = dlmread([dir_.GroundTruthLogs,'map_con.dat']);
    
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
    
    % --- concentration scale.
    conc_text = sprintf('Concentration scale: %d -- %.2d ppm',0,max_concentration);
    %{
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





%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             SENSING COVERAGE
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishSensingCoverage
    
    
    file__cellsize_env = sprintf('%s_cellsize_coverage.dat',para_.environment);
    cellsize_env = load([dir_.env,file__cellsize_env]);
    
    file__origin_env = sprintf('%s_origin_coverage.dat',para_.environment);
    origin_env = load([dir_.env,file__origin_env]);

    cellsize_recon = para_.ReconTomographyCellSizeM;

    
    %---------------------------------------------------------------------------
    switch para_.SensingSystem
        
            %=======================================
        case 'robot'
            %=======================================


            %--- robot logs to M matrix for all conf.
            %====================================================
            fprintf('---> Generating measurement matrix (M) for all the conf.\n');

            meas_files = 'measurements_conf*.dat*';
            files_names = dir([dir_.logs,meas_files]);

            visualize = 0;
            %{            
            for i = 2 %1:size(files_names,1)
                measure_file = sprintf('measurements_conf%d.dat',i);
                % visualize = 1;
                [ M_m,M_cell ] = fRobotLogs2M( measure_file,...
                                               para_,...
                                               dir_,...
                                               visualize );

                Mfile = sprintf('M_conf%02d.mat',i);
                save([dir_.MeasurementLogs,Mfile],'M_m','M_cell');
            end
            pause
            %}


            %--- M matrix for all conf.
            %====================================================

            % -- measurement matrix
            fprintf('---> Combining all measurements (M).\n');

            M = [];
            for num = 1:size(files_names,1)

                Mfile = sprintf('M_conf%02d.mat',num);
                M_current = load([dir_.MeasurementLogs,Mfile]);

                M = [M; M_current.M_m];
                %M = [M; M_current.M_cell];
            end
            
            M(:,1:4)   = M(:,1:4)./cellsize_env;        
            M(:,[1,3]) = M(:,[1,3])+origin_env(1);
            M(:,[2,4]) = M(:,[2,4])+origin_env(2);
            
            
            if visualize == 1

                %--- environment map
                figure('name','Optical beams'); 
                imshow(map_env','InitialMagnification',800); hold on;
                set(gca,'YDir','normal'); colormap('gray'); caxis([-2 0.5])
                plot(origin_env(1),origin_env(2),'ob');

                %{
                %--- optical beams (with color proportional to concentrations)
                beam_colormap = flipud(bone(512));
                beam_weight = M(:,5)/max(M(:,5));
                for i = 1:size(M,1)        
                    plot([M(i,1),M(i,3)],[M(i,2),M(i,4)],...
                        '-','Color',beam_colormap(round(255*beam_weight(i))+100,:));
                end
                %--- coverage area
                x_pts = [M(1,1);M(:,3);M(1,1)];
                y_pts = [M(1,2);M(:,4);M(1,2)];
                plot(x_pts-0.0,y_pts-0.0,'-','LineWidth',1.5,'color','b');
                %}


                M_env = M(:,1:4).*(cellsize_recon/cellsize_env);

                %--- optical beams (with color proportional to concentrations)
                beam_colormap = flipud(bone(512));
                beam_weight = M(:,5)/max(M(:,5));
                for i = 1:size(M_env,1)
                    plot([M_env(i,1),M_env(i,3)],[M_env(i,2),M_env(i,4)],...
                        '-','Color',beam_colormap(round(255*beam_weight(i))+100,:));
                end
                %--- coverage area
                x_pts = [M_env(1,1);M_env(:,3);M_env(1,1)];
                y_pts = [M_env(1,2);M_env(:,4);M_env(1,2)];
                plot(x_pts-0.0,y_pts-0.0,'-','LineWidth',1.5,'color','b');

            end




            %=======================================
        case 'simulation'
            %=======================================


            % -- measurement matrix
            fprintf('- measurement matrix for line model.\n');

            % meas_files = dir([dir_.logs,'measurement*_timestmp*_hotspot*_adaptive*']);
            meas_files = 'measurement_conf*.dat*';
            files_names = dir([dir_.logs,meas_files])

            M = [];
            for i = 1:size(files_names,1)
                files_names(i).name;

                M_current = dlmread([dir_.logs,files_names(i).name]);
                M = [M;M_current];
            end

            %-- final M matrix (new)
            M = M(:,[3,4,5,6,8]);
    end


    start_x__all = M(:,1);
    start_y__all = M(:,2);
    end___x__all = M(:,3);
    end___y__all = M(:,4);
    weight___all = M(:,5)/max(M(:,5));
    
        
    max(M(:,5))
    
    
    
    %max_concentration_optical_beams = max(M(:,8))

    beam_colormap = flipud(bone(512));
    for i = 1:4:size(start_x__all,1)
        plot([start_x__all(i),end___x__all(i)]-0.0,...
             [start_y__all(i),end___y__all(i)]-0.0,...
             '-','Color',beam_colormap(round(255*weight___all(i))+100,:)); %100
    end
    %}
    
    %-- sensing positions
    %start_x__unq = unique(start_x__all)    
    [~,unq_ind,~] = unique(start_x__all)
    
    start_x__unq = start_x__all(unq_ind)
    start_y__unq = start_y__all(unq_ind)
    
    plot(start_x__unq,start_y__unq,...
        'o',...
        'Color',color_conf,...
        'MarkerFaceColor',color_conf,...
        'MarkerSize',marker_size_conf_position);
    
    pause
    
    % --- concentration scale.
    %{
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
    
    
    
    
    
    %----------------------------------------------------------------------------
    
  

end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             CANDIDATE CONFIGURATIONS MAP
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishCandConf
    
    map_candConf = dataPublish.map_candConf;
    
    candConf_ind = find(map_candConf);
    [candConf_r,candConf_c] = ind2sub(size(map_env),candConf_ind);
    plot(candConf_r,candConf_c,...
        'x',...
        'Color',[102,178,255]/255,...
        'MarkerSize',marker_size_cand_conf);
end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             HOTSPOTS:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishHotspots
    
    hotspots = dataPublish.hotspots;
    
    [hot_r,hot_c] = ind2sub(size(map_env),hotspots);
    h_hot = plot(hot_r,hot_c,...
        'ok',...
        'MarkerSize',marker_size_hotspot,...
        'LineWidth',line_width_hotspot,...
        'MarkerFaceColor',[1,0,0]);
    
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
        'Color',[0,0,0.75],...        'EdgeColor',[0.25,0.25,1.0],...
        'Interpreter','latex',...
        'Margin',1,...
        'FontWeight','bold');
end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                        MARK WEAK/REDUNDANT CONFIGURATIONS:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishWeakConf
    
    RedundantConf_ind = dataPublish.RedundantConfNum;
    RedundantConf_orn = dataPublish.RedundantConfOrn;
    
    redConfCell_num = RedundantConf_ind;
    [redConfCell_r,redConfCell_c] = ind2sub(size(map_env),redConfCell_num);
    
    % -- position
    plot(redConfCell_r,redConfCell_c,...
        'o',...
        'Color',color_conf_redundatn,...
        'MarkerFaceColor',color_conf_redundatn,...
        'MarkerSize',marker_size_conf_position+1);
    
    % -- arrow
    for i = 1:numel(RedundantConf_ind)
        
        quiver(redConfCell_r(i),redConfCell_c(i),...
               arrow_length*cosd(RedundantConf_orn(i)),...
               arrow_length*sind(RedundantConf_orn(i)),...
               'LineWidth',line_width_arrow+1,...
               'color',color_conf_redundatn,...
               'MaxHeadSize',arrow_head_size,...
               'MarkerSize',marker_size_arrow);
    end
    
    
    FoV = dataPublish.FoV;
    
    %arrow_length = 2;
    %line_length  = 2;    
    
    % -- symbolic field of view
    for i = 1:numel(RedundantConf_ind)
        
        % -- drat sumbolic field of view
        start_angle  = RedundantConf_orn(i)-FoV/2; %FoVfaces.ang(nCONF(i),1);
        end___angle  = RedundantConf_orn(i)+FoV/2; %start_angle+FoV;
        r            = symbolic_fov_length;
        x1           = r*cosd(start_angle);
        y1           = r*sind(start_angle);
        x2           = r*cosd(end___angle);
        y2           = r*sind(end___angle);
        h_wsfov = plot([redConfCell_r(i),x1+redConfCell_r(i)],...
             [redConfCell_c(i),y1+redConfCell_c(i)],...
              'color',color_conf_redundatn,...
              'LineWidth',line_width_symbolic_fov+1); 
        plot([redConfCell_r(i),x2+redConfCell_r(i)],...
             [redConfCell_c(i),y2+redConfCell_c(i)],...
              'color',color_conf_redundatn,...
              'LineWidth',line_width_symbolic_fov+1);
    end
    
end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                        ADAPTIVE (INITIALLY) PLANNED CONFIGURATIONS:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishAdaptivePlannedConf
    
    plannedConf_ind = dataPublish.AdaptivePlannedConf_ind;
    plannedConf_orn = dataPublish.AdaptivePlannedConf_orn;
    
    plannedConfCell_num = plannedConf_ind;
    [plannedConfCell_r,plannedConfCell_c] = ind2sub(size(map_env),plannedConfCell_num);
    
    plannedConf_seq = dataPublish.AdaptivePlannedConf_seq;
    
    
    % -- position
    plot(plannedConfCell_r,plannedConfCell_c,...
        'o',...
        'Color',color_conf_adaptive_planned,...
        'MarkerFaceColor',color_conf_adaptive_planned,...
        'MarkerSize',marker_size_conf_position);
    
    % -- arrow
    for i = 1:numel(plannedConf_ind)
        
        quiver(plannedConfCell_r(i),plannedConfCell_c(i),...
               arrow_length*cosd(plannedConf_orn(i)),...
               arrow_length*sind(plannedConf_orn(i)),...
               'LineWidth',line_width_arrow,...
               'color',color_conf_adaptive_planned,...
               'MaxHeadSize',arrow_head_size,...
               'MarkerSize',marker_size_arrow);
    end
    
    
    FoV = dataPublish.FoV;
    
    %arrow_length = 2;
    %line_length  = 2;    
    
    % -- symbolic field of view
    for i = 1:numel(plannedConf_ind)
        
        % -- drat sumbolic field of view
        start_angle  = plannedConf_orn(i)-FoV/2; %FoVfaces.ang(nCONF(i),1);
        end___angle  = plannedConf_orn(i)+FoV/2; %start_angle+FoV;
        r            = symbolic_fov_length;
        x1           = r*cosd(start_angle);
        y1           = r*sind(start_angle);
        x2           = r*cosd(end___angle);
        y2           = r*sind(end___angle);
        plot([plannedConfCell_r(i),x1+plannedConfCell_r(i)],...
             [plannedConfCell_c(i),y1+plannedConfCell_c(i)],...
              'color',color_conf_adaptive_planned,...
              'LineWidth',line_width_symbolic_fov); 
        plot([plannedConfCell_r(i),x2+plannedConfCell_r(i)],...
             [plannedConfCell_c(i),y2+plannedConfCell_c(i)],...
              'color',color_conf_adaptive_planned,...
              'LineWidth',line_width_symbolic_fov);
    end
    
    
    
    % -- sequence
    
    for i = 1:numel(plannedConf_ind)

        % conf number.
        x3 = plannedConfCell_r(i)-(distance_conf_seq_adaptive_planned*cosd(plannedConf_orn(i)));
        y3 = plannedConfCell_c(i)-(distance_conf_seq_adaptive_planned*sind(plannedConf_orn(i)));
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
                
        if     plannedConf_orn(i) >= 315 || plannedConf_orn(i) < 045
            HoriAlign = 'right';
            VertAlign = 'middle';
        elseif plannedConf_orn(i) >= 045 && plannedConf_orn(i) < 135
            HoriAlign = 'center';
            VertAlign = 'top'; 
        elseif plannedConf_orn(i) >= 135 && plannedConf_orn(i) < 255
            HoriAlign = 'left';
            VertAlign = 'middle';
        elseif plannedConf_orn(i) >= 255 && plannedConf_orn(i) < 315
            HoriAlign = 'center';
            VertAlign = 'bottom';
        else
            error('WTF');
        end
        
        text(x3,y3,num2str(plannedConf_seq(i)),...
            'HorizontalAlignment',HoriAlign,...
            'VerticalAlignment',VertAlign,...
            'FontWeight','demi',...
            'FontSize',font_size_conf_num,...
            'Color',color_conf_adaptive_planned,...
            'Interpreter','latex',...
            'EdgeColor',color_conf_adaptive_planned,...
            'Margin',1,...
            'FontWeight','bold');
        
    end
    
    % - hotspots
    hotspots_executed = dataPublish.AdaptiveHotspots;
    plot(hotspots_executed(:,1),hotspots_executed(:,2),...
        'ok',...
        'MarkerSize',marker_size_hotspot,...
        'MarkerFaceColor',[1,0,0]);
    
end



%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             MARK CONFIGURATIONS:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishConfPosition
    
    ConfNum = dataPublish.ConfNum;
    ConfCell_num = ConfNum;
    [ConfCell_r,ConfCell_c] = ind2sub(size(map_env),ConfCell_num);
    plot(ConfCell_r,ConfCell_c,...
        'o',...
        'Color',color_conf,...
        'MarkerFaceColor',color_conf,...
        'MarkerSize',marker_size_conf_position);
end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             DRAW FIELD OF VIEW:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishFoV
    
    FoV = dataPublish.FoV;
    SensorRange = dataPublish.SensorRange;
    
    ConfNum = dataPublish.ConfNum;
    ConfCell_num = ConfNum;
    [ConfCell_r,ConfCell_c] = ind2sub(size(map_env),ConfCell_num);
    ConfOrn = dataPublish.ConfOrn;
    
    for i = 1:length(ConfNum)
        
        % Draw FoV
        start_angle  = ConfOrn(i)-FoV/2; % first/start angle.
        sector_angle = FoV;              % increment in the first/start angle.
        MaxPts       = 1000;             % maximum points to plot the FoV.
        [xx,yy]      = SecDraw(start_angle,sector_angle,SensorRange.cell,MaxPts); % x,y points.
        plot(xx+ConfCell_r(i),yy+ConfCell_c(i),...
            'color',color_conf,...
            'LineWidth',line_width_fov);
        
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
    
    ConfNum = dataPublish.ConfNum;
    ConfCell_num = ConfNum;
    [ConfCell_r,ConfCell_c] = ind2sub(size(map_env),ConfCell_num);
    ConfOrn = dataPublish.ConfOrn;
    
    for i = 1:numel(ConfNum)
        
        % -- draw sumbolic field of view
        start_angle  = ConfOrn(i)-FoV/2; %FoVfaces.ang(nCONF(i),1);
        end___angle  = ConfOrn(i)+FoV/2; %start_angle+FoV;
        r            = symbolic_fov_length;
        x1           = r*cosd(start_angle);
        y1           = r*sind(start_angle);
        x2           = r*cosd(end___angle);
        y2           = r*sind(end___angle);
        h_sfov = plot([ConfCell_r(i),x1+ConfCell_r(i)],...
             [ConfCell_c(i),y1+ConfCell_c(i)],...
              'color',color_conf,'LineWidth',line_width_symbolic_fov); 
        plot([ConfCell_r(i),x2+ConfCell_r(i)],...
             [ConfCell_c(i),y2+ConfCell_c(i)],...
              'color',color_conf,'LineWidth',line_width_symbolic_fov); 
    end
end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             DRAW CONF ARROWS:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishArrow
    
    ConfNumbers      = dataPublish.ConfNum;
    ConfOrn = dataPublish.ConfOrn;
    
    [ConfCell_r,ConfCell_c] = ind2sub(size(map_env),ConfNumbers);
    
    for i = 1:length(ConfNumbers)
        
        conf_to_plot = [ConfNumbers(i),ConfCell_r(i),ConfCell_c(i),ConfOrn(i)];
        
        %FoVfaces.lgt
        %nCONF(i)
        %FoVfaces.lgt(nCONF(i))
        
        % Draw arrow
        quiver(ConfCell_r(i),ConfCell_c(i),...
               arrow_length*cosd(ConfOrn(i)),...
               arrow_length*sind(ConfOrn(i)),...
               'LineWidth',line_width_arrow,...
               'color',color_conf,...
               'MaxHeadSize',arrow_head_size,... %'Marker','o',...
               'MarkerSize',marker_size_arrow);
        %pause       
        
        %text(ConfCell_r(i),ConfCell_c(i),'\rightarrow','FontSize',15,'Rotation',ConfOrn(i));        
        %pause
        
        %arrow([ConfCell_r(i),ConfCell_c(i)],...
        %    [ConfCell_r(i)+arrow_length*cosd(ConfOrn(i)),ConfCell_c(i)+arrow_length*sind(ConfOrn(i))])
        %pause
        
    end
end

%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             PRINT CONF SEQUENCE NUMBERS:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishConfNum
    
    
    confSequenceNum = dataPublish.confSequenceNum;
    
    ConfNum = dataPublish.ConfNum;
    ConfOrn = dataPublish.ConfOrn;
    
    [ConfCell_r,ConfCell_c] = ind2sub(size(map_env),ConfNum);
    
    %all_conf = zeros(length(iC),4); 
    
    for i = 1:length(ConfNum)

        % conf number.
        x3 = ConfCell_r(i)-(distance_conf_seq*cosd(ConfOrn(i)));
        y3 = ConfCell_c(i)-(distance_conf_seq*sind(ConfOrn(i)));
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
                
        if     ConfOrn(i) >= 315 || ConfOrn(i) < 045
            HoriAlign = 'right';
            VertAlign = 'middle';
        elseif ConfOrn(i) >= 045 && ConfOrn(i) < 135
            HoriAlign = 'center';
            VertAlign = 'top'; 
        elseif ConfOrn(i) >= 135 && ConfOrn(i) < 255
            HoriAlign = 'left';
            VertAlign = 'middle';
        elseif ConfOrn(i) >= 255 && ConfOrn(i) < 315
            HoriAlign = 'center';
            VertAlign = 'bottom';
        else
            error('WTF');
        end
        
        if confSequenceNum(i) == 12
            HoriAlign = 'right';
            VertAlign = 'middle';
        end
           
        
        text(x3,y3,num2str(confSequenceNum(i)),...
            'HorizontalAlignment',HoriAlign,...
            'VerticalAlignment',VertAlign,...
            'FontWeight','demi',...
            'FontSize',font_size_conf_num,...
            'Color',[0,0,0.75],...
            'Interpreter','latex',...            'EdgeColor',[0.25,0.25,1.0],...
            'Margin',1,...
            'FontWeight','bold');
        
                
        %this_conf = [confSequenceNum(i),nCELLr(i), nCELLc(i), FoVfaces.lgt(nCONF(i))];
        %all_conf(i,:) = this_conf;
    end
    
    %pause
    % all_conf
    %all_conf = sortrows(all_conf,1)

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
%                             REDUNDANT CONF ARROWS:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishRedundantConfArrow
    
    
    redundantConfNum = dataPublish.RedundantConfNum;
    redundantConfOnt = dataPublish.RedundantConfOrn;
    
    [redundantConfCell_r,redundantConfCell_c] = ind2sub(size(map_env),redundantConfNum);
    
    for i = 1:length(redundantConfNum)
        
        % Draw arrow
        quiver(redundantConfCell_r(i),redundantConfCell_c(i),...
               2*cosd(redundantConfOnt(i)),...
               2*sind(redundantConfOnt(i)),...
               'LineWidth',line_width_arrow,...
               'color','m',...
               'MaxHeadSize',arrow_head_size,...
               'Marker','o',...
               'MarkerSize',marker_size_arrow);
        %pause
    end
end


%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             FUSED CONF ARROWS:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishFusedConfArrow
    
    
    fusedNum = dataPublish.FusedConfNum;
    fusedOnt = dataPublish.FusedConfOrn;
    
    [fusedCell_r,fusedCell_c] = ind2sub(size(map_env),fusedNum);
    
    for i = 1:length(fusedNum)
        
        % Draw arrow
        quiver(fusedCell_r(i),fusedCell_c(i),...
               2.5*cosd(fusedOnt(i)),...
               2.5*sind(fusedOnt(i)),...
               'LineWidth',line_width_arrow,...
               'color','g',...
               'MaxHeadSize',arrow_head_size,...
               'Marker','o',...
               'MarkerSize',marker_size_arrow);
        %pause
    end
end

%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                            IN-OUT TSP ANGLES FOR HOTSPOTS:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
if para_.PublishHotspotsInOutTSPAngles
    
    
    
    InOutTSPAnglesMean = dataPublish.InOutTSPAnglesMean;
    
    
    
    for i = 1:size(InOutTSPAnglesMean,1)
        
        hotspotCell_ind = hotspots(i);
        [hotspotCell_r,hotspotCell_c] = ind2sub(size(map_env),hotspotCell_ind);
        
        % entrance mark
        entranceMark_r = hotspotCell_r+(4*cosd(InOutTSPAnglesMean(i,1)));
        entranceMark_c = hotspotCell_c+(4*sind(InOutTSPAnglesMean(i,1)));
        
        % angle from enterance mark to hotspot
        angleMarkToHotspot = rad2deg(atan2(hotspotCell_c-entranceMark_c,hotspotCell_r-entranceMark_r));
        
        % exit mark
        exitMark_r = hotspotCell_r+(0.25*cosd(InOutTSPAnglesMean(i,2)));
        exitMark_c = hotspotCell_c+(0.25*sind(InOutTSPAnglesMean(i,2)));
        
       
        % plot entrance arrow
        quiver(entranceMark_r,entranceMark_c,...
            4*cosd(angleMarkToHotspot),...
            4*sind(angleMarkToHotspot),...
            'LineWidth',line_width_arrow,...
            'color','g',...
            'MaxHeadSize',arrow_head_size,...
            'MarkerSize',marker_size_arrow);
        
        % plot entrance arrow
        quiver(exitMark_r,exitMark_c,...
            4*cosd(InOutTSPAnglesMean(i,2)),...
            4*sind(InOutTSPAnglesMean(i,2)),...
            'LineWidth',line_width_arrow,...
            'color','g',...
            'MaxHeadSize',arrow_head_size,...
            'MarkerSize',marker_size_arrow);
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
    plot(refDistPoint(:,1),refDistPoint(:,2),...
        'db','MarkerSize',marker_size_hotspot,...
        'MarkerFaceColor',[0,1,0]);
end




%/////////////////////////////////////////////////////////////////
%                           LEGENDS
%/////////////////////////////////////////////////////////////////
% hLegend = legend([h_hot,h_wsfov,h_sfov],...
%     'Hotspots',...
%     'Local solutions',...
%     'Global (fused) solution',...
%     'position',[360,090,.25,.25],'Orientation','vertical');
% % 'location','SouthEast','Orientation','vertical');
% legend('boxoff');
% set(hLegend,'Interpreter','latex');
% %legend('boxoff');
% set([hLegend,gca],'FontSize',12);




%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                       EXPLORATION PLAN
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤

if para_.PublishExploration
    
    
    %/////////////////////////////////////////////////////////////////
    % APPROXIMATED TRAVELING ROUTE
    %/////////////////////////////////////////////////////////////////
    tRoute = dataPublish.tRoute;
    TRoute = dataPublish.TRoute;
    
    f = 4
    
    % --- traveling path
    % tRoute(end) = [];
    iC = tRoute;

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


    % --- detailed traveling path ---
    % TRoute(end) = [];
    jC = TRoute;
    nCELLj = fix(jC/f)+(~~mod(jC,f));   % cell num for each conf num.
    % nCONFj = mod(jC,f)+(~mod(jC,f)*f);  % conf num within a cell.
    nCELLj(diff(nCELLj)==0) = []; % remove repeated 
    [nCELLrj,nCELLcj] = ind2sub(size(map_env),nCELLj); % row and col of cell num.
    nCell_rj = nCELLrj; % keep another variable (row num).
    nCell_cj = nCELLcj; % keep another variable (col num).

    % --- little detailed traveling path ---
    % tRoute_in_TRoute = [TRoute, zeros(size(TRoute))];
    % for i = 1:numel(tRoute)
    %     ind = find(tRoute(i)==TRoute);
    %     tRoute_in_TRoute(ind,2) = 1;
    % end

    sample_rate = 0.25;%0.05;
    tRoute_interpol = [];
    prev_ind2 = 1;
    for i = 1:numel(tRoute)-1
        i
        tRoute_interpol = [tRoute_interpol; tRoute(i)];

        %ind1 = find(tRoute(i+0)==TRoute)
        %ind1 = find(tRoute(i+0)==TRoute(prev_ind2))
        ind1 = prev_ind2;

        %ind2 = find(tRoute(i+1)==TRoute)
        ind2 = min(find(tRoute(i+1)==TRoute(ind1:end,1)))+ind1-1
        prev_ind2 = ind2;

        delta_ind = ind2-ind1

        num_samples = round(delta_ind*sample_rate)
        if (num_samples) > 2

            %interpols = randi([tRoute(i+0)+1,tRoute(i+1)-1],num_samples,1);
            interpol_ind = randi([ind1+1,ind2-1],num_samples,1);
            interpol_ind = unique(interpol_ind);
            tRoute_interpol = [tRoute_interpol; TRoute(interpol_ind,1)];
        end
        tRoute_interpol = [tRoute_interpol; tRoute(i+1)];

    end

    kC = tRoute_interpol;
    nCELLk = fix(kC/f)+(~~mod(kC,f));   % cell num for each conf num.
    % nCONFj = mod(jC,f)+(~mod(jC,f)*f);  % conf num within a cell.
    nCELLk(diff(nCELLk)==0) = []; % remove repeated 
    [nCELLrk,nCELLck] = ind2sub(size(map_env),nCELLk); % row and col of cell num.
    nCell_rk = nCELLrk; % keep another variable (row num).
    nCell_ck = nCELLck; % keep another variable (col num).
    
    
    % traveling route
    % ------------------
    %{
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
        plot(x-0.5,y-0.5,'-','color',[0.4 0.4 0.4],'LineWidth',2);
    end
    %}
    
    nCell_rk = nCell_rk+rand(numel(nCell_rk),1)*0.1;
    nCell_ck = nCell_ck+rand(numel(nCell_ck),1)*0.1;

    %pt = interparc(1000,nCell_r,nCell_c,'spline');
    %pt = interparc(numel(nCell_r)*10,[nCell_r;nCell_r(1)],[nCell_c;nCell_c(1)],'spline');
    
    %pt = interparc(numel(nCell_rj)*2,[nCell_rj;nCell_rj(1)],[nCell_cj;nCell_cj(1)],'spline');
    pt = interparc(numel(nCell_rk)*4,...
                   [nCell_rk;nCell_rk(1)],[nCell_ck;nCell_ck(1)],...
                   'spline');
    %pt = interparc(numel(nCell_rk)*15,...
    %               [nCell_rk;nCell_rk(1)],[nCell_ck;nCell_ck(1)],...
    %               'csape');
    
    % Note: interparc interpolates between first and the last points,
    % intermediate points are used for the interpolation but the result may
    % not contain the exact intermedidate points. 
    
    %h_travel = plot(pt(:,1),pt(:,2),'ob','MarkerFaceColor',[119,136,153]/255,'MarkerSize',1.5); %2    
    h_travel = plot(pt(:,1),pt(:,2),...
        'o',...
        'Color',[0,0.5,0],...
        'MarkerFaceColor',[0,1,0],...
        'MarkerSize',1.5); %2
    
    
    
    
    %/////////////////////////////////////////////////////////////////
    %                           OPTICAL RAYS
    %/////////////////////////////////////////////////////////////////
    
    files = dir([dir_.logs,para_.MeasurementFiles])
    M = [];
    for i = 1:size(files,1)
        
        file_name = files(i).name
        current_M = dlmread([dir_.logs,file_name]);
        
        max(current_M(:,8))

        %M = [M;current_M];
        
        
        start_x__all = current_M(:,3);
        start_y__all = current_M(:,4);
        end___x__all = current_M(:,5);
        end___y__all = current_M(:,6);
        %weight___all = M(:,8)/max(M(:,8));
        
        %beam_colormap = flipud(bone(512));
        for j = 1:5:size(start_x__all,1)
            h_beam = plot([start_x__all(j),end___x__all(j)]-0.5,...
                 [start_y__all(j),end___y__all(j)]-0.5,...
                 '-','Color','r'); %beam_colormap(round(255*weight___all(i))+100,:));
        end
        
    end
    
    %/////////////////////////////////////////////////////////////////
    %                       SENSING POSITIONS
    %/////////////////////////////////////////////////////////////////
    
    ConfNum = dataPublish.ConfNum;
    ConfCell_num = ConfNum;
    [ConfCell_r,ConfCell_c] = ind2sub(size(map_env),ConfCell_num);
    h_pos = plot(ConfCell_r-0.5,ConfCell_c-0.5,...
        'o',...
        'Color',color_conf,...
        'MarkerFaceColor',color_conf,...
        'MarkerSize',marker_size_conf_position);
    
    
    ConfOrn = dataPublish.ConfOrn;
    for i = 1:length(ConfOrn)       
        
        % Draw arrow
        %h_pos = quiver(ConfCell_r(i)-0.5,ConfCell_c(i)-0.5,...
        %       (arrow_length-0.5)*cosd(ConfOrn(i)),...
        %       (arrow_length-0.5)*sind(ConfOrn(i)),...
        %       'LineWidth',line_width_arrow,...
        %       'color',color_conf,...
        %       'MaxHeadSize',arrow_head_size,... 
        %       'Marker','o',...
        %       'MarkerSize',marker_size_conf_position);
        %pause
    end
    
    
    %/////////////////////////////////////////////////////////////////
    %                           HOTSPOTS
    %/////////////////////////////////////////////////////////////////
    hotspots = dataPublish.hotspots;
    
    [hot_r,hot_c] = ind2sub(size(map_env),hotspots);
    h_hot = plot(hot_r,hot_c,...
        'ok',...
        'MarkerSize',marker_size_hotspot,...
        'MarkerFaceColor',[1,0,0]);
    
    
    %/////////////////////////////////////////////////////////////////
    %                           LEGENDS
    %/////////////////////////////////////////////////////////////////
    %hLegend = legend([h_travel,h_pos,h_hot,h_beam],...
    %    'Traveling path',...
    %    'Sensing positions',...
    %    'Hotspots',...
    %    'Optical beams',...
    %    'location','SouthEast','Orientation','vertical');
    hLegend = legend([h_hot,h_pos,h_travel,h_beam],...
        'Hotspots',...
        'Sensing positions',...
        'Traveling path',...        
        'Optical beams',...
        'location','SouthEast','Orientation','vertical');
    set(hLegend,'Interpreter','latex');
    %legend('boxoff');
    set([hLegend,gca],'FontSize',7);
    
    
end
    

% hold off
% alpha(.1)

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(2), pos(4)])

end

% ----------------------- End of the document ----------------------------------

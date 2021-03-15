function fVideo_1trms( dataPublish,para_,dir_ )
%fPublish.
% 
% 

% FoV = dataPublish.FoV;
% SensorRange = dataPublish.SensorRange;

%============================================
%           PUBLISH PARAMETERS
%============================================
if strcmp(para_.PublishTextWidth,'single-column') % -- single column.
    
    paraPub_.zoom_in_environment = 750;
    paraPub_.cell_marker_size = 4.85; %4
    paraPub_.color_scale = 'linear'; %
    %beam_color = [255,204,153]/256; % [204,229,255]/256
    paraPub_.marker_size_cand_conf = 2.5;
    paraPub_.marker_size_hotspot = 5; %7;
    %marker_size_weak_conf = 4;
    %marker_size_cond_conf = 5;
    paraPub_.marker_size_conf_position = 2; %5; % icra-submission size: 3
    paraPub_.line_width_fov = 1; %2.5
    paraPub_.symbolic_fov_length = 2.0; %2.0; %3;
    paraPub_.line_width_symbolic_fov = 1.0; %0.5; %1.25;

    paraPub_.arrow_length = 2.0; %3;
    paraPub_.line_width_arrow = 1.0; %0.5; %0.25
    paraPub_.arrow_head_size = 20;
    paraPub_.marker_size_arrow = 1; %2

    paraPub_.font_size_conf_num = 6;%16
    paraPub_.font_size_conc_scale = 10;
    %conf_num_dist = 0.5; % 0.3
    %conf_num_dist_adaptive_planned = 0.75;
end

paraPub_.colorConf         = [0.0 0.0 0.6];
paraPub_.colorExploration  = [065,105,225]/256; %[0.0,0.0,0.7];
paraPub_.colorExploitation = [210,105,030]/256; %[0.7,0.0,0.0];
paraPub_.colorPath         = [0.5,0.5,0.5];



%============================================
% initialize VIDEO RECORDING:
%============================================
file_name = sprintf('%s-%s-video.avi',...
    para_.ExplorationStrategy,para_.ExperimentTitle);
% thisVideo = VideoWriter([dir_.Video,file_name]);
thisVideo = VideoWriter(file_name);
thisVideo.FrameRate = 30;  % Default 30
thisVideo.Quality   = 100;
open(thisVideo);

%============================================
%-- environment map
%============================================
% map_env = dataPublish.map_env;
map_filename = sprintf('%s_map_coverage.dat',para_.environment);
map_env = load([dir_.env,map_filename]);
para_.PublishFigureTitle = "Some";
fPublish_Environment(map_env,para_,paraPub_)

%-- write this frame to video
%---------------------------------
frame = getframe; writeVideo(thisVideo,frame);

% pause video
%----------------------
% pauseVideo = 1;
% for dummyFrames = 1:round(pauseVideo*thisVideo.FrameRate)
%     frame = getframe;
%     writeVideo(thisVideo,frame);
% end
fPauseVideo(1,thisVideo)



%============================================
% 
%============================================
org_env_filename = sprintf('%s_origin_coverage.dat',para_.environment);
initial_position = load([dir_.env,org_env_filename]);
conf_filename = 'executed_confs_final.txt'; 

%-------
% dir
%-------
% dir_.Solutions = dir_.DetectionResults;
% dir_.logs = dir_.DetectionLogs;
confs = load([dir_.Solutions,conf_filename]);


%------------------------------------------------
% TRAVELING PATH INITIAL POSITION TO FIRST CONF
%------------------------------------------------
start_position = initial_position; 
end___position = confs(1,1:2);
[thisVideo] = fPublish_PathTwoPoints( start_position,...
                                      end___position,...
                                      thisVideo,...
                                      map_env,...
                                      paraPub_);
    
% -- optical beams param
% ----------------------------    
beam_colormap = flipud(bone(512));
beam_maxcolor = 300;
beam_mincolor = 100;
beam_maxconce = 5000;

hSymbol = [];

for i = 1:size(confs,1)-1
    
    %--------------------------
    % CONFIGURATION
    %--------------------------
    
    % -- retrive integral measurments
    % ----------------------------
    measu_file = sprintf('measurement_conf%02d.dat',i);
    M = dlmread([dir_.logs,measu_file]);    
    start_x = M(:,3);
    start_y = M(:,4);
    end___x = M(:,5);
    end___y = M(:,6);
    ppm___m = M(:,8);
    
    % -- coverage area
    % ----------------------------    
    %x_pts = [start_x(1),(end___x)',start_x(1)];
    %y_pts = [start_y(1),(end___y)',start_y(1)];
    %
    [xx,yy] = SecDraw( confs(i,3)-(confs(i,4)/2),...
                       confs(i,4),...
                       confs(i,5),...
                       100);
    x_pts = xx+confs(i,1); y_pts = yy+confs(i,2);
    
    %-- plot coverage area
    %-----------------------------
    if  confs(i,7) == 0
        conf_color = paraPub_.colorExploration;
        planType   = 'Exploration';
        planPosition = [65,05];
    elseif confs(i,7) == 1
        conf_color = paraPub_.colorExploitation;
        planType   = 'Exploitation';
        planPosition = [80,05];
    end
    hFoV = plot(x_pts-0.5,y_pts-0.5,'-','LineWidth',1,'color',conf_color);    
    hPosition = plot(start_x(1)-0.5,start_y(1)-0.5,...
        'o','color',conf_color,'MarkerFaceColor',conf_color);
    hPlanType = text(planPosition(1),planPosition(2),planType,...
                    'FontWeight','demi',...
                    'FontSize',12,...
                    'Color','k',...
                    'Interpreter','latex',...
                    'FontWeight','bold');
    
    
    %-- write this frame to video
    %---------------------------------
    writeVideo(thisVideo,getframe);
    fPauseVideo(1,thisVideo)
    
    % -- optical beams
    % ----------------------------
    % {
    for j = 1:10:size(start_x,1)
        this_colo_ind=...
            round((min([beam_maxconce,ppm___m(j)])/beam_maxconce)*...
            beam_maxcolor)+beam_mincolor;
        %hBeam(j) = plot([start_x(j)-0.5,end___x(j)-0.5],...
        %                [start_y(j)-0.5,end___y(j)-0.5],...
        %                '-','color',beam_colormap(this_colo_ind,:));
        hBeam(j) = plot([start_x(j)-0.5,end___x(j)-0.5],...
                        [start_y(j)-0.5,end___y(j)-0.5],...
                        '-','color','r');
        %h_beam.Color(4) = 0.5;
        
        %-- write this frame to video
        %---------------------------------
        writeVideo(thisVideo,getframe);
    end
    %} 
    
    %-- 
    %---------------------------------
    delete(hFoV);
    delete(hPosition);
    delete(hBeam);
    delete(hPlanType);
    
    %========================================
    %       RECONSTRUCTION
    %========================================
    
    [ mean_map,...
      cell_coords_x,...
      cell_coords_y ] = fReconstructionVideo( i,...
                                              para_,...
                                              dir_);
                                          
    %x_coord_file = sprintf('x_coord_conf%02d.dat',i);
    %y_coord_file = sprintf('y_coord_conf%02d.dat',i);
    %recon_file   = sprintf('reconstruction_conf%02d.dat',i);
    
    %cell_coords_x = load([dir_.recon,x_coord_file]);
    %cell_coords_y = load([dir_.recon,y_coord_file]);
    %mean_map = load([dir_.recon,recon_file]);    
    %map_recon = mean_map';
    if i>1; delete(hRecon); end
    map_recon = mean_map';
    [ hRecon ] = fPublish_Reconstruction( map_recon,...
                                     cell_coords_x,...
                                     cell_coords_y,...
                                     map_env,...
                                     para_,...
                                     paraPub_ );
    %--- conf symbol
    %-------------------------------------------
    if i>1
        delete(hSymbol.line1);
        delete(hSymbol.line2);
        delete(hSymbol.posit);
        delete(hSymbol.arrow);
    end
    [ hSymbol ] = fConfSymobls(confs,i,paraPub_);
 
    %-- write video
    writeVideo(thisVideo,getframe);    
    
    %----------------------------
    % TRAVELING PATH
    %----------------------------    
    %[thisVideo] = fPubPath(confs,i,thisVideo,map_env,paraPub_);    
    start_position = confs(i+0,1:2); 
    end___position = confs(i+1,1:2);
    [thisVideo] = fPublish_PathTwoPoints(start_position,...
                                    end___position,...
                                    thisVideo,...
                                    map_env,...
                                    paraPub_);
end


%--------------------------
% LAST CONFIGURATION
%--------------------------

% -- retrive integral measurments
% ----------------------------
measu_file = sprintf('measurement_conf%02d.dat',i+1);
M = dlmread([dir_.logs,measu_file]);    
start_x = M(:,3);
start_y = M(:,4);
end___x = M(:,5);
end___y = M(:,6);
ppm___m = M(:,8);
    
% -- coverage area
% ----------------------------    
[xx,yy] = SecDraw(confs(i+1,3)-(confs(i+1,4)/2),...
                  confs(i+1,4),...
                  confs(i,5),...
                  100);
x_pts = xx+confs(i+1,1); y_pts = yy+confs(i+1,2);
    
%-- plot coverage area
%-----------------------------
if  confs(i+1,7) == 0
    conf_color = paraPub_.colorExploration;
    planType   = 'Exploration';
    planPosition = [65,05];
elseif confs(i+1,7) == 1
    conf_color = paraPub_.colorExploitation;
    planType   = 'Exploitation';
    planPosition = [80,05];
end
hFoV = plot(x_pts-0.5,y_pts-0.5,'-','LineWidth',1,'color',conf_color);    
hPosition = plot(start_x(1)-0.5,start_y(1)-0.5,...
    'o','color',conf_color,'MarkerFaceColor',conf_color);
hPlanType = text(planPosition(1),planPosition(2),planType,...
                'FontWeight','demi',...
                'FontSize',12,...
                'Color','k',...
                'Interpreter','latex',...
                'FontWeight','bold');
                
%-- write this frame to video
%---------------------------------
writeVideo(thisVideo,getframe);
fPauseVideo(1,thisVideo)
    
% -- optical beams
% ----------------------------
% {
for j = 1:10:size(start_x,1)
    hBeam(j) = plot([start_x(j)-0.5,end___x(j)-0.5],...
                    [start_y(j)-0.5,end___y(j)-0.5],...
                    '-','color','r');
    %-- write this frame to video
    %---------------------------------
    writeVideo(thisVideo,getframe);
end
   
delete(hFoV);
delete(hPosition);
delete(hBeam);
delete(hPlanType);

%========================================
%       RECONSTRUCTION
%========================================
[ mean_map,...
  cell_coords_x,...
  cell_coords_y ] = fReconstructionVideo( i+1,...
                                          para_,...
                                          dir_);
delete(hRecon);
map_recon = mean_map';
[ hRecon ] = fPublish_Reconstruction( map_recon,...
                                 cell_coords_x,...
                                 cell_coords_y,...
                                 map_env,...
                                 para_,...
                                 paraPub_ );
%--- conf symbol
%-------------------------------------------
delete(hSymbol.line1);
delete(hSymbol.line2);
delete(hSymbol.posit);
delete(hSymbol.arrow);
[ hSymbol ] = fConfSymobls(confs,i+1,paraPub_);
    
%-- write video
writeVideo(thisVideo,getframe);    

%----------------------------
% TRAVELING PATH
%----------------------------    
start_position = confs(i+1,1:2); 
end___position = initial_position;
[thisVideo] = fPublish_PathTwoPoints(start_position,...
                                end___position,...
                                thisVideo,...
                                map_env,...
                                paraPub_);
%------------------------------                            
%--- conf symbol
%------------------------------
delete(hSymbol.line1);
delete(hSymbol.line2);
delete(hSymbol.posit);
delete(hSymbol.arrow);
[ hSymbol ] = fConfSymobls(confs,i+1,paraPub_);

fPauseVideo(40,thisVideo)
end








% ----------------------- End of the document ----------------------------------

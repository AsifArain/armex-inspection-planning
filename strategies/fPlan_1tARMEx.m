function [ ] = fPlan_1tARMEx( FoV,...
                              cellsize_env,...
                              retrieve1SE_pts,...
                              dir_,...
                              para_,...
                              visualize )

%
%

para_.Publish_1tARMEx_CoverageArea      = 0;
para_.Publish_1tARMEx_OpticalBeams      = 0;
para_.Publish_1tARMEx_ScannedArea       = 0;
para_.Publish_1tARMEx_OrientationArrow  = 1;
para_.Publish_1tARMEx_SymbolicFOV       = 1;
para_.Publish_1tARMEx_ConfSequenceNum   = 1;
para_.Publish_1tARMEx_ConfPosition      = 1;
para_.Publish_1tARMEx_Hotspots          = 1;

%////////////////////////////////////////////////////////
% 
%    RETRIVE DATA:
% 
%////////////////////////////////////////////////////////


file__map_env = sprintf('%s_map_coverage.dat',para_.environment);
map_env = load([dir_.env,file__map_env]);


%///////// preprocessing /////////////
load([dir_.sppdetection,'preprocess_detection.mat'],...
        'SensorRange',... 
        'OBSenv',...
        'OBSconf',...
        'V',...
        'oV',...
        'FoVfaces',...
        'confKept',...
        'T',...
        'E',...
        'o');

%/////////////// detection solution //////////////////
load([dir_.sppdetection,'detection_spp_initial.mat'],'C')

%/////////////// TSP detection solution //////////////////
% load([dir_.Solutions,'detection_initial_tsp.mat'],'confSequenceNum');

switch para_.TravelingType1SE 
    case 'TSP2'
        load([dir_.sppdetection,'detection_tsp_initial.mat'],'confSequenceNum');
    case 'ExploitationEnabledSort'
        load([dir_.sppdetection,'detection_ees_initial.mat'],'confSequenceNum');
end
    
%//////////////// DIRECTORIES
% dir_.logs = dir_.MeasurementLogs;
dir_.Recon = dir_.recon;



f = FoVfaces.num;

%////////////////////////////////////////////////////////
% 
%   INITIAL PLAN FOR GAS DETECTION
% 
%////////////////////////////////////////////////////////

[confs_detection_initial] = fExecutableConfsGD( FoVfaces,...
                                                SensorRange,...
                                                FoV,...
                                                oV,...
                                                C,...
                                                confSequenceNum,...
                                                confKept,...
                                                map_env,...
                                                dir_,...
                                                para_ );
                                                
%//////////////////////////////////////////////////////////
% 
% PLOT MAPS
% 
%//////////////////////////////////////////////////////////

fPlotMaps(map_env,dir_);
% fPlotDetectionPlan(confs_detection_initial,FoV)
coverageCells = [];
fPublish_GDPlan1t( confs_detection_initial,...
                   FoV,...
                   SensorRange,...
                   coverageCells,...
                   map_env,...
                   para_ );


%//////////////////////////////////////////////////////////
%
% INITIALIZE THE EXPLORATION
%
%////////////////////////////////////////////////////////// 
confs_planned = confs_detection_initial;
confs_executed = [];
conf_sequence_num = 1;
hcs_num = 1;
hcs_list_sub = [];
hc_prev = [];
cellsToBeScanned = find(sum(oV,2));
hc_block_num = 1;


while size(confs_planned,1) >= 1
    
    
    %-- Data backup/retrieval point
    %--------------------------------------
    pt_num = 1;
    if any(retrieve1SE_pts == pt_num)
        %fprintf("\nData is being retrived at point# %02d....\n",pt_num);
        fprintf("\nData retrieval point %02d....\n",pt_num);
        load([dir_.str,sprintf('LastSave1SE_Pnt%02d',pt_num)]);
        retrieve1SE_pts(retrieve1SE_pts == pt_num) = [];
    else
        %{
        %- delete file if exist, and create a new
        if exist([dir_.Solutions,'BackupPtsSequence1SE.dat'],'file')==2
            delete([dir_.Solutions,'BackupPtsSequence1SE.dat']);
        end
        %}
        fileID = fopen([dir_.str,'BackupPtsSequence1SE.dat'],'wt');
        fclose(fileID);        
        
        %fprintf("\nData is being backedup at point# %02d....\n",pt_num);
        fprintf("\nData backup point %02d....\n",pt_num);
        save([dir_.str,sprintf('LastSave1SE_Pnt%02d',pt_num)],...
            '-regexp', '^(?!(retrieve1SE_pts)$).')
        
        %dlmwrite([dir_.Solutions,'BackupPtsSequence1SE.dat'],...
        %     pt_num,'-append','delimiter',' ');
        
        dlmwrite([dir_.str,'BackupPtsSequence1SE.dat'],...
             pt_num,'delimiter',' ');
         
        fileID = fopen([dir_.str,'BackupPtsSequence1SE.dat'],'a'); 
        fclose(fileID);        
    end
    %pause
    
    
    
    
    
    %/////////////////////////////////////////////////////////////
    %
    % MOVE TO THE FIRST (PLANNED) CONF AND SAMPLE THE ENVIRONMENT
    %
    %/////////////////////////////////////////////////////////////
    sensing_conf = confs_planned(1,:);
    
    % conf for
    sensing_conf = [sensing_conf 0]
    
    [ cellsToBeScanned,...
      cellsScanned ] = fExecuteThisConf( sensing_conf,...
                                         conf_sequence_num,...
                                         map_env,...
                                         FoV,...
                                         SensorRange,...
                                         cellsToBeScanned,...
                                         para_,...
                                         dir_ );
    %pause                                  
    % -- update confs list 
    % ----------------------------
    %confs_executed     = [confs_executed;confs_planned(1,:)];
    confs_executed
    confs_executed     = [confs_executed;sensing_conf]
    %confs_executed     = [confs_executed;[confs_planned(1,:) 0]];
    confs_executed_hc  = confs_planned(1,:); 
    confs_planned(1,:) = [];


    % -- plot 
    % ----------------------------
    confFor = 'detection';
    [ ppmm ] = fPublish_ExecutedConf( conf_sequence_num,...
                                      sensing_conf,...
                                      cellsScanned,...
                                      SensorRange,...
                                      map_env,...
                                      confFor,...
                                      FoV,...
                                      para_,...
                                      dir_);
    
    
    %pause
    fprintf("\n      Max integral concentration beam: %d ppmm\n",max(ppmm));
    % -- initializations for this hc
    % ----------------------------
    conf_sequence_list_hc = conf_sequence_num;
    
    %/////////////////////////////////////////////////////////////
    %
    %	HIGH CONCENTRATION IS DETECTED
    %
    %/////////////////////////////////////////////////////////////
    if max(ppmm) > para_.ConcentrationThreshould4TomographyPPM %500
        
        
        %//////////////////////////////////////////////////////////
        %
        %               PLANNING FOR TOMOGRAPHY 
        %
        %//////////////////////////////////////////////////////////
        
        
        fprintf(1,'\n ********************************************************\n')
        fprintf(1,'    ADAPTIVE REPLANNING FOR GAS DETECTION AND GDM \n')
        fprintf(1,' ********************************************************\n\n')
        
                
        
        [  hcs_num,...
           hc_prev,...
           hc_block_num,...
           confs_executed,...
           conf_sequence_num,...
           confs_planned,...
           cellsToBeScanned,...
           retrieve1SE_pts,...
           hcs_list_sub ] = fAdaptiveReplanning1t( FoV,...
                                                   map_env,...
                                                   f,...
                                                   o,...
                                                   FoVfaces,...
                                                   E,...
                                                   T,...
                                                   oV,...
                                                   confKept,...
                                                   hcs_num,...
                                                   hc_prev,...
                                                   hc_block_num,...
                                                   SensorRange,...
                                                   confs_executed,...
                                                   confs_executed_hc,...
                                                   conf_sequence_list_hc,...
                                                   conf_sequence_num,...
                                                   confs_planned,...
                                                   cellsToBeScanned,...
                                                   hcs_list_sub,...
                                                   OBSenv,...
                                                   OBSconf,...
                                                   cellsize_env,...
                                                   retrieve1SE_pts,...
                                                   dir_,...
                                                   para_,...
                                                   visualize );
        
        %}
        %clc;
    else
        
        %//////////////////////////////////////////////////////////
        %
        %   ELSE, INCREASE THE CONF SEQUENCE NUMBER 
        %
        %//////////////////////////////////////////////////////////
        conf_sequence_num = conf_sequence_num+1;
        
    end
    
    
    %-- Data backup/retrieval point
    %--------------------------------------
    pt_num = 11;
    if any(retrieve1SE_pts == pt_num)
        fprintf("\nData retrieval point %02d....\n",pt_num);
        load([dir_.str,sprintf('LastSave1SE_Pnt%02d',pt_num)]);
        retrieve1SE_pts(retrieve1SE_pts == pt_num) = [];
    else
        fprintf("\nData backup point %02d....\n",pt_num);
        save([dir_.str,sprintf('LastSave1SE_Pnt%02d',pt_num)],...
            '-regexp', '^(?!(retrieve1SE_pts)$).')
        
        dlmwrite([dir_.str,'BackupPtsSequence1SE.dat'],...
             pt_num,'-append','delimiter',' ');
        fileID = fopen([dir_.str,'BackupPtsSequence1SE.dat'],'a'); 
        fclose(fileID);
        
    end
    %pause
    
end


%-- save plot
pause(2)
filename = sprintf('%sone-tour-strategy',dir_.str);
export_fig(filename, '-pdf','-r500')
pause(1)


%====================================================
%--- save executed conf to file
%====================================================
% filename = 'executed_confs_final.txt';        
% fileID = fopen([dir_.Solutions,filename],'wt'); fclose(fileID);
% % -- write integral concentration data to the file
% dlmwrite([dir_.Solutions,filename],confs_executed,'-append','delimiter',' ');
% fileID = fopen([dir_.Solutions,filename],'a'); fclose(fileID);

filename = 'executed_confs_final.txt';        
fileID = fopen([dir_.str,filename],'wt'); fclose(fileID);
% -- write integral concentration data to the file
dlmwrite([dir_.str,filename],confs_executed(:,[1:6 8]),...
    'delimiter',' ','precision','%06.2f');
fileID = fopen([dir_.str,filename],'a'); fclose(fileID);

%====================================================
%--- save executed conf to file
%====================================================
filename = 'hc_list_final.txt';        
fileID = fopen([dir_.str,filename],'wt'); fclose(fileID);
% -- write integral concentration data to the file
dlmwrite([dir_.str,filename],hcs_list_sub,'-append','delimiter',' ');
fileID = fopen([dir_.str,filename],'a'); fclose(fileID);





%profile off;
%profile off;
%profsave(profile('info'),['profile_tmp/prof']);





end





function fPlotMaps(map_env,dir_)
%
%

% -- environment
% -----------------------
figure; hold on;
% [map_row,map_col] = size(map_env);
% for i = 1:map_row
%     for j = 1:map_col
%         % if obstacle/free
%         if ~map_env(i,j) % occupied
%             col = [0.4 0.4 0.4];
%             faceCol = [0.4 0.4 0.4];
%         elseif map_env(i,j) % free
%             col = [0.9 0.9 0.9];
%             faceCol = [1 1 1];
%             %col = [1 1 1];
%         end
%         plot(i+0.5,j+0.5,'S','Color',col,...
%                          'MarkerFaceColor',faceCol,...
%                          'MarkerSize',20);
%     end
% end

imshow(map_env','InitialMagnification',350); %750 %1400
hold on; set(gca,'YDir','normal'); 
colormap('gray'); 
% caxis([-1 1.25])
caxis([-0.5 1])
%pause


% -- ground truth
% -------------------------

 map_con = dlmread([dir_.gt,'map_con.dat']);
 [map_row,map_col] = size(map_con);
% {
% gdm_colormap = (flipud(colormap('hot'))); % for main fig. 
gdm_colormap = (flipud((hot))); % for main fig. 
%gdm_colormap = (flipud((summer))); % for main fig. 
%gdm_colormap = flipud(autumn(512));

max_concentration = max(map_con(:)); % main fig.
delta = max_concentration/(size(gdm_colormap,1)-1);    
for i = 1:map_row
    for j = 1:map_col
        % if positive concentration
        if map_con(i,j)>100
            linear_color_number = round( min(map_con(i,j),max_concentration)/delta)+1;
            col = gdm_colormap( linear_color_number, :);
            if sum(col) == 3
                col = col-(1e-10);
            end
            % ----------------------------------------------------
            %plot(i+0.5,j+0.5,'s','Color',col,'MarkerFaceColor',col,'MarkerSize',15);
            plot(i+0.0,j+0.0,'s','Color',col,'MarkerFaceColor',col,'MarkerSize',2.95); %10
        end
    end
end
%}
 set(gca,'xtick',1:1:map_row);
 set(gca,'ytick',1:1:map_col);
grid on; axis equal;

%pause

end

    

%{
function fPlotDetectionPlan(confs_detection,FoV)

% --- initial coverage plan
% for i = 1:size(confs_coverage_initial,1)
%     start_angle  = confs_coverage_initial(i,3)-FoV/2; % first/start angle.
%     sector_angle = FoV;              % increment in the first/start angle.
%     MaxPts       = 1000;             % maximum points to plot the FoV.
%     [xx,yy]      = SecDraw(start_angle,sector_angle,SensorRange.cell,MaxPts); % x,y points.
%     plot(xx+confs_coverage_initial(i,1)-0.5,yy+confs_coverage_initial(i,2)-0.5,...
%         'color',[1.0,0.2,1.0],...
%         'LineWidth',2);
% end

for i = 1:size(confs_detection,1)
    % -- drat sumbolic field of view
    start_angle  = confs_detection(i,3)-FoV/2;
    end___angle  = confs_detection(i,3)+FoV/2;
    r            = 5;
    x1           = r*cosd(start_angle);
    y1           = r*sind(start_angle);
    x2           = r*cosd(end___angle);
    y2           = r*sind(end___angle);
    plot([confs_detection(i,1)-0.5,x1+confs_detection(i,1)-0.5],...
         [confs_detection(i,2)-0.5,y1+confs_detection(i,2)-0.5],...
          'color',[0.8,0.8,0.8],'LineWidth',3);      
    plot([confs_detection(i,1)-0.5,x2+confs_detection(i,1)-0.5],...
         [confs_detection(i,2)-0.5,y2+confs_detection(i,2)-0.5],...
          'color',[0.8,0.8,0.8],'LineWidth',3);
    plot(confs_detection(i,1)-0.5,confs_detection(i,2)-0.5,...
        'o','color',[0.8,0.8,0.8],'MarkerFacecolor',[0.8,0.8,0.8],'MarkerSize',8)
end

end
%}

% function [sensing_confs] = fExecutableConf( FoVfaces,...
%                                             SensorRange,...
%                                             FoV,...
%                                             oV,...
%                                             C,...
%                                             confSequenceNum,...
%                                             confKept,...
%                                             map_env,...
%                                             para_ )
% %
% % INITIAL COVERAGE PLAN FOR GAS DETECTION
% 
% f = FoVfaces.num;
% Conf = zeros(size(oV,2),1);
% Conf(confKept) = C;
% Conf_ind = find(Conf);                      % IDs of selected conf.
% Conf_cell = fix(Conf_ind/f)+(~~mod(Conf_ind,f));   % cell num for each conf num.
% [Conf_r,Conf_c] = ind2sub(size(map_env),Conf_cell); % row and col of cell num.
% Conf_o = FoVfaces.lgt(mod(Conf_ind,f)+(~mod(Conf_ind,f)*f));  % conf num within a cell.
% 
% % -- pose
% sensing_confs = [Conf_r+0.5,Conf_c+0.5,Conf_o];
% 
% % -- other parameters
% conf_fov = FoV;
% conf_r   = SensorRange.cell;
% conf_frq = FoV*2;
% conf_tim = para_.GroundTruthTimeStamp;
% 
% sensing_confs(:,4) = conf_fov*ones(size(sensing_confs,1),1);
% sensing_confs(:,5) = conf_r*  ones(size(sensing_confs,1),1);
% sensing_confs(:,6) = conf_frq*ones(size(sensing_confs,1),1);
% sensing_confs(:,7) = conf_tim*ones(size(sensing_confs,1),1);
% 
% 
% % -- conf by sequence
% sensing_confs = [confSequenceNum,sensing_confs];
% sensing_confs = sortrows(sensing_confs,1);
% sensing_confs = sensing_confs(:,2:end);
% 
% 
% 
% 
% 
% 
% 
% 
% end


%{
function [ dists ] = fPathDist(start_position,goal_positions,map_env )
%


Connecting_Distance = 1;
dists = zeros(size(goal_positions,1),1);

for i = 1:size(dists,1)
        
    StartX = round(start_position(2));
    StartY = round(start_position(1));
    MAP    = 1-map_env;
    
    GoalRegister = int8(zeros(size(MAP)));
    GoalRegister(round(goal_positions(i,1)),round(goal_positions(i,2))) = 1;
    
    %goal_positions(i,1)
    %goal_positions(i,2)    
    %GoalRegister    
    %[goal_r,goal_c] = find(GoalRegister)
    
    
    if pdist2(round(start_position),round(goal_positions(i,:)),'euclidean') >= 1
        
        %disp('im here')
    
        OptimalPath = ASTARPATH(StartX,StartY,MAP,GoalRegister,Connecting_Distance);
        dists(i) = size(OptimalPath,1)-1;
        
        %{
        hold on
        plot(OptimalPath(1,1),OptimalPath(1,2),'o','color','k')
        plot(OptimalPath(end,1),OptimalPath(end,2),'o','color','b')
        plot(OptimalPath(:,1),OptimalPath(:,2),'r')
        legend('Goal','Start','Path')
        %}
    
    else
        dists(i) = 0;
    end
    
    


    
    %pause
end



% %Start Positions
% StartX=15;
% StartY=15;
% 
% %Generating goal nodes, which is represented by a matrix. Several goals can be speciefied, in which case the pathfinder will find the closest goal. 
% %a cell with the value 1 represent a goal cell
% GoalRegister=int8(zeros(128,140));
% GoalRegister(110,80)=1;
% 
% %Number of Neighboors one wants to investigate from each cell. A larger
% %number of nodes means that the path can be alligned in more directions. 
% %Connecting_Distance=1-> Path can  be alligned along 8 different direction.
% %Connecting_Distance=2-> Path can be alligned along 16 different direction.
% %Connecting_Distance=3-> Path can be alligned along 32 different direction.
% %Connecting_Distance=4-> Path can be alligned along 56 different direction.
% %ETC......
% 
% Connecting_Distance=1; %Avoid to high values Connecting_Distances for reasonable runtimes. 
% 
% % Running PathFinder
% OptimalPath=ASTARPATH(StartX,StartY,MAP,GoalRegister,Connecting_Distance)


end
%}


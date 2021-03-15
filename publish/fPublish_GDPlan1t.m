function fPublish_GDPlan1t( confs_detection,...
                            FoV,...
                            SensorRange,...
                            coverageCells,...
                            map_env,...
                            para_ )

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

%col = [0.8,0.8,0.8];
% col = [0.0,0.75,0.0];
% col = [0.5,0.5,0.5];
%col = [0.5,0.5,0.5];
%col = [0.0,1.5,0.0];
% col = [102,255,102]/256;
col =[0, 0.5, 0];

for i = 1:size(confs_detection,1)
    
    
    if para_.Publish_1tARMEx_SymbolicFOV == 1
    
        % -- draw symbolic field of view
        start_angle  = confs_detection(i,3)-FoV/2;
        end___angle  = confs_detection(i,3)+FoV/2;
        r            = 2; %SensorRange.cell/6;
        x1           = r*cosd(start_angle);
        y1           = r*sind(start_angle);
        x2           = r*cosd(end___angle);
        y2           = r*sind(end___angle);
        

        %-- first line
        %-----------------------
        hLine1(i) = plot([confs_detection(i,1)-0.5,x1+confs_detection(i,1)-0.5],...
             [confs_detection(i,2)-0.5,y1+confs_detection(i,2)-0.5],...
              'color',col,'LineWidth',1.5); %3     
        %-- second line
        %-----------------------
        hLine2(i) = plot([confs_detection(i,1)-0.5,x2+confs_detection(i,1)-0.5],...
             [confs_detection(i,2)-0.5,y2+confs_detection(i,2)-0.5],...
              'color',col,'LineWidth',1.5); %3
          
    end
        
    %-- position
    %-----------------------
    if para_.Publish_1tARMEx_ConfPosition == 1
    
        hPoint(i) = plot(confs_detection(i,1)-0.5,confs_detection(i,2)-0.5,...
            'o','color',col,'MarkerFacecolor',col,'MarkerSize',2); %8
    end
    
    
    %-- orientation arrow
    %-----------------------
    if para_.Publish_1tARMEx_OrientationArrow == 1
        arrow_length = 3;
        line_width_arrow = 1.5; %2; %0.25
        arrow_head_size = 40; %20;
        marker_size_arrow = 2;  %1; %2
        hArrow(i) = quiver(confs_detection(i,1)-0.5,confs_detection(i,2)-0.5,...
               arrow_length*cosd(confs_detection(i,3)),...
               arrow_length*sind(confs_detection(i,3)),...
               'LineWidth',line_width_arrow,...
               'color',col,...
               'MaxHeadSize',arrow_head_size,... %'Marker','o',...
               'MarkerSize',marker_size_arrow);
           
           
           
%        quiver(confs_detection(i,1)-0.5,confs_detection(i,2)-0.5,...
%                arrow_length*cosd(confs_detection(i,3)),...
%                arrow_length*sind(confs_detection(i,3)),...
%             0,... % automatically scales the arrows to fit within the grid and then stretches them by the factor scale. scale = 2 doubles their relative length, and scale = 0.5 halves the length. Use scale = 0 to plot the velocity vectors without automatic scaling
%             'AutoScale','off',...
%             'MaxHeadSize',arrow_head_size*(3)/norm(arrow_length*cosd(confs_detection(i,3)),arrow_length*sind(confs_detection(i,3))),... headSize is the head size,but it cannot be longer than 1/3 of the vector length (for some reason)
%             'LineWidth',line_width_arrow,...
%             'LineStyle','-'); 
        
        %pause
    end
end

%-- coverage cells
%------------------------------
[cell_r,cell_c] = ind2sub(size(map_env),coverageCells);
hCoverage = plot(cell_r,cell_c,'xg','MarkerSize',2);




if strcmp(para_.SensingSystem,'simulation')
    %pause
end

% % pause(3)
% delete(hLine1);
% delete(hLine2);
% delete(hPoint);
delete(hCoverage);
% delete(hArrow);
% pause(1)


file_name = sprintf('fig/phd-final-seminar/%s-1t-armex-detection-plan',para_.ExperimentTitle);
export_fig([file_name], '-pdf','-r500')
% pause


if para_.Publish_1tARMEx_SymbolicFOV == 1
    if numel(hLine1) > 0
        delete(hLine1)
    end
    if numel(hLine2) > 0
        delete(hLine2)
    end
end
if para_.Publish_1tARMEx_ConfPosition == 1
    if numel(hPoint) > 0
        delete(hPoint)
    end
end
if para_.Publish_1tARMEx_OrientationArrow == 1
    if numel(hArrow) > 0
        delete(hArrow)
    end
end


end

function [ M_m,M_cell ] = fRobotLogs2M( measure_file,...
                                        para_,...
                                        dir_,...
                                        visulize )
%fExperimentLogs2M creates M matrix from measurements logs collected during
%experiments for either gas detection or gas tomography.
%   
% INPUTS:
% 
% OUTPUTS:
% M: M matrix to be used for generating gas distribution map.
% each row vector of M matrix contains initial point of an optical ray
% (x and y values in meters w.r.t the robot origin in the map), end point
% of the optical ray (x and y values in meters w.r.t the robot origin in
% the map), and gas concentration value (in ppm).

% measure_file

% m_logs = load([dir_.MeasurementLogs,measure_file]);
m_logs = load([dir_.logs,measure_file]);


map_filename      = sprintf('%s_map_raycast.dat',para_.environment);
origin_filename   = sprintf('%s_origin_raycast.dat',para_.environment);
cellsize_filename = sprintf('%s_cellsize_raycast.dat',para_.environment);


map_raycast      = load([dir_.env,map_filename]);
origin_raycast   = load([dir_.env,origin_filename]);
cellsize_raycast = load([dir_.env,cellsize_filename]);


% robot_origin_in_map_x = origin_raycast(1);
% robot_origin_in_map_y = origin_raycast(2);
x_origin_cell = origin_raycast(1);
y_origin_cell = origin_raycast(2);







if visulize == 1
    figure('name','xxx'); 
    imshow(map_raycast','InitialMagnification',800); hold on;
    set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
    plot(x_origin_cell,y_origin_cell,'ob');
end


conf_num         = m_logs(:,1);
conf_datetime    = m_logs(:,2);
robot_x          = m_logs(:,3);
robot_y          = m_logs(:,4);
robot_w          = m_logs(:,5);
pan_angle        = m_logs(:,6);
tilt_angle       = m_logs(:,7);
measurement      = m_logs(:,8);
% sensing_rng_cell = para_.SensingRangeM/dataYAML.resolution;
sensing_rng_cell = para_.SensingRangeM/cellsize_raycast;



M_m    = inf(size(m_logs,1),5);
M_cell = inf(size(m_logs,1),5);

for i = 1:size(m_logs,1)
    
    %-- start point (m)
    x_start_m    = robot_x(i);
    y_start_m    = robot_y(i);
    
    %-- start point (cell)
    %x_start_cell = (robot_x(i)/dataYAML.resolution)+x_origin_cell;
    %y_start_cell = (robot_y(i)/dataYAML.resolution)+y_origin_cell;
    x_start_cell = (robot_x(i)/cellsize_raycast)+x_origin_cell;
    y_start_cell = (robot_y(i)/cellsize_raycast)+y_origin_cell;
    
    
    %-- orientation
    %rad2deg(deg2rad(robot_w(i))+pan_angle(i));    
    robot_orientation = deg2rad(robot_w(i));
    pan_angle_this = pan_angle(i);
    
    %-- end points (m/cell)
    [ x_end_m,...
      y_end_m,...
      x_end_cell,...
      y_end_cell ] = fRayEndPointWithObstacles( x_start_cell,...
                                                y_start_cell,...
                                                sensing_rng_cell,...
                                                robot_orientation,...
                                                pan_angle_this,...
                                                x_origin_cell,...
                                                y_origin_cell,...
                                                cellsize_raycast,...
                                                map_raycast );
    %-- measurements
    ppm = measurement(i);
    
    %-- update M matrix
    M_m(i,:)    = [x_start_m,y_start_m,x_end_m,y_end_m,ppm];
    M_cell(i,:) = [x_start_cell,y_start_cell,x_end_cell,y_end_cell,ppm];
    
    
    if visulize == 1
        plot([x_start_cell,x_end_cell],...
             [y_start_cell,y_end_cell],...
             '-','color',rand(1,3))
        %pause
    end
    
end


% ***********************************************************************
% Post-approximation:
% ----------------------------
% integrate concentrations of beams with same intial and final positions
% due to approximation.
% ***********************************************************************

M_m_uniqueXY = unique(M_m(:,1:4),'rows');
M_cell_uniqueXY = unique(M_cell(:,1:4),'rows');

if size(M_m_uniqueXY,1) ~= size(M_cell_uniqueXY,1)
    warning("Unequal M matrix sizes in meters and cells");
end

for i = 1: size(M_m_uniqueXY,1)
    %{    
    M_m_uniqueXY(i,5) = ...
        sum( M_m( M_m(:,1) == M_m_uniqueXY(i,1) &...
                  M_m(:,2) == M_m_uniqueXY(i,2) &...
                  M_m(:,3) == M_m_uniqueXY(i,3) &...
                  M_m(:,4) == M_m_uniqueXY(i,4) ,    5) )
    %}        
    M_m_uniqueXY(i,5)    = sum( M_m(    all(M_m   (:,1:4) == M_m_uniqueXY   (i,1:4),2),5));
    M_cell_uniqueXY(i,5) = sum( M_cell( all(M_cell(:,1:4) == M_cell_uniqueXY(i,1:4),2),5));
end

% M_m(1:20,:)
% M_m_uniqueXY(1:20,:)
% 
% fileID = fopen('M_m_test.dat','wt'); fclose(fileID);
% dlmwrite('M_m_test.dat',M_m,'delimiter',' ');
% fileID = fopen('M_m_test.dat','a'); fclose(fileID);
% 
% 
% fileID = fopen('M_m_Utest.dat','wt'); fclose(fileID);
% dlmwrite('M_m_Utest.dat',M_m_uniqueXY,'delimiter',' ');
% fileID = fopen('M_m_Utest.dat','a'); fclose(fileID);
% 
% 
% fileID = fopen('M_cell_test.dat','wt'); fclose(fileID);
% dlmwrite('M_cell_test.dat',M_cell,'delimiter',' ');
% fileID = fopen('M_cell_test.dat','a'); fclose(fileID);
% 
% 
% fileID = fopen('M_cell_Utest.dat','wt'); fclose(fileID);
% dlmwrite('M_cell_Utest.dat',M_cell_uniqueXY,'delimiter',' ');
% fileID = fopen('M_cell_Utest.dat','a'); fclose(fileID);
% 
% pause


% Final M
%-----------------------------------
M_m    = M_m_uniqueXY;
M_cell = M_cell_uniqueXY;



% if visulize == 1
%     x_pts = [M_cell(1,1);M_cell(:,3);M_cell(1,1)];
%     y_pts = [M_cell(1,2);M_cell(:,4);M_cell(1,2)];
%     % y_pts = [start_y__all(1),(end___y__all)',start_y__all(1)];
%     plot(x_pts-0.0,y_pts-0.0,'-','LineWidth',1.5,'color','b');
% end

% SaveFile = ([M_LOC,M_FILE,'.mat']);
% save(SaveFile,'M');

end

function [ x_end_m,...
           y_end_m,...
           x_end_cell,...
           y_end_cell ] = fRayEndPointWithObstacles( x_start_cell,...
                                                     y_start_cell,...
                                                     sensing_rng_cell,...
                                                     robot_orientation,...
                                                     pan_angle_this,...
                                                     x_origin_cell,...
                                                     y_origin_cell,...
                                                     cellsize_raycast,...
                                                     map_raycast )
% fRayEndPointWithObstacles finds the point of projection for an optical
% ray for the given environment map (obstacle information), sensing range
% and angle of projection.
% 
% INPUTS:
% x_start, y_start: x and y start points of optical ray in meters -- w.r.t
% the robot origion in the map.
% sensing_rng: sensing range in meter.
% robot_orientation: robot orientation in radians
% pan_angle: pan angle in radians.
% robot_origin_in_map_x, robot_origin_in_map_y: cell number along x-axis
% and y-axis for the robot origion in the map.
% dataYAML.resolution: the real environment length (in meter) represented by a single
% pixel in the map.
% map_raycast: environment map contains obstacle information.
% 
% OUTPUTS:
% x_end,y_end: x and y end points of optical ray (in meters, w.r.t. the
% robot origion in the map).


% ------------------------------------------------------------------------------
% RAY CASTING
% ------------------------------------------------------------------------------

% -- angle of ray projection w.r.t to the map -- normalized between [0,360] degree
ProjectionAngle = rad2deg(normalizeAngle(robot_orientation+pan_angle_this));


% -- end limit point of projection of optical ray
x_endLimitPoint = x_start_cell + sensing_rng_cell*cosd(ProjectionAngle);
y_endLimitPoint = y_start_cell + sensing_rng_cell*sind(ProjectionAngle);
% endLimitPoint = [x_endLimitPoint,y_endLimitPoint];


% -- conditions to sample intermediate points along x or y-axis.
cond_x1 = ProjectionAngle>=000 && ProjectionAngle<=045;
cond_x2 = ProjectionAngle>=135 && ProjectionAngle<=225;
cond_x3 = ProjectionAngle>=315 && ProjectionAngle<=360;

cond_y1 = ProjectionAngle>=045 && ProjectionAngle<=135;
cond_y2 = ProjectionAngle>=225 && ProjectionAngle<=315;

% -- sampling length.
delta = 0.01;
% Note: fine/coarse sampling changes the results for intersection points
% and hence changes the list of cells intersected. TODO: Find the optimal
% value for delta.

if cond_x1 || cond_x3    
    
    % intermediate "delta xy points" from robot position up to the sensing
    % range in the direction of optical ray.
    deltaX = [0:delta:(x_endLimitPoint-x_start_cell),x_endLimitPoint-x_start_cell];
    deltaY = deltaX*tand(ProjectionAngle);    
    
elseif cond_x2    
    
    % intermediate "delta xy points" from robot position up to the sensing
    % range in the direction of optical ray.
    deltaX = [0:-delta:(x_endLimitPoint-x_start_cell),x_endLimitPoint-x_start_cell];
    deltaY = deltaX*tand(ProjectionAngle);
    
elseif cond_y1
    
    % intermediate "delta xy points" from robot position up to the sensing
    % range in the direction of optical ray.
    deltaY = [0:delta:(y_endLimitPoint-y_start_cell),y_endLimitPoint-y_start_cell];
    deltaX = deltaY/tand(ProjectionAngle);    
    
elseif cond_y2
    
    % intermediate "delta xy points" from robot position up to the sensing
    % range in the direction of optical ray.
    deltaY = [0:-delta:(y_endLimitPoint-y_start_cell),y_endLimitPoint-y_start_cell];
    deltaX = deltaY/tand(ProjectionAngle);
    
end

% intermediate xy points corresponding to intermediate delta xy points.
x = x_start_cell+deltaX;
y = y_start_cell+deltaY;

% 4 decimal precision
x = round(x*1e+4)/1e+4;
y = round(y*1e+4)/1e+4;

% list of xy points.
xy_points = [x',y'];

% xy points that belong to the intersection point (x.5,y.5)
ind_x  = find(abs(fix(x)-x)==0.5);
ind_yy = find(abs(fix(y(ind_x))-y(ind_x))==0.5);
ind_y  = ind_x(ind_yy);

for i = 1:numel(ind_y)
    xy_points = [xy_points;...
                 ceil(x(ind_y(i))),  ceil(y(ind_y(i)))  ;...
                 ceil(x(ind_y(i))),  floor(y(ind_y(i))) ;...
                 floor(x(ind_y(i))), ceil(y(ind_y(i)))  ;...
                 floor(x(ind_y(i))), floor(y(ind_y(i)))];
end

% -- convert xy_list from meter to cell numbers
% xy_cell = round(xy_points/dataYAML.resolution);
xy_cell = round(xy_points);


% convert xy_cell w.r.t robot origin in the map to cell numbers in the map
% xy_cell_map(:,1) = xy_cell(:,1)+robot_origin_in_map_x;
% xy_cell_map(:,2) = robot_origin_in_map_y-xy_cell(:,2);
% % xy_cell_map(:,2) = xy_cell(:,2)+robot_origin_in_map_y;

xy_cell_map(:,1) = xy_cell(:,1);
xy_cell_map(:,2) = xy_cell(:,2);

%---- truncated beam cells for the map size only
%=================================================


xy_cell_map = xy_cell_map(xy_cell_map(:,1)>0,:);
xy_cell_map = xy_cell_map(xy_cell_map(:,2)>0,:);
xy_cell_map = xy_cell_map(xy_cell_map(:,1)<=size(map_raycast,1),:);
xy_cell_map = xy_cell_map(xy_cell_map(:,2)<=size(map_raycast,2),:);

% plot(xy_cell_map(:,1),xy_cell_map(:,2),'xr')
% ------------------------------------------------------------------------------
% CHECKING FOR OBSTACLES TO FIND THE END POINT
% ------------------------------------------------------------------------------

% -- initialize end point cell number
x_end_cell = xy_cell_map(1,1);
y_end_cell = xy_cell_map(1,2);

% -- iterate untill you find an obstacle for sensing range end point.
% Note the row and col for the map.
for i = 1:size(xy_cell_map,1)
    xy_cell_map(i,:);
    map_raycast(xy_cell_map(i,1),xy_cell_map(i,2));
    if map_raycast(xy_cell_map(i,1),xy_cell_map(i,2))
        x_end_cell = xy_cell_map(i,1);
        y_end_cell = xy_cell_map(i,2);
    else
        break;
    end
end

% plot([x_start_cell,x_end_cell],[y_start_cell,y_end_cell],'-r')

% -- convert cell number back to meters from the robot origin in the map.
% x_end = (x_end_cell-robot_origin_in_map_x)*dataYAML.resolution
% y_end = (robot_origin_in_map_y-y_end_cell)*dataYAML.resolution
% % y_end = (y_end_cell-robot_origin_in_map_y)*PixSize;

% x_end_m = (x_end_cell-x_origin_cell)*dataYAML.resolution;
% y_end_m = (y_end_cell-y_origin_cell)*dataYAML.resolution;

x_end_m = (x_end_cell-x_origin_cell)*cellsize_raycast;
y_end_m = (y_end_cell-y_origin_cell)*cellsize_raycast;



end

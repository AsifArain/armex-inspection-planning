
function eSensorDom = fEffectiveSensorDomain( sensor,obs,SensorRange )
%
%
%

%-- PARAMETERS
%================================
visualize = 0;
SensorRange.cell = SensorRange.cell+1;
aperture_growth = 2/(distancePoints(sensor,obs)/SensorRange.cell);
beams_per_cell = 4;


% -- OBS ORIENTATIONS
%================================
% obs_angle = rad2deg(angle2Points(sensor,obs));

%-- obs corners
obs_corners = [ obs(1)+0.5,obs(2)+0.5;...
                obs(1)+0.5,obs(2)-0.5;...
                obs(1)-0.5,obs(2)+0.5;...
                obs(1)-0.5,obs(2)-0.5 ];

%-- angles to obs corners
obs_corners_angles = rad2deg(angle2Points(sensor,obs_corners));

%-- angles range (min, max)
obs_minmax_angles = [min(obs_corners_angles),max(obs_corners_angles)];

%-- aperture angle
aperture_angle = rad2deg(angleAbsDiff(deg2rad(obs_minmax_angles(1)),...
                                      deg2rad(obs_minmax_angles(2)) ));
                                  
%-- aperture circumference
aperture_cicum = 2*pi*SensorRange.cell*((2*aperture_growth*aperture_angle)/360);

%-- non-geometric aperture circumference (think about 350 deg to 10 deg)
aperture_cicum_non_geometric = abs(diff(obs_minmax_angles));


%-- to deal with special cases
if round(aperture_cicum_non_geometric*1e3)/1e3 == round(aperture_angle*1e3)/1e3

    obs_angles = linspace(obs_minmax_angles(1)-aperture_growth,...
                          obs_minmax_angles(2)+aperture_growth,...
                          ceil(aperture_cicum*beams_per_cell));
    
else
    
    obs_angles = linspace(obs_minmax_angles(2)-aperture_growth,...
                          obs_minmax_angles(2)+aperture_angle+aperture_growth,...
                          ceil(aperture_cicum*beams_per_cell));
    
end
   
%-- 360 deg
obs_angles = wrapTo360(obs_angles);


%-- beam end point
beam_end_x = sensor(1) + SensorRange.cell*cosd(obs_angles);
beam_end_y = sensor(2) + SensorRange.cell*sind(obs_angles);



% -- visulization for troubleshooting
if visualize == 1
    figure; hold on; 
    % vis_min_x = round(min([beam_initial_x,beam_end_____x]))-1
    % vis_max_x = round(max([beam_initial_x,beam_end_____x]))+1
    % 
    % vis_min_y = round(min([beam_initial_y,beam_end_____y]))-1
    % vis_max_y = round(max([beam_initial_y,beam_end_____y]))+1
    plot([beam_initial_x,beam_end_x],[beam_initial_y,beam_end_y],'-k');
    % set(gca,'xtick',[-round(beam_____range)-2:1:round(beam_____range)+2]);
    % set(gca,'ytick',[-round(beam_____range)-2:1:round(beam_____range)+2]);
    set(gca,'xtick',[round(beam_initial_x-beam_____range)-2:1:round(beam_initial_x+beam_____range)+2]);
    set(gca,'ytick',[round(beam_initial_y-beam_____range)-2:1:round(beam_initial_y+beam_____range)+2]);
    grid on; axis equal; 
    % xlim([-round(beam_____range)-2,round(beam_____range)+2]); 
    % ylim([-round(beam_____range)-2,round(beam_____range)+2]); 
    xlim([round(beam_initial_x-beam_____range)-2,round(beam_initial_x+beam_____range)+2]); 
    ylim([round(beam_initial_y-beam_____range)-2,round(beam_initial_y+beam_____range)+2]); 
end

% -- beam length along x-axis and y-axis
beam_dx = SensorRange.cell*cosd(obs_angles);
beam_dy = SensorRange.cell*sind(obs_angles);



%--- HERE WE GO
%===========================================


eSensorDom = [];

for num = 1:numel(obs_angles)



angle_this = obs_angles(num);
beam_dx_this = beam_dx(num);
beam_dy_this = beam_dy(num);
    
% -- conditions to project the line along (+/-x,+/-y).
cond_quad1 = angle_this > 000 && angle_this < 090;
cond_quad2 = angle_this > 090 && angle_this < 180;
cond_quad3 = angle_this > 180 && angle_this < 270;
cond_quad4 = angle_this > 270 && angle_this < 360;
cond____90 = angle_this == 090;
cond___270 = angle_this == 270;
cond_____0 = angle_this == 000 || angle_this == 360;
cond___180 = angle_this == 180;


% -- modeling optical beam on a regular grid
if cond_quad1
    
    %disp('--- at quadrant 1')
    
    % -- steps along x-axis
    dx1 = [0,ceil(sensor(1))-sensor(1):1:beam_dx_this,beam_dx_this];
    dy1 = dx1*tand(angle_this);
    
    % -- steps along y-axis
    dy2 = ceil(sensor(2))-sensor(2):1:beam_dy_this;
    dx2 = dy2/tand(angle_this);
    
    % -- all steps
    dx = ([dx1,dx2]);
    dy = ([dy1,dy2]);
    
    % -- <x,y> points
    x = sensor(1)+dx;
    y = sensor(2)+dy;
    
    % -- sorted list
    xy = [x',y'];
    %xy = unique(xy,'rows');
    xy = sortrows(xy,1);
    
   
elseif cond_quad2
    
    %disp('--- at quadrant 2')
    
    % -- steps along x-axis
    dx1 = [0,floor(sensor(1))-sensor(1):-1:beam_dx_this,beam_dx_this];
    dy1 = dx1*tand(angle_this);
    
    % -- steps along y-axis
    dy2 = ceil(sensor(2))-sensor(2):1:beam_dy_this;
    dx2 = dy2/tand(angle_this);
    
    % -- all steps
    dx = ([dx1,dx2]);
    dy = ([dy1,dy2]);
    
    % -- <x,y> points
    x = sensor(1)+dx;
    y = sensor(2)+dy;
    
    % -- sorted list
    xy = [x',y'];
    %xy = unique(xy,'rows');
    xy = sortrows(xy,-1);
    
    
elseif cond_quad3
    
    %disp('--- at quadrant 3')
    
    % -- steps along x-axis
    dx1 = [0,floor(sensor(1))-sensor(1):-1:beam_dx_this,beam_dx_this];
    dy1 = dx1*tand(angle_this);
    
    % -- steps along y-axis
    dy2 = floor(sensor(2))-sensor(2):-1:beam_dy_this;
    dx2 = dy2/tand(angle_this);
    
    % -- all steps
    dx = ([dx1,dx2]);
    dy = ([dy1,dy2]);
    
    % -- <x,y> points
    x = sensor(1)+dx;
    y = sensor(2)+dy;
    
    % -- sorted list
    xy = [x',y'];
    %xy = unique(xy,'rows');
    xy = sortrows(xy,-1);
    
    
    
elseif cond_quad4
    
    %disp('--- at quadrant 4')
    
    % -- steps along x-axis
    dx1 = [0,ceil(sensor(1))-sensor(1):1:beam_dx_this,beam_dx_this];
    dy1 = dx1*tand(angle_this);
    
    % -- steps along y-axis
    dy2 = floor(sensor(2))-sensor(2):-1:beam_dy_this;
    dx2 = dy2/tand(angle_this);
   
    % -- all steps
    dx = ([dx1,dx2]);
    dy = ([dy1,dy2]);
    
    % -- <x,y> points
    x = sensor(1)+dx;
    y = sensor(2)+dy;
    
    % -- sorted list
    xy = [x',y'];
    %xy = unique(xy,'rows');
    xy = sortrows(xy,1);
    
    
elseif cond____90
    
    %disp('--- at 90 deg')
    
    % -- steps along x-axis
    dx1 = []; 
    dy1 = []; 
    
    % -- steps along y-axis    
    dy2 = [0,ceil(sensor(2))-sensor(2):1:beam_dy_this,beam_dy_this];
    dx2 = dy2/tand(angle_this);
    
    % -- all steps
    dx = ([dx1,dx2]);
    dy = ([dy1,dy2]);
    
    % -- <x,y> points
    x = sensor(1)+dx;
    y = sensor(2)+dy;
    
    % -- sorted list
    xy = [x',y'];
    %xy = unique(xy,'rows');
    xy = sortrows(xy,2);
    
    
elseif cond___270
    
    %disp('--- at 270 deg')
    
    % -- steps along x-axis
    dx1 = []; 
    dy1 = []; 
    
    % -- steps along y-axis
    %dy2 = ceil(sensor(2))-sensor(2):1:beam________dy;
    dy2 = [0,floor(sensor(2))-sensor(2):-1:beam_dy_this,beam_dy_this];
    dx2 = dy2/tand(angle_this);
    
    % -- all steps
    dx = ([dx1,dx2]);
    dy = ([dy1,dy2]);
    
    % -- <x,y> points
    x = sensor(1)+dx;
    y = sensor(2)+dy;
    
    % -- sorted list
    xy = [x',y'];
    %xy = unique(xy,'rows');
    xy = sortrows(xy,-2);
    

elseif cond_____0
    
    %disp('--- at quadrant 1')
    
    % -- steps along x-axis
    dx1 = [0,ceil(sensor(1))-sensor(1):1:beam_dx_this,beam_dx_this];
    dy1 = dx1*tand(angle_this);
    
    % -- steps along y-axis
    dy2 = []; %ceil(sensor(2))-sensor(2):1:beam________dy;
    dx2 = []; %dy2/tand(beam_direction);
    
    % -- all steps
    dx = ([dx1,dx2]);
    dy = ([dy1,dy2]);
    
    % -- <x,y> points
    x = sensor(1)+dx;
    y = sensor(2)+dy;
    
    % -- sorted list
    xy = [x',y'];
    %xy = unique(xy,'rows');
    xy = sortrows(xy,1);
    
    
elseif cond___180
    
    %disp('--- at quadrant 2')
    
    % -- steps along x-axis
    dx1 = [0,floor(sensor(1))-sensor(1):-1:beam_dx_this,beam_dx_this];
    dy1 = dx1*tand(angle_this);
    
    % -- steps along y-axis
    dy2 = []; %ceil(sensor(2))-sensor(2):1:beam________dy;    
    dx2 = []; %dy2/tand(beam_direction);
    
    % -- all steps
    dx = ([dx1,dx2]);
    dy = ([dy1,dy2]);
    
    % -- <x,y> points
    x = sensor(1)+dx;
    y = sensor(2)+dy;
    
    % -- sorted list
    xy = [x',y'];
    %xy = unique(xy,'rows');
    xy = sortrows(xy,-1);
    
           
    
end




% -- visualize for troubleshooting
% if visualize == 1
%     plot(xy(:,1),xy(:,2),'or');
%     plot(xy(:,1),xy(:,2),'.r','MarkerSize',10);
% end

% -- voxels list
% find mid points along steps
d_midpoints = zeros(2,size(xy,1)-1);
for i = 1:size(d_midpoints,2)
    d_midpoints(1,i) = (xy(i+1,1)-xy(i,1))/2;
    d_midpoints(2,i) = (xy(i+1,2)-xy(i,2))/2;
end
% mid points for the voxels
beam_voxel_midpts = (xy(1:end-1,:))'+d_midpoints;

% -- visualize for troubleshooting
if visualize == 1
    plot(beam_voxel_midpts(1,:),beam_voxel_midpts(2,:),'xk');
    %pause
end

% -- voxels list and end points
beam_voxel_list = floor(beam_voxel_midpts);
% beam_voxel_points = (xy(2:end,:))';

% -- visualize for troubleshooting
if visualize == 1
    plot(beam_voxel_list(1,:)+0.5,beam_voxel_list(2,:)+0.5,'sg','MarkerSize',20);
    for i = 1:size(beam_voxel_list,2)
        plot(beam_voxel_list(1,i)+0.5,beam_voxel_list(2,i)+0.5,'ob','MarkerSize',25);
        %pause
    end
end

eSensorDom = [eSensorDom;beam_voxel_list'];


end

eSensorDom = unique(eSensorDom,'rows');

% -- beam length through voxels
% beam_voxel_length = zeros(1,size(xy,1)-1);
% for i = 1:size(beam_voxel_length,2)
%     beam_voxel_length(i) = sqrt((xy(i+1,2)-xy(i,2))^2+(xy(i+1,1)-xy(i,1))^2);
% end

    
        

end
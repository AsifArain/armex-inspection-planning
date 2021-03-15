% function [beam_end_____x,...
%           beam_end_____y,...
%           beam_voxel___list,...
%           beam_voxel_length] = fRayCasting(beam_initial_x,...
%                                            beam_initial_y,...
%                                            beam_direction,...
%                                            beam_____range)

close all; clear all; clc;


map = ones(15,15);
% map_num: 11
%map([48,49,64,65,71,79,80,85,86,94,100,101,114:116,129:131,145:147,161,169:171,184:186]) = 0;
% map_num: 12
%map([48:50,63:65,71,72,78:80,85:87,93:95,100:102,114:117,129:132,145:147,161,162,169:171,184:186]) = 0;
% map_num: 13
%map([34,48:50,63:65,71,72,78:80,85:87,93:95,99,100:102,109,114:117,129:132,145:147,161,162,169:171,184:187]) = 0;
map(10,3) = 0;
map(11,4) = 0;
[map_row,map_col] = size(map)

map_concentration = 1000*rand(size(map))

% -- beam initial point
beam_initial_x = 6.5
beam_initial_y = 6.5
% beam_direction = 16.7826 %cond_x1
% beam_direction = 163.2174 %cond_x2
% beam_direction = 360-16.7826 %cond_x3

% -- beam direction
beam_direction = 180 
% beam_direction = atan2((beam_end_____y-beam_initial_y),(beam_end_____x-beam_initial_x))*(180/pi);
beam_direction = beam_direction+(360*(beam_direction<0));

% -- add a random number between -1 to 1 if beam direction is 45, 135, 225, or 315 deg
beam_direction = beam_direction+((2*rand(1)-1)*(...
                                                beam_direction== 45||... %  45 deg
                                                beam_direction==135||... % 135 deg
                                                beam_direction==225||... % 225 deg
                                                beam_direction==315  ... % 315 deg
                                                ));

% -- beam range
beam_____range = 6.55

% -- beam end point
beam_end_____x = beam_initial_x + beam_____range*cosd(beam_direction);
beam_end_____y = beam_initial_y + beam_____range*sind(beam_direction);

% ---
figure; hold on; 

for i = 1:map_row
    for j = 1:map_col
        % if obstacle/free
        if ~map(i,j) % occupied
            %col = [0.4 0.4 0.4];
            %faceCol = [0.4 0.4 0.4];
            col = [0 0 0];
            faceCol = [0.4 0.4 0.4];
        elseif map(i,j) % free
            %col = [0.9 0.9 0.9];
            col = [0 0 0];
            faceCol = [1 1 1];
            %col = [1 1 1];
        end
        plot(i-0.5,j-0.5,'S','Color',col,'MarkerFaceColor',faceCol,'MarkerSize',25); hold on;
    end
end
axis equal

%pause
% -- visulization for troubleshooting
plot([beam_initial_x,beam_end_____x],[beam_initial_y,beam_end_____y],'-k');
% set(gca,'xtick',[-10:1:10]);
% set(gca,'ytick',[-10:1:10]);
% grid on; axis equal; 
% xlim([-8,8]); ylim([-8,8]); 

% -- beam length along x-axis and y-axis
beam________dx = beam_____range*cosd(beam_direction);
beam________dy = beam_____range*sind(beam_direction);


% -- conditions to project the line along (+/-x,+/-y).
cond_quad1 = beam_direction>=0  && beam_direction<90;
cond_quad2 = beam_direction>90  && beam_direction<=180;
cond_quad3 = beam_direction>180 && beam_direction<270;
cond_quad4 = beam_direction>270 && beam_direction<=360;
cond____90 = beam_direction==90;
cond___270 = beam_direction==270;

% -- modeling optical beam on a regular grid
if cond_quad1
    
    disp('--- at quadrant 1')
    
    % -- steps along x-axis
    dx1 = [0,ceil(beam_initial_x)-beam_initial_x:1:beam________dx,beam________dx];
    dy1 = dx1*tand(beam_direction);
    
    % -- steps along y-axis
    dy2 = ceil(beam_initial_y)-beam_initial_y:1:beam________dy;
    dx2 = dy2/tand(beam_direction);
    
    % -- all steps
    dx = ([dx1,dx2]);
    dy = ([dy1,dy2]);
    
    % -- <x,y> points
    x = beam_initial_x+dx;
    y = beam_initial_y+dy;
    
    % -- sorted list
    xy = [x',y'];
    xy = unique(xy,'rows');
    xy = sortrows(xy,1);
    
    % -- visualize for troubleshooting
    plot(xy(:,1),xy(:,2),'or');
    plot(xy(:,1),xy(:,2),'.r','MarkerSize',10);
    
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
    plot(beam_voxel_midpts(1,:),beam_voxel_midpts(2,:),'xk');
    
    % -- voxels list
    beam_voxel___list = ceil(beam_voxel_midpts);
    
    % -- visualize for troubleshooting
    plot(beam_voxel___list(1,:)-0.5,beam_voxel___list(2,:)-0.5,'sg','MarkerSize',20);
    for i = 1:size(beam_voxel___list,2)
        plot(beam_voxel___list(1,i)-0.5,beam_voxel___list(2,i)-0.5,'ob','MarkerSize',25);
        %pause
    end
    
    % -- beam length through voxels
    beam_voxel_length = zeros(1,size(xy,1)-1);
    %for i = 1:size(beam_voxel_length,2)
    %    beam_voxel_length(i) = sqrt((xy(i+1,2)-xy(i,2))^2+(xy(i+1,1)-xy(i,1))^2);
    %end
    
    
elseif cond_quad2
    
    disp('--- at quadrant 2')
    
    % -- steps along x-axis
    dx1 = [0,floor(beam_initial_x)-beam_initial_x:-1:beam________dx,beam________dx];
    dy1 = dx1*tand(beam_direction);
    
    % -- steps along y-axis
    dy2 = ceil(beam_initial_y)-beam_initial_y:1:beam________dy;
    dx2 = dy2/tand(beam_direction);
    
    % -- all steps
    dx = ([dx1,dx2]);
    dy = ([dy1,dy2]);
    
    % -- <x,y> points
    x = beam_initial_x+dx;
    y = beam_initial_y+dy;
    
    % -- sorted list
    xy = [x',y'];
    xy = unique(xy,'rows');
    xy = sortrows(xy,-1);
    
    % -- visualize for troubleshooting
    plot(xy(:,1),xy(:,2),'or');
    plot(xy(:,1),xy(:,2),'.r','MarkerSize',10);
    
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
    plot(beam_voxel_midpts(1,:),beam_voxel_midpts(2,:),'xk');
    
    % -- voxels list
    beam_voxel___list = ceil(beam_voxel_midpts);
        
    % -- visualize for troubleshooting
    plot(beam_voxel___list(1,:)-0.5,beam_voxel___list(2,:)-0.5,'sg','MarkerSize',20);
    %for i = 1:size(beam_voxel___list,2)
    %    plot(beam_voxel___list(1,i)-0.5,beam_voxel___list(2,i)-0.5,'ob','MarkerSize',25);
    %    pause
    %end
    
    % -- beam length through voxels
    beam_voxel_length = zeros(1,size(xy,1)-1);
    for i = 1:size(beam_voxel_length,2)
        beam_voxel_length(i) = sqrt((xy(i+1,2)-xy(i,2))^2+(xy(i+1,1)-xy(i,1))^2);
    end
        
    
elseif cond_quad3
    
    disp('--- at quadrant 3')
    
    % -- steps along x-axis
    dx1 = [0,floor(beam_initial_x)-beam_initial_x:-1:beam________dx,beam________dx];
    dy1 = dx1*tand(beam_direction);
    
    % -- steps along y-axis
    dy2 = floor(beam_initial_y)-beam_initial_y:-1:beam________dy;
    dx2 = dy2/tand(beam_direction);
    
    % -- all steps
    dx = ([dx1,dx2]);
    dy = ([dy1,dy2]);
    
    % -- <x,y> points
    x = beam_initial_x+dx;
    y = beam_initial_y+dy;
    
    % -- sorted list
    xy = [x',y'];
    xy = unique(xy,'rows');
    xy = sortrows(xy,-1);
    
    % -- visualize for troubleshooting
    plot(xy(:,1),xy(:,2),'or');
    plot(xy(:,1),xy(:,2),'.r','MarkerSize',10);
    
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
    plot(beam_voxel_midpts(1,:),beam_voxel_midpts(2,:),'xk');
    
    % -- voxels list
    beam_voxel___list = ceil(beam_voxel_midpts);
    
    % -- visualize for troubleshooting
    plot(beam_voxel___list(1,:)-0.5,beam_voxel___list(2,:)-0.5,'sg','MarkerSize',20);
    %for i = 1:size(beam_voxel___list,2)
    %    plot(beam_voxel___list(1,i)-0.5,beam_voxel___list(2,i)-0.5,'ob','MarkerSize',25);
    %    pause
    %end
    
    % -- beam length through voxels
    beam_voxel_length = zeros(1,size(xy,1)-1);
    for i = 1:size(beam_voxel_length,2)
        beam_voxel_length(i) = sqrt((xy(i+1,2)-xy(i,2))^2+(xy(i+1,1)-xy(i,1))^2);
    end
    
    
elseif cond_quad4
    
    disp('--- at quadrant 4')
    
    % -- steps along x-axis
    dx1 = [0,ceil(beam_initial_x)-beam_initial_x:1:beam________dx,beam________dx];
    dy1 = dx1*tand(beam_direction);
    
    % -- steps along y-axis
    dy2 = floor(beam_initial_y)-beam_initial_y:-1:beam________dy;
    dx2 = dy2/tand(beam_direction);
   
    % -- all steps
    dx = ([dx1,dx2]);
    dy = ([dy1,dy2]);
    
    % -- <x,y> points
    x = beam_initial_x+dx;
    y = beam_initial_y+dy;
    
    % -- sorted list
    xy = [x',y'];
    xy = unique(xy,'rows');
    xy = sortrows(xy,1);
    
    % -- visualize for troubleshooting
    plot(xy(:,1),xy(:,2),'or');
    plot(xy(:,1),xy(:,2),'.r','MarkerSize',10);
    
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
    plot(beam_voxel_midpts(1,:),beam_voxel_midpts(2,:),'xk');
    
    % -- voxels list
    beam_voxel___list = ceil(beam_voxel_midpts);
        
    % -- visualize for troubleshooting
    plot(beam_voxel___list(1,:)-0.5,beam_voxel___list(2,:)-0.5,'sg','MarkerSize',20);
    %for i = 1:size(beam_voxel___list,2)
    %    plot(beam_voxel___list(1,i)-0.5,beam_voxel___list(2,i)-0.5,'ob','MarkerSize',25);
    %    %pause
    %end
    
    % -- beam length through voxels
    beam_voxel_length = zeros(1,size(xy,1)-1);
    for i = 1:size(beam_voxel_length,2)
        beam_voxel_length(i) = sqrt((xy(i+1,2)-xy(i,2))^2+(xy(i+1,1)-xy(i,1))^2);
    end
    
    
elseif cond____90
    
    disp('--- at 90 deg')
    
    % -- steps along x-axis
    dx1 = []; %[0,ceil(beam_initial_x)-beam_initial_x:1:beam________dx,beam________dx];
    dy1 = []; %dx1*tand(beam_direction);
    
    % -- steps along y-axis
    %dy2 = ceil(beam_initial_y)-beam_initial_y:1:beam________dy;
    dy2 = [0,ceil(beam_initial_y)-beam_initial_y:1:beam________dy,beam________dy];
    dx2 = dy2/tand(beam_direction);
    
    % -- all steps
    dx = ([dx1,dx2]);
    dy = ([dy1,dy2]);
    
    % -- <x,y> points
    x = beam_initial_x+dx;
    y = beam_initial_y+dy;
    
    % -- sorted list
    xy = [x',y'];
    xy = unique(xy,'rows');
    xy = sortrows(xy,2);
    
    % -- visualize for troubleshooting
    plot(xy(:,1),xy(:,2),'or');
    plot(xy(:,1),xy(:,2),'.r','MarkerSize',10);
    
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
    plot(beam_voxel_midpts(1,:),beam_voxel_midpts(2,:),'xk');
    
    % -- voxels list
    beam_voxel___list = ceil(beam_voxel_midpts);
    
    % -- visualize for troubleshooting
    plot(beam_voxel___list(1,:)-0.5,beam_voxel___list(2,:)-0.5,'sg','MarkerSize',20);
    %for i = 1:size(beam_voxel___list,2)
    %    plot(beam_voxel___list(1,i)-0.5,beam_voxel___list(2,i)-0.5,'ob','MarkerSize',25);
    %    %pause
    %end
    
    % -- beam length through voxels
    beam_voxel_length = zeros(1,size(xy,1)-1);
    for i = 1:size(beam_voxel_length,2)
        beam_voxel_length(i) = sqrt((xy(i+1,2)-xy(i,2))^2+(xy(i+1,1)-xy(i,1))^2);
    end
    
elseif cond___270
    
    disp('--- at 270 deg')
    
    % -- steps along x-axis
    dx1 = []; %[0,ceil(beam_initial_x)-beam_initial_x:1:beam________dx,beam________dx];
    dy1 = []; %dx1*tand(beam_direction);
    
    % -- steps along y-axis
    %dy2 = ceil(beam_initial_y)-beam_initial_y:1:beam________dy;
    dy2 = [0,floor(beam_initial_y)-beam_initial_y:-1:beam________dy,beam________dy];
    dx2 = dy2/tand(beam_direction);
    
    % -- all steps
    dx = ([dx1,dx2]);
    dy = ([dy1,dy2]);
    
    % -- <x,y> points
    x = beam_initial_x+dx;
    y = beam_initial_y+dy;
    
    % -- sorted list
    xy = [x',y'];
    xy = unique(xy,'rows');
    xy = sortrows(xy,-2);
    
    % -- visualize for troubleshooting
    plot(xy(:,1),xy(:,2),'or');
    plot(xy(:,1),xy(:,2),'.r','MarkerSize',10);
    
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
    plot(beam_voxel_midpts(1,:),beam_voxel_midpts(2,:),'xk');
    
    % -- voxels list
    beam_voxel___list = ceil(beam_voxel_midpts);
    
    % -- visualize for troubleshooting
    plot(beam_voxel___list(1,:)-0.5,beam_voxel___list(2,:)-0.5,'sg','MarkerSize',20);
    %for i = 1:size(beam_voxel___list,2)
    %    plot(beam_voxel___list(1,i)-0.5,beam_voxel___list(2,i)-0.5,'ob','MarkerSize',25);
    %    %pause
    %end
    
    % -- beam length through voxels
    beam_voxel_length = zeros(1,size(xy,1)-1);
    for i = 1:size(beam_voxel_length,2)
        beam_voxel_length(i) = sqrt((xy(i+1,2)-xy(i,2))^2+(xy(i+1,1)-xy(i,1))^2);
    end    
    
end

% -- measurement of integral concentration on grid map with obstacles
ppm_m = 0;
for i = 1:size(beam_voxel___list,2)
    
    % if voxel is out of map
    if beam_voxel___list(1,i) <= 0 || beam_voxel___list(1,i) > map_row ||...
       beam_voxel___list(2,i) <= 0 || beam_voxel___list(2,i) > map_col
       
        break;
    end
    
    % if voxel is an occupied cell
    if map(beam_voxel___list(1,i),beam_voxel___list(2,i)) == 0
        
        break;
    end
        
    this_ppm_m = beam_voxel_length(i)*map_concentration(beam_voxel___list(1,i),beam_voxel___list(2,i));
    
    ppm_m = ppm_m + this_ppm_m;
    
    %ppm_m = ppm_m + beam_voxel_length(i)*...
    %                map_concentration(beam_voxel___list(1,i),beam_voxel___list(2,i));
    
end
size(beam_voxel___list,2)
i
ppm_m


% % remove repeated cells
% cell_xy = unique(cell_xy,'rows');
% 
% % validation for the visibility
% ind  = find(cell_xy(:,1)==obs(1)&cell_xy(:,2)==obs(2));
% % if size(ind,1) > 0
% %     validation = 0;
% % end
% validation = (size(ind,1) < 1);


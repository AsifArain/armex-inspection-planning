function [beam_end_____x,...
          beam_end_____y,...
          beam_voxel___list,...
          beam_voxel_points,...
          beam_voxel_midpts,...
          beam_voxel_length] = fBeamOnSimpleGridA(beam_initial_x,...
                                                  beam_initial_y,...
                                                  beam_direction,...
                                                  beam_____range,...
                                                  map___environment)
%

visualize = 0;
% -- beam initial point
%beam_initial_x = 6.5
%beam_initial_y = 6.5

% -- beam direction
%beam_direction = 180 
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
%beam_____range = 6.55

% -- beam end point
beam_end_____x = beam_initial_x + beam_____range*cosd(beam_direction);
beam_end_____y = beam_initial_y + beam_____range*sind(beam_direction);

% ---
% fsamina
% figure; hold on; 
% for i = 1:map_row
%     for j = 1:map_col
%         % if obstacle/free
%         if ~map(i,j) % occupied
%             %col = [0.4 0.4 0.4];
%             %faceCol = [0.4 0.4 0.4];
%             col = [0 0 0];
%             faceCol = [0.4 0.4 0.4];
%         elseif map(i,j) % free
%             %col = [0.9 0.9 0.9];
%             col = [0 0 0];
%             faceCol = [1 1 1];
%             %col = [1 1 1];
%         end
%         plot(i-0.5,j-0.5,'S','Color',col,'MarkerFaceColor',faceCol,'MarkerSize',25); hold on;
%     end
% end
% axis equal

% -- visulization for troubleshooting
if visualize == 1
    figure; hold on; 
    % vis_min_x = round(min([beam_initial_x,beam_end_____x]))-1
    % vis_max_x = round(max([beam_initial_x,beam_end_____x]))+1
    % 
    % vis_min_y = round(min([beam_initial_y,beam_end_____y]))-1
    % vis_max_y = round(max([beam_initial_y,beam_end_____y]))+1
    plot([beam_initial_x,beam_end_____x],[beam_initial_y,beam_end_____y],'-k');
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
beam________dx = beam_____range*cosd(beam_direction);
beam________dy = beam_____range*sind(beam_direction);


% -- conditions to project the line along (+/-x,+/-y).
cond_quad1 = beam_direction>=000 && beam_direction< 090;
cond_quad2 = beam_direction> 090 && beam_direction<=180;
cond_quad3 = beam_direction> 180 && beam_direction< 270;
cond_quad4 = beam_direction> 270 && beam_direction<=360;
cond____90 = beam_direction==090;
cond___270 = beam_direction==270;


% -- modeling optical beam on a regular grid
if cond_quad1
    
    %disp('--- at quadrant 1')
    
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
    if visualize == 1
        plot(xy(:,1),xy(:,2),'or');
        plot(xy(:,1),xy(:,2),'.r','MarkerSize',10);
    end
    
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
    beam_voxel___list = floor(beam_voxel_midpts);
    beam_voxel_points = (xy(2:end,:))';
    
    % -- visualize for troubleshooting
    if visualize == 1
        plot(beam_voxel___list(1,:)+0.5,beam_voxel___list(2,:)+0.5,'sg','MarkerSize',20);
        for i = 1:size(beam_voxel___list,2)
            plot(beam_voxel___list(1,i)+0.5,beam_voxel___list(2,i)+0.5,'ob','MarkerSize',25);
            %pause
        end
    end
    
    % -- beam length through voxels
    beam_voxel_length = zeros(1,size(xy,1)-1);
    for i = 1:size(beam_voxel_length,2)
        beam_voxel_length(i) = sqrt((xy(i+1,2)-xy(i,2))^2+(xy(i+1,1)-xy(i,1))^2);
    end
    
    
elseif cond_quad2
    
    %disp('--- at quadrant 2')
    
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
    if visualize == 1
        plot(xy(:,1),xy(:,2),'or');
        plot(xy(:,1),xy(:,2),'.r','MarkerSize',10);
        %pause
    end
    
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
    beam_voxel___list = floor(beam_voxel_midpts);
    beam_voxel_points = (xy(2:end,:))';
        
    % -- visualize for troubleshooting
    if visualize == 1
        plot(beam_voxel___list(1,:)+0.5,beam_voxel___list(2,:)+0.5,'sg','MarkerSize',20);
        for i = 1:size(beam_voxel___list,2)
            plot(beam_voxel___list(1,i)+0.5,beam_voxel___list(2,i)+0.5,'ob','MarkerSize',25);
            %pause
        end
    end
    
    % -- beam length through voxels
    beam_voxel_length = zeros(1,size(xy,1)-1);
    for i = 1:size(beam_voxel_length,2)
        beam_voxel_length(i) = sqrt((xy(i+1,2)-xy(i,2))^2+(xy(i+1,1)-xy(i,1))^2);
    end
        
    
elseif cond_quad3
    
    %disp('--- at quadrant 3')
    
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
    if visualize == 1
        plot(xy(:,1),xy(:,2),'or');
        plot(xy(:,1),xy(:,2),'.r','MarkerSize',10);
        %pause
    end
    
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
    beam_voxel___list = floor(beam_voxel_midpts);
    beam_voxel_points = (xy(2:end,:))';
    
    % -- visualize for troubleshooting
    if visualize == 1
        plot(beam_voxel___list(1,:)+0.5,beam_voxel___list(2,:)+0.5,'sg','MarkerSize',20);
        for i = 1:size(beam_voxel___list,2)
            plot(beam_voxel___list(1,i)+0.5,beam_voxel___list(2,i)+0.5,'ob','MarkerSize',25);
            %pause
        end
    end
    
    % -- beam length through voxels
    beam_voxel_length = zeros(1,size(xy,1)-1);
    for i = 1:size(beam_voxel_length,2)
        beam_voxel_length(i) = sqrt((xy(i+1,2)-xy(i,2))^2+(xy(i+1,1)-xy(i,1))^2);
    end
    
    
elseif cond_quad4
    
    %disp('--- at quadrant 4')
    
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
    if visualize == 1
        plot(xy(:,1),xy(:,2),'or');
        plot(xy(:,1),xy(:,2),'.r','MarkerSize',10);
        %pause
    end
    
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
    if visualize
        plot(beam_voxel_midpts(1,:),beam_voxel_midpts(2,:),'xk');
    end
    
    % -- voxels list and end points
    beam_voxel___list = floor(beam_voxel_midpts);
    beam_voxel_points = (xy(2:end,:))';
        
    % -- visualize for troubleshooting
    if visualize == 1
        plot(beam_voxel___list(1,:)+0.5,beam_voxel___list(2,:)+0.5,'sg','MarkerSize',20);
        for i = 1:size(beam_voxel___list,2)
            plot(beam_voxel___list(1,i)+0.5,beam_voxel___list(2,i)+0.5,'ob','MarkerSize',25);
            %pause
        end
    end
    
    % -- beam length through voxels
    beam_voxel_length = zeros(1,size(xy,1)-1);
    for i = 1:size(beam_voxel_length,2)
        beam_voxel_length(i) = sqrt((xy(i+1,2)-xy(i,2))^2+(xy(i+1,1)-xy(i,1))^2);
    end
    
    
elseif cond____90
    
    %disp('--- at 90 deg')
    
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
    if visualize == 1
        plot(xy(:,1),xy(:,2),'or');
        plot(xy(:,1),xy(:,2),'.r','MarkerSize',10);
    end
    
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
    end
    
       
    % -- voxels list and end points
    beam_voxel___list = floor(beam_voxel_midpts);
    beam_voxel_points = (xy(2:end,:))';
    
    % -- visualize for troubleshooting
    if visualize == 1
        plot(beam_voxel___list(1,:)+0.5,beam_voxel___list(2,:)+0.5,'sg','MarkerSize',20);
        for i = 1:size(beam_voxel___list,2)
            plot(beam_voxel___list(1,i)+0.5,beam_voxel___list(2,i)+0.5,'ob','MarkerSize',25);
            %pause
        end
    end
    
    % -- beam length through voxels
    beam_voxel_length = zeros(1,size(xy,1)-1);
    for i = 1:size(beam_voxel_length,2)
        beam_voxel_length(i) = sqrt((xy(i+1,2)-xy(i,2))^2+(xy(i+1,1)-xy(i,1))^2);
    end
    
elseif cond___270
    
    %disp('--- at 270 deg')
    
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
    if visualize == 1
        plot(xy(:,1),xy(:,2),'or');
        plot(xy(:,1),xy(:,2),'.r','MarkerSize',10);
    end
    
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
    end

    
    % -- voxels list and end points
    beam_voxel___list = floor(beam_voxel_midpts);
    beam_voxel_points = (xy(2:end,:))';
    
    % -- visualize for troubleshooting
    if visualize == 1
        plot(beam_voxel___list(1,:)+0.5,beam_voxel___list(2,:)+0.5,'sg','MarkerSize',20);
        for i = 1:size(beam_voxel___list,2)
            plot(beam_voxel___list(1,i)+0.5,beam_voxel___list(2,i)+0.5,'ob','MarkerSize',25);
            %pause
        end
    end
    
    % -- beam length through voxels
    beam_voxel_length = zeros(1,size(xy,1)-1);
    for i = 1:size(beam_voxel_length,2)
        beam_voxel_length(i) = sqrt((xy(i+1,2)-xy(i,2))^2+(xy(i+1,1)-xy(i,1))^2);
    end    
    
end



end
clear;
clc;
close all;

% map_spec_name='MapSpec_201510211300';
% occ_map_loc='fr-campus-20040714.carmen.gfs.png';
occ_map_loc='freiburg.png';
% obstacles_loc=load_occupancy_map_ASIF('figure',occ_map_loc,'spec_file',map_spec_name,'plot_map',0);


offset=[0 0];
resolution=0.05;
rotation=0;
markercolor='k';
markersize=1.5;
only_data=1;




[I,map_color]=imread(occ_map_loc);
I = im2bw(I,map_color,0.29); % binary map
% I = imbinarize(I,map_color,0.95); % binary map


% [X,map_color] = imread(paramPreprocess.MapFig);
% Binary/BW map.
%map_BW = im2bw(X,map_color,0.95); % binary map
% map_BW = X;
        
% I_gray=rgb2gray(I)'; figure; image(I_gray)
% map_input=double(I_gray);
map_input=double(I); figure; image(I);
offset_x=offset(1);
offset_y=offset(2);
[rows,cols]=find(map_input==0);
% obstacle_mat=[cols(:),rows(:)];
obstacle_mat=[rows(:),cols(:)];
obstacle_mat=obstacle_mat.*resolution;
obstacle_mat(:,1)=obstacle_mat(:,1)+offset_x;
obstacle_mat(:,2)=obstacle_mat(:,2)+offset_y;
RA=[cosd(rotation) -sind(rotation); sind(rotation) cosd(rotation)];
obstacle_mat=obstacle_mat*RA;

if only_data==1
    plot(obstacle_mat(:,1),obstacle_mat(:,2),...
        'sk','markerfacecolor',markercolor,'markersize',markersize);
end


obstacles_loc = obstacle_mat;

figure;plot(obstacles_loc(:,2),obstacles_loc(:,1),'bx')
hold on;
plot(0,0,'ro','markerfacecolor','r','markersize',10)


%/////////////////////////////////////////////////////////////////////////7

obstacles_loc2 = unique(round(obstacles_loc),'rows');

figure;plot(obstacles_loc2(:,1),obstacles_loc2(:,2),'bs','MarkerFaceColor','b')
axis equal


obstacles_loc2(:,1) = obstacles_loc2(:,1)-min(obstacles_loc2(:,1))+1;
obstacles_loc2(:,2) = obstacles_loc2(:,2)-min(obstacles_loc2(:,2))+1;

figure;plot(obstacles_loc2(:,1),obstacles_loc2(:,2),'bs','MarkerFaceColor','b')
axis equal



map = ones(max(obstacles_loc2(:,2)),max(obstacles_loc2(:,1)));

obs_ind = sub2ind(size(map),obstacles_loc2(:,2),obstacles_loc2(:,1));
map(obs_ind) = 0;


figure; imshow(map','InitialMagnification',500); hold on;
set(gca,'YDir','normal');
colormap('gray'); 

filename = 'freiburg.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,map,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);


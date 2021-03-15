clear;
clc;
close all;

% map_spec_name='MapSpec_201510211300';
occ_map_loc='sample_environment.png';
% obstacles_loc=load_occupancy_map_ASIF('figure',occ_map_loc,'spec_file',map_spec_name,'plot_map',0);

offset=[0 0];
resolution=0.07;
rotation=0;
markercolor='k';
markersize=1.5;
only_data=1;




I=imread(occ_map_loc);
I_gray=rgb2gray(I);

% I_gray = imrotate(I_gray,180);
% I_gray = flipdim(I_gray,2);
I_gray = flipdim(I_gray,1);

map_input=double(I_gray);

% map_input=double(I);
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
axis equal



%/////////////////////////////////////////////////////////////////////////7

% obstacles_loc2 = unique(round(obstacles_loc),'rows');
obstacles_loc2 = unique(ceil(obstacles_loc),'rows');

figure;plot(obstacles_loc2(:,1),obstacles_loc2(:,2),'bs','MarkerFaceColor','b')
axis equal

% pause

obstacles_loc2(:,1) = obstacles_loc2(:,1)-min(obstacles_loc2(:,1))+1;
obstacles_loc2(:,2) = obstacles_loc2(:,2)-min(obstacles_loc2(:,2))+1;

figure;plot(obstacles_loc2(:,1),obstacles_loc2(:,2),'bs','MarkerFaceColor','b')
axis equal



map = ones(max(obstacles_loc2(:,2)),max(obstacles_loc2(:,1)));

obs_ind = sub2ind(size(map),obstacles_loc2(:,2),obstacles_loc2(:,1));
map(obs_ind) = 0;


% -- fill holes
map = 1-map;
map_remove = imfill(map,4,'holes');
figure; imshow(map','InitialMagnification',500); hold on;
ind_remove = find(map_remove(:)==0);
% map = bwmorph(map,'fill',5000);
% map = imfill(map,0,4);
map = 1-map;
figure; imshow(map','InitialMagnification',500); hold on;

map(ind_remove) = 0;

figure; imshow(map','InitialMagnification',500); hold on;
% set(gca,'YDir','normal');
colormap('gray'); 

% figure; imshow(map','InitialMagnification',500); hold on;
% set(gca,'YDir','normal');
% colormap('gray'); 


% -------------- fill holes ------------------
% map = 1-map;
% % https://se.mathworks.com/matlabcentral/newsreader/view_thread/173142
% map = ~bwareaopen(~map, 100);
% map = 1-map;


% ---------- fill holes using 4-connectivity ------------------------
CC = bwconncomp(map,4);
conn_sizes = zeros(CC.NumObjects,1);
for i = 1:CC.NumObjects
    conn_sizes(i) = size(CC.PixelIdxList{i},1);
end
[val,ind] = max(conn_sizes);
map = zeros(size(map));
map(CC.PixelIdxList{ind}) = 1;

% ------------ morphological operations ---------------------------
% map = bwmorph(map,'remove');
% map = bwmorph(map,'fill',5000);
% map = bwmorph(map,'majority',5);

close all

figure; imshow(map','InitialMagnification',500); hold on;
set(gca,'YDir','normal');
colormap('gray'); 


filename = 'sample-environment.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,map,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);


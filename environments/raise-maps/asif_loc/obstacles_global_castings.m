clear;
clc;
close all;


occ_offset=[-66.5 -7.5];
occ_rot=0;
occ_resolution=0.07;



occ_map_loc='GHG_global2.png';


obstacles_loc=load_occupancy_map_ASIF('figure',occ_map_loc,'offset',occ_offset,'resolution',occ_resolution,'rotation',occ_rot,0);

figure;
figure;plot(obstacles_loc(:,1),obstacles_loc(:,2),'bs')
hold on;
plot(0,0,'ro','markerfacecolor','r','markersize',10)



% 

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


% -- fill holes
map = 1-map;

% https://se.mathworks.com/matlabcentral/newsreader/view_thread/173142
map = ~bwareaopen(~map, 100);

map = 1-map;


% map = bwmorph(map,'remove');
% map = bwmorph(map,'fill',5000);
% map = bwmorph(map,'majority',5);


figure; imshow(map','InitialMagnification',500); hold on;


filename = 'global_castings.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,map,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

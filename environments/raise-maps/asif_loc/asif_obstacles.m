clear;
clc;
close all;



map_spec_name='MapSpec_201510211300';
occ_map_loc='map_johnsson2.png';
obstacles_loc=load_occupancy_map_ASIF('figure',occ_map_loc,'spec_file',map_spec_name,'plot_map',0);

figure;plot(obstacles_loc(:,1),obstacles_loc(:,2),'bx')
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

% -- fll holes
map = 1-map;
% https://se.mathworks.com/matlabcentral/newsreader/view_thread/173142
map = ~bwareaopen(~map, 100);
map = 1-map;

figure; imshow(map','InitialMagnification',500); hold on;
set(gca,'YDir','normal');
colormap('gray'); 



filename = 'johnsson_metal.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,map,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);


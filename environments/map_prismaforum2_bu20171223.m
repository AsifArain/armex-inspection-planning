clear; clc; close all;

%--- parameters
%=====================================
occ_map_loc          = 'prismaforum2.pgm';
offset               = [0,0];
resolution           = 0.05;
rotation             = 0;
markercolor          = 'k';
markersize           = 1.5;
only_data            = 1;
desired_resolution_m = 0.50; %1.0; %0.25;
yaml_filename        = 'prismaforum2.yaml';

%--- map file
%=====================================
% [X,map_color] = imread(occ_map_loc);
[map_gray,map_color] = imread(occ_map_loc);

% Binary/BW map.
%map_BW = im2bw(X,map_color,0.95); % binary map
% map_BW = X;
% size(map_BW)
% figure('name','map_BW'); 
% imshow(map_BW);

figure('name','gray map'); imshow(map_gray); hold on;
plot(2000,2000,'*r')


%---- Map Cut
%=====================================
% [map_BW,mapCut] = imcrop;
% [X,Y,I2,mapCut] = imcrop;
% mapCut = 1.0e+03 * [1.9245    1.2365    2.6080 1.1120];
mapCut = 1.0e+03 * [0.6495    1.2755    2.8940    1.0320];
% Lets cut the map.
mapCut = round(mapCut);
map_gray(mapCut(2)+mapCut(4):end,:)   = [];
map_gray(1:mapCut(2),:)               = [];
map_gray(:,mapCut(1)+mapCut(3):end)   = [];
map_gray(:,1:mapCut(1))               = [];

figure('name','truncated gray map'); imshow(map_gray)



%-- read YAML file
%=====================================
dataYAML = ReadYamlRaw(yaml_filename);
dataYAML_origin = cell2mat(dataYAML.origin);
robot_origin_x  = abs(dataYAML_origin(1)/dataYAML.resolution);
robot_origin_y  = abs(dataYAML_origin(2)/dataYAML.resolution);





%-- map for raycasting (measurements)
%=====================================
pix_m = 1/dataYAML.resolution;
% Number of pixels in a cell for desired cell size.
pix_cell = desired_resolution_m*pix_m;

map_gray = flipdim(map_gray,1);


map_raycasting = map_gray';
map_raycasting = imbinarize(map_raycasting,0.9);


origin_raycasting = [round((robot_origin_x-mapCut(1))),...
                     round((mapCut(4)-(robot_origin_y-mapCut(2))))];

figure('name','ray casting map'); 
imshow(map_raycasting','InitialMagnification',500); hold on;
set(gca,'YDir','normal');
colormap('gray'); 
hold on;
plot(origin_raycasting(1),origin_raycasting(2),'*r')

filename = 'prismaforum2_map_raycasting.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,map_raycasting,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum2_origin_raycasting.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,origin_raycasting,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum2_mapsize_raycasting.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,size(map_raycasting),'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum2_cellsize_raycasting.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,resolution,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);


%==============================================================
%                       COVERAGE MAP
%==============================================================

%-- grow obstacles
r_cov  = 0.20/resolution;
SE_coverage = strel('disk',r_cov);
map_coverage = imerode(map_raycasting,SE_coverage);	% or BW2 = imdilate(BW, SE);

% figure('name','coverage map (without connectivity check)'); 
% imshow(map_coverage','InitialMagnification',500); hold on;
% set(gca,'YDir','normal');
% colormap('gray'); 
% hold on;
% plot(origin_raycasting(1),origin_raycasting(2),'*r')

%-- connectivity check
map_conf_labels = bwlabel(map_coverage);
% label number to be selected
label_cell_counts = zeros(numel(unique(map_conf_labels(:)))-1,1);
for i = 1:numel(unique(map_conf_labels(:)))-1
    label_cell_counts(i) = numel(find(map_conf_labels(:)==i));
end
[~,label_concerned] = max(label_cell_counts);    

label_ind = find(map_conf_labels(:)==label_concerned);
map_coverage = zeros(size(map_coverage));
map_coverage(label_ind) = 1;


figure('name','coverage map'); 
imshow(map_coverage','InitialMagnification',500); hold on;
set(gca,'YDir','normal');
colormap('gray'); 
hold on;
plot(origin_raycasting(1),origin_raycasting(2),'*r')


%-- save data

filename = 'prismaforum2_map_coverageL.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,map_coverage,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum2_origin_coverageL.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,origin_raycasting,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum2_mapsize_coverageL.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,size(map_coverage),'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum2_cellsize_coverageL.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,resolution,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

%-- rescale
resize_scale = resolution/desired_resolution_m;
map_coverage = imresize(map_coverage,resize_scale,'bilinear');
map_coverage = imbinarize(map_coverage,0.70); %
% map_coverage = imbinarize(map_coverage,'adaptive'); %

%-- connectivity check of the final map
map_conf_labels = bwlabel(map_coverage);
% label number to be selected
label_cell_counts = zeros(numel(unique(map_conf_labels(:)))-1,1);
for i = 1:numel(unique(map_conf_labels(:)))-1
    label_cell_counts(i) = numel(find(map_conf_labels(:)==i));
end
[~,label_concerned] = max(label_cell_counts);
label_ind = find(map_conf_labels(:)==label_concerned);
map_coverage = zeros(size(map_coverage));
map_coverage(label_ind) = 1;

%-- origin
origin_coverage = round(origin_raycasting*resize_scale);


figure('name','coverage map'); 
imshow(map_coverage','InitialMagnification',500); hold on;
set(gca,'YDir','normal');
colormap('gray'); 
hold on;
plot(origin_coverage(1),origin_coverage(2),'*r');


%-- save data

filename = 'prismaforum2_map_coverage.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,map_coverage,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum2_origin_coverage.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,origin_coverage,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum2_mapsize_coverage.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,size(map_coverage),'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum2_cellsize_coverage.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,desired_resolution_m,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

%{
%---- rescale
resize_scale = resolution/desired_resolution_m;

map_coverage1 = imresize(map_coverage,resize_scale,'bilinear');
map_coverage1 = imbinarize(map_coverage1,0.90); %

figure('name','coverage map 1'); 
imshow(map_coverage1','InitialMagnification',500); hold on;
set(gca,'YDir','normal');
colormap('gray'); 



[rows,cols]         = find(map_coverage~=1); %find(map_input==0);
obstacle_mat        = [rows(:),cols(:)];
obstacle_mat        = obstacle_mat.*resolution/desired_resolution_m;
offset_x            = offset(1);
offset_y            = offset(2);
obstacle_mat(:,1)   = obstacle_mat(:,1)+offset_x;
obstacle_mat(:,2)   = obstacle_mat(:,2)+offset_y;
RA                  =[cosd(rotation) -sind(rotation); sind(rotation) cosd(rotation)];
obstacle_mat        = obstacle_mat*RA;


origin_cov_map = [round((robot_origin_x-mapCut(1))/pix_cell),...
                  round((mapCut(4)-(robot_origin_y-mapCut(2)))/pix_cell)];
figure; hold on;
if only_data==1
    plot(obstacle_mat(:,1),obstacle_mat(:,2),...
        'sk','markerfacecolor',markercolor,'markersize',markersize);
end
axis equal


obstacles_loc = obstacle_mat;

figure;
plot(obstacles_loc(:,2),obstacles_loc(:,1),'bx')
hold on;
plot(origin_cov_map(1),origin_cov_map(2),'*r')

obstacles_loc2 = unique(round(obstacles_loc),'rows');

figure;
plot(obstacles_loc2(:,1),obstacles_loc2(:,2),'bs','MarkerFaceColor','b')
hold on;
plot(origin_cov_map(1),origin_cov_map(2),'*r')

axis equal

%pause
% obstacles_loc2 = obstacles_loc;
obstacles_loc2(:,1) = obstacles_loc2(:,1)-min(obstacles_loc2(:,1))+1;
obstacles_loc2(:,2) = obstacles_loc2(:,2)-min(obstacles_loc2(:,2))+1;

% figure;
% plot(obstacles_loc2(:,1),obstacles_loc2(:,2),'bs','MarkerFaceColor','b')
% hold on;
% plot(origin_cov_map(2),origin_cov_map(1),'*r')
% axis equal
% pause

% map_coverage = ones(max(obstacles_loc2(:,2)),max(obstacles_loc2(:,1)));
% obs_ind = sub2ind(size(map_coverage),obstacles_loc2(:,2),obstacles_loc2(:,1));
% map_coverage(obs_ind) = 0;

map_coverage2 = ones(max(obstacles_loc2(:,1)),max(obstacles_loc2(:,2)));
obs_ind = sub2ind(size(map_coverage2),obstacles_loc2(:,1),obstacles_loc2(:,2));
map_coverage2(obs_ind) = 0;

figure('name','coverage map2'); 
imshow(map_coverage2','InitialMagnification',500); hold on;
set(gca,'YDir','normal');
colormap('gray'); 
hold on;
plot(origin_cov_map(1),origin_cov_map(2),'*r')



figure('name','coverage map1-2'); 
imshow((map_coverage2-map_coverage1)','InitialMagnification',500); hold on;
set(gca,'YDir','normal');
colormap('gray'); 
hold on;
plot(origin_cov_map(1),origin_cov_map(2),'*r')

%}



%==============================================================
%                       NAVIGATION MAP
%==============================================================

%-- grow obstacles
r_nav  = 0.50/resolution;
SE_navigation = strel('disk',r_nav);
map_navigation = imerode(map_raycasting,SE_navigation);	% or BW2 = imdilate(BW, SE);

% figure('name','navigation map (without connectivity check)'); 
% imshow(map_navigation','InitialMagnification',500); hold on;
% set(gca,'YDir','normal');
% colormap('gray'); 
% hold on;
% plot(origin_raycasting(1),origin_raycasting(2),'*r')


%-- connectivity check
map_conf_labels = bwlabel(map_navigation);

% label number to be selected
label_cell_counts = zeros(numel(unique(map_conf_labels(:)))-1,1);
for i = 1:numel(unique(map_conf_labels(:)))-1
    label_cell_counts(i) = numel(find(map_conf_labels(:)==i));
end
[~,label_concerned] = max(label_cell_counts);    

label_ind = find(map_conf_labels(:)==label_concerned);
map_navigation = zeros(size(map_navigation));
map_navigation(label_ind) = 1;


figure('name','navigation map'); 
imshow(map_navigation','InitialMagnification',500); hold on;
set(gca,'YDir','normal');
colormap('gray'); 
hold on;
plot(origin_raycasting(1),origin_raycasting(2),'*r')


%-- rescale
resize_scale = resolution/desired_resolution_m;
map_navigation = imresize(map_navigation,resize_scale,'bilinear');
map_navigation = imbinarize(map_navigation,0.50); %

% %-- connectivity check of the final map
% map_conf_labels = bwlabel(map_navigation);
% % label number to be selected
% label_cell_counts = zeros(numel(unique(map_conf_labels(:)))-1,1);
% for i = 1:numel(unique(map_conf_labels(:)))-1
%     label_cell_counts(i) = numel(find(map_conf_labels(:)==i));
% end
% [~,label_concerned] = max(label_cell_counts);
% label_ind = find(map_conf_labels(:)==label_concerned);
% map_navigation = zeros(size(map_navigation));
% map_navigation(label_ind) = 1;

%-- origin
origin_navigation = round(origin_raycasting*resize_scale);


figure('name','navigation map'); 
imshow(map_navigation','InitialMagnification',500); hold on;
set(gca,'YDir','normal');
colormap('gray'); 
hold on;
plot(origin_navigation(1),origin_navigation(2),'*r');


%-- save data

filename = 'prismaforum2_map_navigation.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,map_navigation,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum2_origin_navigation.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,origin_navigation,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum2_mapsize_navigation.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,size(map_navigation),'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum2_cellsize_navigation.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,desired_resolution_m,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

%==============================================================
%                  SENSOR PLACEMENT MAP
%==============================================================

%-- grow obstacles
r_conf = 0.75/resolution;
SE_conf = strel('disk',r_conf);
map_conf = imerode(map_raycasting,SE_conf);	% or BW2 = imdilate(BW, SE);

% figure('name','conf map (without connectivity check)'); 
% imshow(map_conf','InitialMagnification',500); hold on;
% set(gca,'YDir','normal');
% colormap('gray'); 
% hold on;
% plot(origin_raycasting(1),origin_raycasting(2),'*r')

%-- connectivity check
map_conf_labels = bwlabel(map_conf);

% label number to be selected
label_cell_counts = zeros(numel(unique(map_conf_labels(:)))-1,1);
for i = 1:numel(unique(map_conf_labels(:)))-1
    label_cell_counts(i) = numel(find(map_conf_labels(:)==i));
end
[~,label_concerned] = max(label_cell_counts);    

label_ind = find(map_conf_labels(:)==label_concerned);
map_conf = zeros(size(map_conf));
map_conf(label_ind) = 1;

figure('name','conf map'); 
imshow(map_conf','InitialMagnification',500); hold on;
set(gca,'YDir','normal');
colormap('gray'); 
hold on;
% plot(124,335,'*r')


%---- rescale
resize_scale = resolution/desired_resolution_m;
map_conf = imresize(map_conf,resize_scale,'bilinear');
map_conf = imbinarize(map_conf,0.10); %

% %-- connectivity check of the final map
% map_conf_labels = bwlabel(map_conf);
% % label number to be selected
% label_cell_counts = zeros(numel(unique(map_conf_labels(:)))-1,1);
% for i = 1:numel(unique(map_conf_labels(:)))-1
%     label_cell_counts(i) = numel(find(map_conf_labels(:)==i));
% end
% [~,label_concerned] = max(label_cell_counts);
% label_ind = find(map_conf_labels(:)==label_concerned);
% map_conf = zeros(size(map_conf));
% map_conf(label_ind) = 1;

%-- origin
origin_conf = round(origin_raycasting*resize_scale);

figure('name','conf map'); 
imshow(map_conf','InitialMagnification',500); hold on;
set(gca,'YDir','normal');
colormap('gray'); 
hold on;
plot(origin_conf(1),origin_conf(2),'*r');

%-- save data

filename = 'prismaforum2_map_conf.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,map_conf,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum2_origin_conf.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,origin_conf,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum2_mapsize_conf.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,size(map_conf),'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum2_cellsize_conf.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,desired_resolution_m,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

% pause
% 
% %===================================== 
% %-- environment map for planning (OLD)
% %=====================================
% 
% % I=imread(occ_map_loc);
% % I_gray=map_gray;
% 
% % I_gray = imrotate(I_gray,180);
% % I_gray = flipdim(I_gray,2);
% % I_gray = flipdim(I_gray,1);
% 
% map_input           = double(map_gray);
% % map_input         = double(I);
% offset_x            = offset(1);
% offset_y            = offset(2);
% [rows,cols]         = find(map_input<210); %find(map_input==0);
% % obstacle_mat=[cols(:),rows(:)];
% obstacle_mat        = [rows(:),cols(:)];
% obstacle_mat        = obstacle_mat.*resolution/desired_resolution_m;
% obstacle_mat(:,1)   = obstacle_mat(:,1)+offset_x;
% obstacle_mat(:,2)   = obstacle_mat(:,2)+offset_y;
% RA                  =[cosd(rotation) -sind(rotation); sind(rotation) cosd(rotation)];
% obstacle_mat        = obstacle_mat*RA;
% 
% 
% % robotStartCellInMap = [round((robot_origin_x-mapCut(1))/pix_cell),...
% %                        round((robot_origin_y-mapCut(2))/pix_cell)];
% 
% % robotStartCellInMap = [round((robot_origin_x-mapCut(1))/pix_cell),...
% %                        round((mapCut(4)-(robot_origin_y-mapCut(2)))/pix_cell)];
%                    
% origin_env_map = [round((robot_origin_x-mapCut(1))/pix_cell),...
%                   round((mapCut(4)-(robot_origin_y-mapCut(2)))/pix_cell)];
% figure; hold on;                   
% if only_data==1
%     plot(obstacle_mat(:,1),obstacle_mat(:,2),...
%         'sk','markerfacecolor',markercolor,'markersize',markersize);
% end
% axis equal
% 
% obstacles_loc = obstacle_mat;
% 
% figure;
% plot(obstacles_loc(:,2),obstacles_loc(:,1),'bx')
% hold on;
% % plot(0,0,'ro','markerfacecolor','r','markersize',10)
% plot(origin_env_map(1),origin_env_map(2),'*r')
% % pause
% %/////////////////////////////////////////////////////////////////////////7
% % {
% obstacles_loc2 = unique(round(obstacles_loc),'rows');
% figure;
% plot(obstacles_loc2(:,1),obstacles_loc2(:,2),'bs','MarkerFaceColor','b')
% hold on;
% plot(origin_env_map(2),origin_env_map(1),'*r')
% %}
% axis equal
% 
% %pause
% % obstacles_loc2 = obstacles_loc;
% obstacles_loc2(:,1) = obstacles_loc2(:,1)-min(obstacles_loc2(:,1))+1;
% obstacles_loc2(:,2) = obstacles_loc2(:,2)-min(obstacles_loc2(:,2))+1;
% 
% figure;
% plot(obstacles_loc2(:,1),obstacles_loc2(:,2),'bs','MarkerFaceColor','b')
% hold on;
% plot(origin_env_map(2),origin_env_map(1),'*r')
% axis equal
% 
% % pause
% 
% map = ones(max(obstacles_loc2(:,2)),max(obstacles_loc2(:,1)));
% 
% obs_ind = sub2ind(size(map),obstacles_loc2(:,2),obstacles_loc2(:,1));
% map(obs_ind) = 0;
% 
% 
% 
% 
% figure; imshow(map','InitialMagnification',500); hold on;
% set(gca,'YDir','normal');
% colormap('gray'); 
% hold on;
% plot(origin_env_map(1),origin_env_map(2),'*r')
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % map = map
% 
% % ---------------- fll holes ------------------------
% % map = 1-map;
% % % https://se.mathworks.com/matlabcentral/newsreader/view_thread/173142
% % map = ~bwareaopen(~map, 100);
% % map = 1-map;plot(robotStartCellInMap(2),robotStartCellInMap(1),'*r')
% 
% % ---------- fill holes using 4-connectivity ------------------------
% %{
% CC = bwconncomp(map,4);
% conn_sizes = zeros(CC.NumObjects,1);
% for i = 1:CC.NumObjects
%     conn_sizes(i) = size(CC.PixelIdxList{i},1);
% end
% [val,ind] = max(conn_sizes);
% map = zeros(size(map));
% map(CC.PixelIdxList{ind}) = 1;
% %}
% 
% %close all
% 
% figure; imshow(map','InitialMagnification',500); hold on;
% set(gca,'YDir','normal');
% colormap('gray'); 
% plot(origin_env_map(1),origin_env_map(2),'*r')
% 
% 
% filename = 'prismaforum.dat';
% fileID = fopen(filename,'wt'); fclose(fileID);
% dlmwrite(filename,map,'-append','delimiter',' ');
% fileID = fopen(filename,'a'); fclose(fileID);
% 
% 
% filename = 'prismaforum_origin_map.dat';
% fileID = fopen(filename,'wt'); fclose(fileID);
% dlmwrite(filename,origin_env_map,'-append','delimiter',' ');
% fileID = fopen(filename,'a'); fclose(fileID);
% 
% 
% 
% % ------------
% % close all
% 
% map = dlmread('prismaforum.dat');
% figure; imshow(map','InitialMagnification',500); hold on;
% set(gca,'YDir','normal');
% colormap('gray'); 
% plot(origin_env_map(1),origin_env_map(2),'*r')
% 
% 


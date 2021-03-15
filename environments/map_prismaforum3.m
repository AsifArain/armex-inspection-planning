clear; clc; close all;

%--- parameters
%==============================================================
map_filename  = 'prismaforum3.pgm';
yaml_filename = 'prismaforum3.yaml';
%==============================================================


%-- read YAML file
%==============================================================
dataYAML        = ReadYamlRaw(yaml_filename);
dataYAML_origin = cell2mat(dataYAML.origin);
robot_origin_x  = abs(dataYAML_origin(1)/dataYAML.resolution);
robot_origin_y  = abs(dataYAML_origin(2)/dataYAML.resolution);
%==============================================================



%-- resolutions
%==============================================================
cellsize_raycast    = dataYAML.resolution;
cellsize_navigation = 0.50; 
cellsize_conf       = 0.50;
cellsize_coverage   = 0.50;
cellsize_coverageL  = dataYAML.resolution;



%--- map file
%==============================================================
[map_gray,map_color] = imread(map_filename);

figure('name','gray map'); imshow(map_gray); hold on;
plot(2000,2000,'*r')


%---- Map Cut
%==============================================================
% [map_BW,mapCut] = imcrop;
% [X,Y,I2,mapCut] = imcrop;
mapCut = 1.0e+03*[0.6375,1.6725,2.9040,0.5690];
% mapCut = 1.0e+03*[0.6495,1.5755,2.8940,1.0320];
% Lets cut the map.
mapCut = round(mapCut);
map_gray(mapCut(2)+mapCut(4):end,:)   = [];
map_gray(1:mapCut(2),:)               = [];
map_gray(:,mapCut(1)+mapCut(3):end)   = [];
map_gray(:,1:mapCut(1))               = [];

figure('name','truncated gray map'); imshow(map_gray)


%==============================================================
%                   ORIGINAL MAP (BINARY)
%==============================================================

map_gray     = flip(map_gray,1);
map_gray     = map_gray';
map_original = imbinarize(map_gray,0.9);

%-- origin
origin_original = [round((robot_origin_x-mapCut(1))),...
                   round((mapCut(4)-(robot_origin_y-mapCut(2))))];

%-- publish               
figure('name','original map'); 
imshow(map_original','InitialMagnification',500); hold on;
set(gca,'YDir','normal');
colormap('gray'); 
hold on;
plot(origin_original(1),origin_original(2),'*r')

%-- save
filename = 'prismaforum3_map_original.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,map_original,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum3_origin_original.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,origin_original,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum3_mapsize_original.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,size(map_original),'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum3_cellsize_original.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,dataYAML.resolution,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);


%==============================================================
%                     RAYCASTING MAP
%==============================================================

%-- rescale
resize_factor_ray = dataYAML.resolution/cellsize_raycast;
map_raycast       = imresize(map_original,resize_factor_ray,'bilinear');
% map_raycast       = imbinarize(map_raycast,0.50); %

%-- origin
origin_raycast = round(origin_original*resize_factor_ray);

%-- publish
figure('name','ray casting map'); 
imshow(map_raycast','InitialMagnification',500); hold on;
set(gca,'YDir','normal');
colormap('gray'); 
hold on;
plot(origin_raycast(1),origin_raycast(2),'*r')

%-- save
filename = 'prismaforum3_map_raycast.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,map_raycast,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum3_origin_raycast.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,origin_raycast,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum3_mapsize_raycast.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,size(map_raycast),'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum3_cellsize_raycast.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,cellsize_raycast,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);


%==============================================================
%                       COVERAGE MAP
%==============================================================

%-- grow obstacles
r_cov  = 0.20/dataYAML.resolution;
SE_coverage = strel('disk',r_cov);
map_coverage = imerode(map_original,SE_coverage);	% or BW2 = imdilate(BW, SE);

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

label_ind = map_conf_labels(:)==label_concerned;
map_coverage = zeros(size(map_coverage));
map_coverage(label_ind) = 1;

%-- publish
figure('name','coverage map (large)'); 
imshow(map_coverage','InitialMagnification',500); hold on;
set(gca,'YDir','normal');
colormap('gray'); 
hold on;
plot(origin_original(1),origin_original(2),'*r')


%-- save data

filename = 'prismaforum3_map_coverageL.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,map_coverage,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum3_origin_coverageL.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,origin_original,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum3_mapsize_coverageL.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,size(map_coverage),'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum3_cellsize_coverageL.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,dataYAML.resolution,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);



%-- rescale
resize_factor_coverage = dataYAML.resolution/cellsize_coverage;
map_coverage = imresize(map_coverage,resize_factor_coverage,'bilinear');
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
label_ind = map_conf_labels(:)==label_concerned;
map_coverage = zeros(size(map_coverage));
map_coverage(label_ind) = 1;

%-- origin
origin_coverage = round(origin_original*resize_factor_coverage);

%-- publish
figure('name','coverage map'); 
imshow(map_coverage','InitialMagnification',500); hold on;
set(gca,'YDir','normal');
colormap('gray'); 
hold on;
plot(origin_coverage(1),origin_coverage(2),'*r');


%-- save data

filename = 'prismaforum3_map_coverage.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,map_coverage,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum3_origin_coverage.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,origin_coverage,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum3_mapsize_coverage.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,size(map_coverage),'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum3_cellsize_coverage.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,cellsize_coverage,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);


%==============================================================
%                       NAVIGATION MAP
%==============================================================

%-- grow obstacles
r_nav  = 0.50/dataYAML.resolution;
SE_navigation = strel('disk',r_nav);
map_navigation = imerode(map_original,SE_navigation);	% or BW2 = imdilate(BW, SE);

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

label_ind = map_conf_labels(:)==label_concerned;
map_navigation = zeros(size(map_navigation));
map_navigation(label_ind) = 1;

%-- publish
figure('name','navigation map'); 
imshow(map_navigation','InitialMagnification',500); hold on;
set(gca,'YDir','normal');
colormap('gray'); 
hold on;
plot(origin_original(1),origin_original(2),'*r')


%-- rescale
resize_factor_navigation = dataYAML.resolution/cellsize_navigation;
map_navigation = imresize(map_navigation,resize_factor_navigation,'bilinear');
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
origin_navigation = round(origin_original*resize_factor_navigation);


figure('name','navigation map'); 
imshow(map_navigation','InitialMagnification',500); hold on;
set(gca,'YDir','normal');
colormap('gray'); 
hold on;
plot(origin_navigation(1),origin_navigation(2),'*r');


%-- save data

filename = 'prismaforum3_map_navigation.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,map_navigation,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum3_origin_navigation.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,origin_navigation,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum3_mapsize_navigation.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,size(map_navigation),'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum3_cellsize_navigation.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,cellsize_navigation,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

%==============================================================
%                  SENSOR PLACEMENT MAP
%==============================================================

%-- grow obstacles
r_conf = 0.75/dataYAML.resolution;
SE_conf = strel('disk',r_conf);
map_conf = imerode(map_original,SE_conf);	% or BW2 = imdilate(BW, SE);

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
resize_factor_conf = dataYAML.resolution/cellsize_conf;
map_conf = imresize(map_conf,resize_factor_conf,'bilinear');
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
origin_conf = round(origin_original*resize_factor_conf);

figure('name','conf map'); 
imshow(map_conf','InitialMagnification',500); hold on;
set(gca,'YDir','normal');
colormap('gray'); 
hold on;
plot(origin_conf(1),origin_conf(2),'*r');

%-- save data

filename = 'prismaforum3_map_conf.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,map_conf,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum3_origin_conf.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,origin_conf,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum3_mapsize_conf.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,size(map_conf),'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'prismaforum3_cellsize_conf.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,cellsize_conf,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);



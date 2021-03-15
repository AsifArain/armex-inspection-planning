clear; clc; close all;

%--- parameters
%=====================================
occ_map_loc          = 'lab.pgm';
offset               = [0,0];
resolution           = 0.05;
rotation             = 0;
markercolor          = 'k';
markersize           = 1.5;
only_data            = 1;
desired_resolution_m = 0.25;
yaml_filename        = 'lab.yaml';

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
%[map_BW,mapCut] = imcrop;
%[X,Y,I2,mapCut] = imcrop;
mapCut = 1.0e+03 *[1.9385    1.8675    0.3120    0.2220];
% Lets cut the map.
mapCut = round(mapCut);
map_gray(mapCut(2)+mapCut(4):end,:)   = [];
map_gray(1:mapCut(2),:)               = [];
map_gray(:,mapCut(1)+mapCut(3):end)   = [];
map_gray(:,1:mapCut(1))               = [];
figure('name','truncated gray map'); imshow(map_gray)


pause



%map_BW = bwmorph(map_BW,'diag');
%figure('name','map_BW with morpho'); 
%imshow(map_BW);

% If you want to cut the map.
%map = ~map;
%{
=== for Lab1
map_BW(2060:end,:) = [];
figure; imshow(map_BW);
map_BW(1:1680,:)   = [];
figure; imshow(map_BW);
map_BW(:,2180:end) = [];
figure; imshow(map_BW);
map_BW(:,1:1900)   = [];
figure; imshow(map_BW);
%}
%pause

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ---- Map Cut -----
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %[map_BW,mapCut] = imcrop;
% %[X,Y,I2,mapCut] = imcrop;
% mapCut = 1.0e+03 *[1.9385    1.8675    0.3120    0.2220];
% % Lets cut the map.
% mapCut = round(mapCut)
% map_BW(mapCut(2)+mapCut(4):end,:)   = [];
% map_BW(1:mapCut(2),:)               = [];
% map_BW(:,mapCut(1)+mapCut(3):end)   = [];
% map_BW(:,1:mapCut(1))               = [];
% figure; imshow(map_BW)


%-- read YAML file
dataYAML = ReadYamlRaw(yaml_filename);
dataYAML_origin = cell2mat(dataYAML.origin);
robot_origin_x = abs(dataYAML_origin(1)/dataYAML.resolution);
robot_origin_y = abs(dataYAML_origin(2)/dataYAML.resolution);

pix_m = 1/dataYAML.resolution;
% Number of pixels in a cell for desired cell size.
pix_cell = desired_resolution_m*pix_m;

map_BW = flipdim(map_BW,1);

map_raycasting = map_BW';
map_raycasting = imbinarize(map_raycasting,0.9);


origin_raycasting = [round((robot_origin_x-mapCut(1))),...
                     round((mapCut(4)-(robot_origin_y-mapCut(2))))];

figure; imshow(map_raycasting','InitialMagnification',500); hold on;
set(gca,'YDir','normal');
colormap('gray'); 
hold on;
plot(origin_raycasting(1),origin_raycasting(2),'*r')

filename = 'lab_map_raycasting.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,map_raycasting,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'lab_origin_raycasting.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,origin_raycasting,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

pause

% I=imread(occ_map_loc);
I_gray=map_BW;

% I_gray = imrotate(I_gray,180);
% I_gray = flipdim(I_gray,2);
% I_gray = flipdim(I_gray,1);

map_input=double(I_gray);

% map_input=double(I);
offset_x=offset(1);
offset_y=offset(2);
[rows,cols]=find(map_input==0);
% obstacle_mat=[cols(:),rows(:)];
obstacle_mat=[rows(:),cols(:)];
obstacle_mat=obstacle_mat.*resolution/desired_resolution_m;
obstacle_mat(:,1)=obstacle_mat(:,1)+offset_x;
obstacle_mat(:,2)=obstacle_mat(:,2)+offset_y;
RA=[cosd(rotation) -sind(rotation); sind(rotation) cosd(rotation)];
obstacle_mat=obstacle_mat*RA;





% robotStartCellInMap = [round((robot_origin_x-mapCut(1))/pix_cell),...
%                        round((robot_origin_y-mapCut(2))/pix_cell)];

robotStartCellInMap = [round((robot_origin_x-mapCut(1))/pix_cell),...
                       round((mapCut(4)-(robot_origin_y-mapCut(2)))/pix_cell)];

if only_data==1
    plot(obstacle_mat(:,1),obstacle_mat(:,2),...
        'sk','markerfacecolor',markercolor,'markersize',markersize);
end



obstacles_loc = obstacle_mat;

figure;
plot(obstacles_loc(:,2),obstacles_loc(:,1),'bx')
hold on;
% plot(0,0,'ro','markerfacecolor','r','markersize',10)
plot(robotStartCellInMap(1),robotStartCellInMap(2),'*r')


%/////////////////////////////////////////////////////////////////////////7

obstacles_loc2 = unique(round(obstacles_loc),'rows');

figure;
plot(obstacles_loc2(:,1),obstacles_loc2(:,2),'bs','MarkerFaceColor','b')
hold on;
plot(robotStartCellInMap(2),robotStartCellInMap(1),'*r')
axis equal


obstacles_loc2(:,1) = obstacles_loc2(:,1)-min(obstacles_loc2(:,1))+1;
obstacles_loc2(:,2) = obstacles_loc2(:,2)-min(obstacles_loc2(:,2))+1;

figure;
plot(obstacles_loc2(:,1),obstacles_loc2(:,2),'bs','MarkerFaceColor','b')
hold on;
plot(robotStartCellInMap(2),robotStartCellInMap(1),'*r')
axis equal



map = ones(max(obstacles_loc2(:,2)),max(obstacles_loc2(:,1)));

obs_ind = sub2ind(size(map),obstacles_loc2(:,2),obstacles_loc2(:,1));
map(obs_ind) = 0;



% map = map

% ---------------- fll holes ------------------------
% map = 1-map;
% % https://se.mathworks.com/matlabcentral/newsreader/view_thread/173142
% map = ~bwareaopen(~map, 100);
% map = 1-map;plot(robotStartCellInMap(2),robotStartCellInMap(1),'*r')

% ---------- fill holes using 4-connectivity ------------------------
CC = bwconncomp(map,4);
conn_sizes = zeros(CC.NumObjects,1);
for i = 1:CC.NumObjects
    conn_sizes(i) = size(CC.PixelIdxList{i},1);
end
[val,ind] = max(conn_sizes);
map = zeros(size(map));
map(CC.PixelIdxList{ind}) = 1;


close all

figure; imshow(map','InitialMagnification',500); hold on;
set(gca,'YDir','normal');
colormap('gray'); 
plot(robotStartCellInMap(2),robotStartCellInMap(1),'*r')


filename = 'lab.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,map,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);



% ------------
close all

map = dlmread('lab.dat');
figure; imshow(map','InitialMagnification',500); hold on;
set(gca,'YDir','normal');
colormap('gray'); 



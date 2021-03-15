%----------------------------------------------------------------------
%
%       SAMPLE MAP FOR THE RECONSTRUCTION FIGURE TO BE USED IN THE JFR
%       PAPER
%
%-----------------------------------------------------------------------
close all; clear all; clc;

map = ones(100,60);
map(1:end,1)   = 0;
map(1:end,end) = 0;
map(1,1:end)   = 0;
map(end,1:end) = 0;

map(45,20:30) = 0;
map(46,20:30) = 0;

map(45:75,30) = 0;
map(45:75,29) = 0;


cellsize = 1;
origin = [10,10];
mapcut = [0 0 100 60];


% MAPS
%-------------------------------
filename = 'sample-env-paper.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,map,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'sample-env-paper_map_conf.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,map,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'sample-env-paper_map_coverage.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,map,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'sample-env-paper_map_raycast.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,map,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);


% CELL SIZE
%--------------------------------
filename = 'sample-env-paper_cellsize_conf.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,cellsize,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'sample-env-paper_cellsize_coverage.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,cellsize,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'sample-env-paper_cellsize_original.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,cellsize,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'sample-env-paper_cellsize_raycast.dat'; 
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,cellsize,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

% ORIGIN
%----------------------------------
filename = 'sample-env-paper_origin_conf.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,origin,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'sample-env-paper_origin_coverage.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,origin,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'sample-env-paper_origin_raycast.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,origin,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

% MAP CUT
%----------------------------------
filename = 'sample-env-paper_mapcut_original.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,mapcut,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);


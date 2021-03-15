%----------------------------------------------------------------------
%
%       SAMPLE MAP FOR THE RECONSTRUCTION FIGURE TO BE USED IN THE JFR
%       PAPER
%
%-----------------------------------------------------------------------
close all; clear all; clc;

map = ones(50);
cellsize = 1;
origin = [24,6];
mapcut = [0 0 50 50];


% MAPS
%-------------------------------
filename = 'sample-reconstruction.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,map,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'sample-reconstruction_map_conf.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,map,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'sample-reconstruction_map_coverage.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,map,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'sample-reconstruction_map_raycast.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,map,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);


% CELL SIZE
%--------------------------------
filename = 'sample-reconstruction_cellsize_conf.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,cellsize,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'sample-reconstruction_cellsize_coverage.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,cellsize,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'sample-reconstruction_cellsize_original.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,cellsize,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'sample-reconstruction_cellsize_raycast.dat'; 
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,cellsize,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

% ORIGIN
%----------------------------------
filename = 'sample-reconstruction_origin_conf.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,origin,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'sample-reconstruction_origin_coverage.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,origin,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

filename = 'sample-reconstruction_origin_raycast.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,origin,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);

% MAP CUT
%----------------------------------
filename = 'sample-reconstruction_mapcut_original.dat';
fileID = fopen(filename,'wt'); fclose(fileID);
dlmwrite(filename,mapcut,'-append','delimiter',' ');
fileID = fopen(filename,'a'); fclose(fileID);


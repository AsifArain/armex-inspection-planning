close all; clear all; clc;


% map = dlmread('corridor.dat');
% map = 1-map;
% 
% figure; imshow(map','InitialMagnification',500); hold on;
% %set(gca,'YDir','normal');
% colormap('gray'); 
% 
% ------------------------------------

map = dlmread('johnsson-metal.dat');

figure; imshow(map','InitialMagnification',500); hold on;
set(gca,'YDir','normal');
colormap('gray'); 




% ------------------------------------
% map = dlmread('global-castings.dat');
% 
% figure; imshow(map','InitialMagnification',500); hold on;
% set(gca,'YDir','normal');
% colormap('gray'); 
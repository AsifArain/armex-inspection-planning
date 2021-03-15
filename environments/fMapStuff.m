function [map_env,...
          map_conf,...
          origin_env,...
          origin_conf,...
          cellsize_env,...
          cellsize_conf,...
          cellsize_org] = fMapStuff(para_,dir_)

%************************************************
%               MAP AND REALTED:
%************************************************
file__map_env   = sprintf('%s_map_coverage.dat',para_.environment);
file__map_conf  = sprintf('%s_map_conf.dat',para_.environment);
file__org_env   = sprintf('%s_origin_coverage.dat',para_.environment);
file__org_conf  = sprintf('%s_origin_conf.dat',para_.environment);
file__cell_env  = sprintf('%s_cellsize_coverage.dat',para_.environment);
file__cell_conf = sprintf('%s_cellsize_conf.dat',para_.environment);
file__cell_org  = sprintf('%s_cellsize_original.dat',para_.environment);
file__cut_org   = sprintf('%s_mapcut_original.dat',para_.environment);
map_env         = load([dir_.env,file__map_env]);
map_conf        = load([dir_.env,file__map_conf]);
origin_env      = load([dir_.env,file__org_env]);
origin_conf     = load([dir_.env,file__org_conf]);
cellsize_env    = load([dir_.env,file__cell_env]);
cellsize_conf   = load([dir_.env,file__cell_conf]);
cellsize_org    = load([dir_.env,file__cell_org]);
% mapCut          = load([dir_.env,file__cut_org]);


end
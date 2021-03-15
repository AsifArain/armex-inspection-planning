function fPubEnvironment(map_env,para_,pubPara_)
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             ENVIRONMENT MAP:
% 
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% map_env = dlmread([dir_.env,para_.EnvironmentFigFile]);
% map_env = 1-map_env;
% map_env = dataPublish.map_env;
% [map_row,map_col] = size(map_env);

h = figure('name',para_.PublishFigureTitle); 
imshow(map_env','InitialMagnification',pubPara_.zoom_in_environment);
hold on;
set(gca,'YDir','normal');
colormap('gray'); 

%caxis([-5 2]) % for resulting gdm
%caxis([-12 4]) % for sensor placement solution
% caxis([-0.25 1]) % for sensor placement solution
% caxis([-2 2])
caxis([-1 1.25])
%caxis([-10 1.5]) % for gdm insets

end
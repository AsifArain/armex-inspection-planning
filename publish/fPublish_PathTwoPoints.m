function [thisVideo] = fPublish_PathTwoPoints( start_position,...
                                               end___position,...
                                               thisVideo,...
                                               map_env,...
                                               paraPub_ )
% Publish path between two points and write it to a video
% @Asif Arain

%-- A-star path
% round(start_position)
% round(end___position)
[ path_astar ] = fPath2(round(start_position),...
                        round(end___position),...
                        map_env );

% sample_rate = 0.5;%0.05;
% num_samples = round(size(path_astar,1)*sample_rate);
num_astar_path_points = round(size(path_astar,1)*paraPub_.astar_path_sampling_fraction);
if (num_astar_path_points) > 2

    %path_sample = path(unique(randperm(size(path,1),num_samples)),:)
    random_selection = 1+randperm(size(path_astar,1)-2);
    path_sample = path_astar(unique(random_selection(1:num_astar_path_points)),:);
    path_sample = [path_astar(1,:);path_sample;path_astar(end,:)];
    path_pts = interparc(size(path_astar,1)*2,path_sample(:,1),path_sample(:,2),'spline');
else
    path_pts = path_astar;
end


for j = 1:size(path_pts,1)
    hPath(j) = plot(path_pts(j,1)-1.0,path_pts(j,2)-1.0,'.','color',paraPub_.colorPath);
    %-- write video
    writeVideo(thisVideo,getframe);
end

end

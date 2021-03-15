function [ path_astar ] = fPath2(start_position,goal_position,map_env )
%
%*************************************************************************
% Traveling path between two positions.
%*************************************************************************
%
% @AsifArain - 18-Oct-2018

Connecting_Distance = 1;

StartX = round(start_position(1,2));
StartY = round(start_position(1,1));
MAP    = 1-map_env;

GoalRegister = int8(zeros(size(MAP)));
GoalRegister(round(goal_position(1,1)),round(goal_position(1,2))) = 1;

if pdist2(round(start_position(1,:)),round(goal_position(1,:)),'euclidean') >= 1
    OptimalPath = ASTARPATH(StartX,StartY,MAP,GoalRegister,Connecting_Distance);
end

% path = [goal_position;OptimalPath;start_position]
% path = OptimalPath;
path_astar = flipud(OptimalPath);

end

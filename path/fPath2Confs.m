function [ OptimalPath ] = fPath2Confs(start_position,goal_position,map_env )
%
%*************************************************************************
% Traveling path between two positions.
%*************************************************************************
%
% @AsifArain - 07-Apr-2018



Connecting_Distance = 1;

StartX = round(start_position(1,2));
StartY = round(start_position(1,1));
MAP    = 1-map_env;

GoalRegister = int8(zeros(size(MAP)));
GoalRegister(round(goal_position(1,1)),round(goal_position(1,2))) = 1;


if pdist2(round(start_position(1,1:2)),round(goal_position(1,1:2)),'euclidean') >= 1
    
    OptimalPath = ASTARPATH(StartX,StartY,MAP,GoalRegister,Connecting_Distance);    
    
else
    OptimalPath = [];
end




   
end

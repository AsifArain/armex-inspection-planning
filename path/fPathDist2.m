function [ dists ] = fPathDist2(start_positions,goal_positions,map_env )
%
%*************************************************************************
% Traveling path distance matrix between two lists of sensing positions.
% Rows are start poistion and columns are end positions.
%*************************************************************************
%
% @AsifArain - 07-Apr-2018

Connecting_Distance = 1;
%dists = zeros(size(start_positions,1),size(goal_positions,1));
dists = sparse(size(start_positions,1),size(goal_positions,1));


for i = 1:size(dists,1)    
    parfor j = 1:size(dists,2)        
        
        
        %fprintf('this i and j are: %d and %d\n',i,j);
        
        StartX = round(start_positions(i,2));
        StartY = round(start_positions(i,1));
        MAP    = 1-map_env;
        
        %size(MAP)

        GoalRegister = int8(zeros(size(MAP)));
        GoalRegister(round(goal_positions(j,1)),round(goal_positions(j,2))) = 1;

        if pdist2(round(start_positions(i,:)),round(goal_positions(j,:)),'euclidean') >= 1
            OptimalPath = ASTARPATH(StartX,StartY,MAP,GoalRegister,Connecting_Distance);            
            %dists(i,j) = size(OptimalPath,1)-1;
            
            %-- Manhattan dist one is one step length and manhattan dist
            %two is 1.4142 length
            steps = sum(abs(diff(OptimalPath,1,1)),2);            
            dists(i,j) = (numel(find(steps==1))*1) + (numel(find(steps==2))*1.4142);
            
        else
            dists(i,j) = 0;
        end
    end
end
end

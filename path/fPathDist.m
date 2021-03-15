function [ dists ] = fPathDist(start_position,goal_positions,map_env )
%


Connecting_Distance = 1;
dists = zeros(size(goal_positions,1),1);

for i = 1:size(dists,1)
        
    StartX = round(start_position(2));
    StartY = round(start_position(1));
    MAP    = 1-map_env;
    
    GoalRegister = int8(zeros(size(MAP)));
    GoalRegister(round(goal_positions(i,1)),round(goal_positions(i,2))) = 1;
    
    %goal_positions(i,1)
    %goal_positions(i,2)    
    %GoalRegister    
    %[goal_r,goal_c] = find(GoalRegister)
    
    
    if pdist2(round(start_position),round(goal_positions(i,:)),'euclidean') >= 1
        
        %disp('im here')
    
        OptimalPath = ASTARPATH(StartX,StartY,MAP,GoalRegister,Connecting_Distance);
        dists(i) = size(OptimalPath,1)-1;
        
        %{
        hold on
        plot(OptimalPath(1,1),OptimalPath(1,2),'o','color','k')
        plot(OptimalPath(end,1),OptimalPath(end,2),'o','color','b')
        plot(OptimalPath(:,1),OptimalPath(:,2),'r')
        legend('Goal','Start','Path')
        %}
    
    else
        dists(i) = 0;
    end
    
    


    
    %pause
end



% %Start Positions
% StartX=15;
% StartY=15;
% 
% %Generating goal nodes, which is represented by a matrix. Several goals can be speciefied, in which case the pathfinder will find the closest goal. 
% %a cell with the value 1 represent a goal cell
% GoalRegister=int8(zeros(128,140));
% GoalRegister(110,80)=1;
% 
% %Number of Neighboors one wants to investigate from each cell. A larger
% %number of nodes means that the path can be alligned in more directions. 
% %Connecting_Distance=1-> Path can  be alligned along 8 different direction.
% %Connecting_Distance=2-> Path can be alligned along 16 different direction.
% %Connecting_Distance=3-> Path can be alligned along 32 different direction.
% %Connecting_Distance=4-> Path can be alligned along 56 different direction.
% %ETC......
% 
% Connecting_Distance=1; %Avoid to high values Connecting_Distances for reasonable runtimes. 
% 
% % Running PathFinder
% OptimalPath=ASTARPATH(StartX,StartY,MAP,GoalRegister,Connecting_Distance)


end

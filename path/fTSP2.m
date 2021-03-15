function [ tRoute,...
           tCost,...
           tTSP ] = fTSP2( o,...
                           FoVfaces,...
                           C,...
                           E,...
                           T,...
                           map_env,...
                           OBS,...
                           ParamTSP )
%fTSP generates travelling route for THE Art Gallery Problem.
% Date: 2014-01-17, Rev 1
% 
% OUTPUTS
% --------------------------------------------------------------------------------------------------------------------------------
% tRoute    :
% TRoute    :
% tCost     :
% tTSP      :
% 
% INPUTS:
% --------------------------------------------------------------------------------------------------------------------------------
% o         : [scalar] number of outgoing edges 
% FoVfaces  :
% C         :
% E         :
% T         :
% map       :


tTSPi = tic;

% CONSTANT:
n = numel((map_env));%numel(find(map));
f = FoVfaces.num;

% Index numbers of Conf.
Conf = find(C);
for i = 1:numel(OBS.ind)
    Conf(Conf>((OBS.ind(i)-1)*f)) = Conf(Conf>((OBS.ind(i)-1)*f))+f;
end
% Conf
% size(Conf)
% pause
% Cell number of each Conf.
iCELL = fix(Conf/f)+(~~mod(Conf,f));
% 
[rCELL,cCELL] = ind2sub(size(map_env),iCELL);



%% TSP:
xy   = [rCELL cCELL];
dmat = inf(numel(Conf));
Pred = cell(numel(Conf),numel(Conf));

if ParamTSP.EoTinAstar
    Eo             = [];
    ParamAstar.EoT = 0;
    
elseif ~ParamTSP.EoTinAstar
    % Eo & T:
    Eo    = zeros(size(T'));
    i     = 1:length(E);
    Eo(i) = fix(E(i)/o)+(~~mod(E(i),o));
    Eo    = Eo';
    % Fixing Eo and T for not connected Conf(s).
    i = 1:length(Eo(:,1));
    T((fix(Eo(i,1)/f)+(~~mod(Eo(i,1),f)))'==fix(i/f)+(~~mod(i,f)),1)  = inf;
    Eo((fix(Eo(i,1)/f)+(~~mod(Eo(i,1),f)))'==fix(i/f)+(~~mod(i,f)),1) = inf;

    % Insert obstacles.
    for i = 1:numel(OBS.ind)
        T = [T(1:(OBS.ind(i)-1)*f,:); inf(f,length(T(1,:))); T((OBS.ind(i)-1)*f+1:end,:)];
        Eo = [Eo(1:(OBS.ind(i)-1)*f,:); inf(f,length(Eo(1,:))); Eo((OBS.ind(i)-1)*f+1:end,:)];    
    end
    % Rearrange IDs after inserting Obstacles.
    for i = 1:numel(OBS.ind)
        Eo(Eo>((OBS.ind(i)-1)*f))= Eo(Eo>((OBS.ind(i)-1)*f))+f;
    end
    ParamAstar.EoT = 1;
end
ParamAstar.TravCost = ParamTSP.TravCost;
ParamAstar.ConfOBS  = ParamTSP.ConfOBS;

% t_init = tic;
% % for i = 1:numel(Conf)
% %     for j = 1:numel(Conf)        
% %         [ gCost,gPred,~ ] = fAstar( n,f,o,Conf(i),Conf(j),E,Eo,T,map,OBS,ParamAstar );
% %         dmat(i,j) = gCost;
% %         Pred{i,j} = gPred;
% %     end
% % end
% % dmat
% comp_time =  1e-4*(round(toc(t_init)*1e4));
% disp(['Compuation time -- A* old: ',num2str(comp_time),' sec']);


t_init = tic;
% for i = 1:numel(Conf)    
%     xy = [rCELL cCELL];    
%     dists = fPathDist(xy(i,:),xy,map_env);
%     %[ gCost,gPred,~ ] = fAstar( n,f,o,Conf(i),Conf(j),E,Eo,T,map,OBS,ParamAstar );
%     %dmat(i,j) = gCost;
%     %Pred{i,j} = gPred;
%     dmat(i,:) = dists;
% end
dmat = fPathDist2(xy,xy,map_env);
comp_time =  1e-4*(round(toc(t_init)*1e4));
disp(['Compuation time (A*): ',num2str(comp_time),' sec']);


popSize     = numel(Conf);
showProg    = 0;
showResult  = 0;
% xy
% dmat
% popSize
% showProg
% showResult
[optRoute,tCost] = tsp_nn(xy,dmat,popSize,showProg,showResult);

tRoute = Conf(optRoute);
tRoute = [tRoute; tRoute(1)]; % traveling route with conf IDs/positions.
% 
% TRoute = [];
% optRoute = [optRoute optRoute(1)];
% for i = 1:numel(tRoute)-1
%     TRoute = [TRoute fliplr(Pred{optRoute(i),optRoute(i+1)})];
% end
% TRoute = [TRoute tRoute(end)]';
% dmatTSP = dmat;
% PredTSP = Pred;


%**********************************************************************
% Total computation time.
%**********************************************************************
tTSP = 1e-4*(round(toc(tTSPi)*1e4));

fprintf(1,'\nComputation time: %0.4f sec \n',tTSP)

end

%TSP_NN Traveling Salesman Problem (TSP) Nearest Neighbor (NN) Algorithm
%   The Nearest Neighbor algorithm produces different results depending on
%   which city is selected as the starting point. This function determines
%   the Nearest Neighbor routes for multiple starting points and returns
%   the best of those routes
%
% Summary:
%     1. A single salesman travels to each of the cities and completes the
%        route by returning to the city he started from
%     2. Each city is visited by the salesman exactly once
%
% Input:
%     XY (float) is an Nx2 matrix of city locations, where N is the number of cities
%     DMAT (float) is an NxN matrix of point to point distances/costs
%     POPSIZE (scalar integer) is the size of the population (should be <= N)
%     SHOWPROG (scalar logical) shows the GA progress if true
%     SHOWRESULT (scalar logical) shows the GA results if true
%
% Output:
%     OPTROUTE (integer array) is the best route found by the algorithm
%     MINDIST (scalar float) is the cost of the best route
%
% Example:
%     n = 50;
%     xy = 10*rand(n,2);
%     popSize = n;
%     showProg = 1;
%     showResult = 1;
%     a = meshgrid(1:n);
%     dmat = reshape(sqrt(sum((xy(a,:)-xy(a',:)).^2,2)),n,n);
%     [optRoute,minDist] = tsp_nn(xy,dmat,popSize,showProg,showResult);
%
% Example:
%     n = 100;
%     phi = (sqrt(5)-1)/2;
%     theta = 2*pi*phi*(0:n-1);
%     rho = (1:n).^phi;
%     [x,y] = pol2cart(theta(:),rho(:));
%     xy = 10*([x y]-min([x;y]))/(max([x;y])-min([x;y]));
%     popSize = n;
%     showProg = 1;
%     showResult = 1;
%     a = meshgrid(1:n);
%     dmat = reshape(sqrt(sum((xy(a,:)-xy(a',:)).^2,2)),n,n);
%     [optRoute,minDist] = tsp_nn(xy,dmat,popSize,showProg,showResult);
%
% Example:
%     n = 50;
%     xyz = 10*rand(n,3);
%     popSize = n;
%     showProg = 1;
%     showResult = 1;
%     a = meshgrid(1:n);
%     dmat = reshape(sqrt(sum((xyz(a,:)-xyz(a',:)).^2,2)),n,n);
%     [optRoute,minDist] = tsp_nn(xyz,dmat,popSize,showProg,showResult);
%
% See also: tsp_ga, tspo_ga, tspof_ga, tspofs_ga, distmat
%
% Author: Joseph Kirk
% Email: jdkirk630@gmail.com
% Release: 1.3
% Release Date: 11/07/11
function varargout = tsp_nn(xy,dmat,popSize,showProg,showResult)

% Process Inputs and Initialize Defaults
nargs = 5;
for k = nargin:nargs-1
    switch k
        case 0
            xy = 10*rand(100,2);
        case 1
            N = size(xy,1);
            a = meshgrid(1:N);
            dmat = reshape(sqrt(sum((xy(a,:)-xy(a',:)).^2,2)),N,N);
        case 2
            N = size(xy,1);
            popSize = N;
        case 3
            showProg = 1;
        case 4
            showResult = 1;
        otherwise
    end
end

% Verify Inputs
[N,dims] = size(xy);
[nr,nc] = size(dmat);
if N ~= nr || N ~= nc
    error('Invalid XY or DMAT inputs!')
end
n = N;

% Sanity Checks
popSize = max(1,min(n,round(real(popSize(1)))));
showProg = logical(showProg(1));
showResult = logical(showResult(1));

% Initialize the Population
pop = zeros(popSize,n);

% Run the NN
distHistory = zeros(1,popSize);
if showProg
    pfig = figure('Name','TSP_NN | Current Solution','Numbertitle','off');
end
for p = 1:popSize
    d = 0;
    thisRte = zeros(1,n);
    visited = zeros(1,n);
    I = p;
    visited(I) = 1;
    thisRte(1) = I;
    for k = 2:n
        dists = dmat(I,:);
        dists(logical(visited)) = NaN;
        dMin = min(dists(~visited));
        J = find(dists == dMin,1);
        visited(J) = 1;
        thisRte(k) = J;
        d = d + dmat(I,J);
        I = J;
    end
    d = d + dmat(I,p);
    pop(p,:) = thisRte;
    distHistory(p) = d;

    if showProg
        % Plot the Current Route
        figure(pfig);
        rte = thisRte([1:n 1]);
        if dims > 2, plot3(xy(rte,1),xy(rte,2),xy(rte,3),'r.-');
        else plot(xy(rte,1),xy(rte,2),'r.-'); end
        title(sprintf('Total Distance = %1.4f',distHistory(p)));
    end
end

% Find the Minimum Distance Route
[minDist,index] = min(distHistory);
optRoute = pop(index,:);

if showResult
    % Plot the Best Route
    figure(pfig);
    rte = optRoute([1:n 1]);
    if dims > 2, plot3(xy(rte,1),xy(rte,2),xy(rte,3),'r.-');
    else plot(xy(rte,1),xy(rte,2),'r.-'); end
    title(sprintf('Total Distance = %1.4f',minDist));
    % Plots the NN Results
    figure('Name','TSP_NN | Results','Numbertitle','off');
    subplot(2,2,1);
    pclr = ~get(0,'DefaultAxesColor');
    if dims > 2, plot3(xy(:,1),xy(:,2),xy(:,3),'.','Color',pclr);
    else plot(xy(:,1),xy(:,2),'.','Color',pclr); end
    title('City Locations');
    subplot(2,2,2);
    imagesc(dmat(optRoute,optRoute));
    title('Distance Matrix');
    subplot(2,2,3);
    rte = optRoute([1:n 1]);
    if dims > 2, plot3(xy(rte,1),xy(rte,2),xy(rte,3),'r.-');
    else plot(xy(rte,1),xy(rte,2),'r.-'); end
    title(sprintf('Total Distance = %1.4f',minDist));
    subplot(2,2,4);
    plot(sort(distHistory,'descend'),'b','LineWidth',2);
    title('Distances');
    set(gca,'XLim',[0 popSize+1],'YLim',[0 1.1*max([1 distHistory])]);
end

% Return Outputs
if nargout
    varargout{1} = optRoute;
    varargout{2} = minDist;
    varargout{3} = pop;
end
end

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

end

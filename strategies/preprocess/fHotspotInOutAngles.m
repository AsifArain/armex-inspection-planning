
function [ hotspotSEQ,...
           InOutTSPAnglesAll,...
           InOutTSPAnglesMean,...
           FocusTSPAngles ] = fHotspotInOutAngles( hotspots,...
                                                   map_env,...
                                                   o,...
                                                   c,...
                                                   f,...
                                                   E,...
                                                   T,...
                                                   OBS,...
                                                   FoV,...
                                                   FoVfaces,...
                                                   robotStartCellInMap,...
                                                   SensorRange,...
                                                   para_,...
                                                   visualize )
%
%


% idices of centroid hotspots
% centroid_hotspot_cell_ind = find(map_centroid_hotspot)
centroid_hotspot_cell_ind = sub2ind(size(map_env),hotspots(:,1),hotspots(:,2));
centroid_hotspot_cell_ind = sort(centroid_hotspot_cell_ind);

% [centroid_hotspot_cell_r,centroid_hotspot_cell_c] = find(map_centroid_hotspot);

centroid_hotspot_conf = centroid_hotspot_cell_ind*f;

centroid_hotspot_conf_vector = zeros(1,c);
centroid_hotspot_conf_vector(centroid_hotspot_conf) = 1;

% startPosition_cell = [8,64]
% robotStartPosition_ind = sub2ind(size(map_env),paramPreprocess.robotStartPosition_cell(1),...
%                                                paramPreprocess.robotStartPosition_cell(2));
robotStartPosition_ind = sub2ind(size(map_env),robotStartCellInMap(1),...
                                               robotStartCellInMap(2));
robotStartPosition_conf = robotStartPosition_ind*f;
centroid_hotspot_conf_vector(robotStartPosition_conf) = 1;
% size(centroid_hotspot_conf_vector)

% removing conf on occupied cells
% --------------------------------

% obs-conf nums
confOBS = zeros(numel(OBS.ind)*f,1);
for i = 1:numel(OBS.ind)
    confOBS(((i-1)*f)+1:((i-1)*f)+f) = (((OBS.ind(i)-1)*f)+1):(((OBS.ind(i)-1)*f)+f);
end
centroid_hotspot_conf_vector(:,confOBS)  = [];
% size(centroid_hotspot_conf_vector)

% ---------------------------------------------------------------------------
% ---- TSP
% ---------------------------------------------------------------------------

% % TSP PARAMETERS
% % --------------
% paramTSP.TravCost    = 10; %ParamPreprocess.TravCost;
% paramTSP.ConfOBS     = 1;
% paramTSP.EoTinAstar  = 0;


[ tRoute,...
  TRoute,...
  tCost,...
  dmatTSP,...
  PredTSP,...
  tTSP ] = fTSP( o,...
                 FoVfaces,...
                 centroid_hotspot_conf_vector',...
                 E,...
                 T,...
                 map_env,...
                 OBS,...
                 para_ );
             

% ---- convert TSP solution conf into cell numbers
tRouteCELL = fix(tRoute/f)+(~~mod(tRoute,f));   % cell num for each conf num.
tRouteCELL(end) = [];

% find(tRouteCELL==robotStartPosition_ind);

tRouteCELL = circshift(tRouteCELL,-find(tRouteCELL==robotStartPosition_ind)+1);
% OR
% tRouteCELL = circshift(tRouteCELL,numel(tRouteCELL)-find(tRouteCELL==start_position_ind)+1)

hotspotSEQ = tRouteCELL(2:end);

% ---- convert TSP solution (with predecessor) conf into cell numbers
TRouteCELL = fix(TRoute/f)+(~~mod(TRoute,f));   % cell num for each conf num.
% -- remove repeated cells. 
% Note: this method will remove all duplicate entries of incoming and
% outgoing cells of a hotspot -- not good!
% [values,ia,ic]=unique(TRouteCELL);
% TRouteCELL = TRouteCELL(sort(ia),:)
% -- lets do some thing else
TRouteCELL = TRouteCELL([1;diff(TRouteCELL)]~=0);
% TRouteCELL(end) = [];

% -- shifting the vector to the start position
% find(TRouteCELL==robotStartPosition_ind)

TRouteCELL = circshift(TRouteCELL,numel(TRouteCELL)-...
    find(TRouteCELL==robotStartPosition_ind)+1);
incomingAnglesMean = inf(numel(hotspotSEQ),1);
outgoingAnglesMean = inf(numel(hotspotSEQ),1);

% inoutCellsAll = zeros(numel(hotspotSEQ)*10,2);
% inoutCellsAll = [];
inCellsAll  = [];
outCellsAll = [];

incomingAnglesAll = [];
outgoingAnglesAll = [];

% hotspotSEQ
% ANGLE_STEPS = 5;
% paramPreprocess.NumOfCellsForInOutAnglesTSP = 5;

for i = 1:numel(hotspotSEQ)
    
    %disp('----------- new hotspot ---------------')
    
    ind = find(TRouteCELL==hotspotSEQ(i));
    
    % indices for incoming angles.
    inInd = ind-1:-1:ind-para_.NumOfCellsForInOutAnglesTSP;
    if inInd(end) < 1
        inInd = [ind-1:-1:1,numel(TRouteCELL):-1:numel(TRouteCELL)-...
            (para_.NumOfCellsForInOutAnglesTSP-ind)];
    end
    
    % indices for outgoing angles.
    outInd = ind+1:+1:ind+para_.NumOfCellsForInOutAnglesTSP;
            
    if outInd(end)>numel(TRouteCELL)
        
        %outInd = [ind+1:1:numel(TRouteCELL), 1:1:numel(TRouteCELL)-ind+1]
        
        % -- #TODO temp fix
        outInd = [ind+1:1:numel(TRouteCELL), 1:1:numel(find(outInd>numel(TRouteCELL)))];
    end
    
    % incoming and outgoing cells
    %incomingCells = TRouteCELL(ind-1:-1:ind-5);
    %outgoingCells = TRouteCELL(ind+1:+1:ind+5);
    incomingCells = TRouteCELL(inInd);
    outgoingCells = TRouteCELL(outInd);
    
    [hotspot_r,hotspot_c]             = ind2sub(size(map_env),hotspotSEQ(i));
    [incomingCells_r,incomingCells_c] = ind2sub(size(map_env),incomingCells);
    [outgoingCells_r,outgoingCells_c] = ind2sub(size(map_env),outgoingCells);
    
    %inoutCellsAll = [inoutCellsAll; [incomingCells_r,incomingCells_c];...
    %    [outgoingCells_r,outgoingCells_c]];
    
    inCellsAll  = [inCellsAll;  [incomingCells_r,incomingCells_c]];
    outCellsAll = [outCellsAll; [outgoingCells_r,outgoingCells_c]];
    
    incomingAngles = atan2(((incomingCells_c)-(hotspot_c) ),...
        ((incomingCells_r)-(hotspot_r)))*(180/pi);
    incomingAngles = incomingAngles+(360*(incomingAngles<0));
    
    outgoingAngles = atan2(((outgoingCells_c)-(hotspot_c) ),...
        ((outgoingCells_r)-(hotspot_r)))*(180/pi);
    outgoingAngles = outgoingAngles+(360*(outgoingAngles<0));
    
    incomingAnglesAll = [incomingAnglesAll;incomingAngles];
    outgoingAnglesAll = [outgoingAnglesAll;outgoingAngles];
    
    % --- mean angles
    % Note: this is arthimetic mean in Eucledian-space, which is wrong.
    %incomingAngles_mean = mean(incomingAngles)
    %outgoingAngles_mean = mean(outgoingAngles)
    % --- here is the way in non-Eucledian space.
    % reference: [1] https://en.wikipedia.org/wiki/Mean_of_circular_quantities
    % reference: [2] http://catless.ncl.ac.uk/Risks/7.44.html#subj4
    % reference: [3] http://stackoverflow.com/questions/491738/how-do-you-calculate-the-average-of-a-set-of-angles
    incomingAngles_mean = rad2deg(atan2(deg2rad(sum(sind(incomingAngles))/numel(incomingAngles)),...
        deg2rad(sum(cosd(incomingAngles))/numel(incomingAngles))));
    incomingAngles_mean = incomingAngles_mean+(360*(incomingAngles_mean<0));
    
    outgoingAngles_mean = rad2deg(atan2(deg2rad(sum(sind(outgoingAngles))/numel(outgoingAngles)),...
        deg2rad(sum(cosd(outgoingAngles))/numel(outgoingAngles))));
    outgoingAngles_mean = outgoingAngles_mean+(360*(outgoingAngles_mean<0));

        
    incomingAnglesMean(i) = incomingAngles_mean;
    outgoingAnglesMean(i) = outgoingAngles_mean;
    
end

% incomingAnglesMean
% outgoingAnglesMean

% size(incomingAnglesAll)
% size(outgoingAnglesAll)
% 
% size(incomingAnglesMean)
% size(outgoingAnglesMean)

InOutTSPAnglesAll  = [incomingAnglesAll,outgoingAnglesAll];
InOutTSPAnglesMean = [incomingAnglesMean,outgoingAnglesMean];

FocusTSPAngles = rad2deg(atan2(deg2rad(sum(sind(InOutTSPAnglesMean),2)/size(InOutTSPAnglesMean,2)),...
        deg2rad(sum(cosd(InOutTSPAnglesMean),2)/size(InOutTSPAnglesMean,2))));
FocusTSPAngles = FocusTSPAngles+(360*(FocusTSPAngles<0));



if visualize == 1
V = [];

% publish parameteres: 
% ----------------------
ParamPublish.MarkVis    = 0;  % mark visibility.
ParamPublish.ConfSymbol = 1;  % draw sensing conf symbol.
ParamPublish.FoV        = 0;  % draw field of view.
ParamPublish.fillFoV    = 0;  % fill field of view with colors.
ParamPublish.Map        = 1;  % publish map.
ParamPublish.CellSize   = 10; % cell size.
ParamPublish.TravRoute  = 1;  % draw travelling route.
ParamPublish.PrintPDF   = 0;  % draw travelling route.
ParamPublish.PrintEPS   = 0;  % draw travelling route.
ParamPublish.PrintTIKZ  = 0;  % draw travelling route.

figure('name','global TSP');
[ tPublishTSP ] = fPublish_TSP( map_env,...
                                FoV,FoVfaces,...
                                V,...
                                centroid_hotspot_conf_vector',...
                                SensorRange,...
                                OBS,...
                                TRoute,...
                                ParamPublish );
hold on;
% for i = 1:size(inoutCellsAll,1)
%     plot(inoutCellsAll(i,1)-0.5,inoutCellsAll(i,2)-0.5,'og','MarkerFaceColor','g','MarkerSize',2);
% end
for i = 1:size(inCellsAll,1)
    plot(inCellsAll(i,1)-0.5,inCellsAll(i,2)-0.5,'og','MarkerFaceColor','g','MarkerSize',4);
end
for i = 1:size(outCellsAll,1)
    plot(outCellsAll(i,1)-0.5,outCellsAll(i,2)-0.5,'xr','MarkerFaceColor','r','MarkerSize',4);
end
end

end

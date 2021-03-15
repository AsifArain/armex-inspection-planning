function [ V,...
           oV,...
           constV,...
           redundantConfGlobal_ind,...
           confKept,...
           confRemoved,...
           cellKept,...
           cellRemoved,...
           FoVfaces,...
           map_candConfFusion,...
           confTheta,...
           confCrossAngles,...
           confCrossAngles_G,...
           infoGain_sppHotspots,...
           Z ] = fPreprocessFusion( FoVfaces,...
                                    hotspotSEQ,...
                                    map_env,...
                                    FoV,...
                                    SolutionLoc,...
                                    map_coverage,...
                                    SensorRange,...
                                    OBS,...
                                    alpha,...
                                    beta,...
                                    gamma,...
                                    NumberOfSelectedConf_sppHotspots,...
                                    paramPreprocess )





%% SENSING PARAMETERS AND CONSTANTS:
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
FoVfaces.num        = paramPreprocess.numFoVs;   % number of conf per cell.
n                   = numel(map_env);         % number of cells in the map.
f                   = FoVfaces.num;       % for simplification.
c                   = n*f;                % total conf.



%% RETRIEVE INFO FROM SPP-HOTSPOTS
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Lets retrieve objective values (G and G+U) of subproblems, reference
% distance points, maximum Ds and 1Vs.
disp('... retrieving data from SPP-Hotspots solutions')

% load preprocessing of first center of mass for initialization.
dataHotspots = load([SolutionLoc,'preprocess_coverage_crossangles_hotspot',num2str(1),'.mat'],...
    'oV','map_candConf','map_coverage');


% initialize.
selectedConfGlobal_vec          = zeros(size(dataHotspots.oV,2),1); % global conf vector
selectedConfHotspots_vec        = zeros(size(dataHotspots.oV,2),numel(hotspotSEQ)); % conf vector for each hotspot
selectedConfThetaHotspots_vec   = zeros(size(dataHotspots.oV,2),numel(hotspotSEQ)); % orientations of conf in each hotspot
referenceDistPointHotspots      = zeros(size(hotspotSEQ,1),2); % reference point used to calculate D of each hotspots [r,c]
sppHotspots_optVal              = zeros(size(hotspotSEQ)); % objective value (G+U) for each spp-hotspot subproblem
sppHotspots_optVal_G            = zeros(size(hotspotSEQ)); % objective value (G) for each spp-hotspot subproblem
maxDs                           = zeros(size(hotspotSEQ)); % maximum distance value of each hotspot to be used to normalize D vector
max1Vs                          = zeros(size(hotspotSEQ)); % maximum 1V value of each hotspot to be used to normalize 1V vector


% for each hotspot
for i = 1:numel(hotspotSEQ)
    
    % load preprocessing
    dataHotspots = load([SolutionLoc,'preprocess_coverage_crossangles_hotspot',num2str(i),'.mat'],...
        'oV','V','confKept','conf_theta','conf_crossAngles_G',...
        'refDistPoint','D');
    
    % load SPP solution for the center of mass
    selectedConf = load([SolutionLoc,'/spp_selectedconf_hotspot',num2str(i),'.dat']);
    
    selectedConf = round(selectedConf);
    
    % conf vector with all cells in the map.
    conf = zeros(size(dataHotspots.oV,2),1);
    conf(dataHotspots.confKept) = selectedConf;
           
    % update overall conf vector.
    selectedConfGlobal_vec = selectedConfGlobal_vec+conf;
    
    % update conf vector for all hotspots
    selectedConfHotspots_vec(:,i) = conf;
    
    % orientations of the selected conf in each spp-hotspots
    selectedConfThetaHotspots_vec(dataHotspots.confKept(find(selectedConf)),i) =...
        dataHotspots.conf_theta(find(selectedConf));
    
    % spp-hotspots solution/objective values
    C = selectedConf;
    G = alpha*dataHotspots.conf_crossAngles_G;
    
    Uc = gamma*((ones(size(dataHotspots.V,1),1)'*dataHotspots.V)/max((ones(size(dataHotspots.V,1),1)'*dataHotspots.V)));
    Ud = (1-gamma)*(1-(dataHotspots.D/max(dataHotspots.D)));
    U  = beta*(Uc'+Ud);

    sppHotspots_optVal(i) = C'*(G*C+U);
    %C'*(alpha*conf_crossAngles_G*C+beta*(  (gamma*((ones(size(V,1),1)'*V)/max((ones(size(V,1),1)'*V))))' +  (1-gamma)*(1-(D/max(D)))  ))
    
    sppHotspots_optVal_G(i) = C'*(G*C);
    
    % max distances in Ds of each hotpsot
    maxDs(i) = max(dataHotspots.D);
    
    % max visibilites in 1Vs of each hotpsot
    max1Vs(i) = max(sum(dataHotspots.V,1));
    
    % reference distance point [r,c] in each hotspot
    referenceDistPointHotspots(i,:) = dataHotspots.refDistPoint;
    
end

% sppHotspots_optVal
% referenceDistPointHotspots
% maxDs
% sppHotspots_optVal_G
% max1Vs

%% REDUNDANT CONFIGURATIONS
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

disp('... finding redundant configurations')


% --- global conf cells numbers
% -----------------------------------------------------------
% indices of global conf
selectedConfGlobal_ind = find(selectedConfGlobal_vec); 
% indices of global conf cells
selectedConfGlobalCell_ind =...
    fix(selectedConfGlobal_ind/f)+(~~mod(selectedConfGlobal_ind,f)); 
 %subscripts of global conf indices
[selectedConfGlobalCell_r,selectedConfGlobalCell_c] =...
    ind2sub(size(map_env),selectedConfGlobalCell_ind);


% --- plot global conf cells
% -----------------------------------------------------------
figure('name','all selected conf from subproblems'); imshow(map_env); hold on;
plot(selectedConfGlobalCell_c,selectedConfGlobalCell_r,'*b');


% --- distance between all selected conf
% -----------------------------------------------------------
dist_ConfGlobal = zeros(numel(selectedConfGlobalCell_ind));

for i = 1:numel(selectedConfGlobalCell_ind)
    for j = 1:numel(selectedConfGlobalCell_ind)
        dist_ConfGlobal(i,j) = sqrt(((selectedConfGlobalCell_c(j)-selectedConfGlobalCell_c(i))^2)+...            
                                    ((selectedConfGlobalCell_r(j)-selectedConfGlobalCell_r(i))^2));
    end
end
% distance to itself is inf
for i = 1:numel(selectedConfGlobalCell_ind)
    dist_ConfGlobal(i,i) = inf;
end


% --- check for redundant conf
% -----------------------------------------------------------
% min distance (in cells) to check redundant conf.
DistForRedConfs_cell = paramPreprocess.DistanceForRedundantConfs_m/paramPreprocess.MapRes;

% --- finding redundant conf from the distance matrix
minDist_vec = zeros(size(dist_ConfGlobal,1),1);

for i = 1:size(dist_ConfGlobal,1)
    [i,find(dist_ConfGlobal(i,:) <= DistForRedConfs_cell)]
    if any(dist_ConfGlobal(i,:) <= DistForRedConfs_cell)
        minDist_vec(i) = 1;
    end
end
minDist_vec


% --- finding sets/pairs of redundant confs.
% -----------------------------------------------------------
redundantConfPairs = zeros(numel(selectedConfGlobalCell_ind));
for i = 1:numel(selectedConfGlobalCell_ind)
    redundantConfPairs(i,i) = 1; % include ith conf vs ith conf
    redundantConfPairs(i,find(dist_ConfGlobal(i,:)<= DistForRedConfs_cell)) = 1;    
end
% remove duplicate rows
redundantConfPairs = unique(redundantConfPairs,'rows');

% QUESTIONABLE
while sum(sum(redundantConfPairs,1),2)~=size(redundantConfPairs,2)
    B = [];
    for i = 1:size(redundantConfPairs,2)
        B = [B;sum(redundantConfPairs(find(redundantConfPairs(:,i)),:),1)];        
    end
    C = unique(B,'rows');
    D = zeros(size(C));
    D(find(C)) = 1;
    redundantConfPairs = D;
end

% finding single conf rows (to be removed from the list)
toRemove_ind = zeros(size(redundantConfPairs,1),1);
for i = 1:size(redundantConfPairs,1)
    if numel(find(redundantConfPairs(i,:)))<=1
        toRemove_ind(i) = 1;
    end
end
% now remove them
redundantConfPairs(find(toRemove_ind),:) = [];
redundantConfPairs







% --- indices of (global) redundant conf.
% -----------------------------------------------------------
redundantConfGlobal_ind = selectedConfGlobal_ind(find(minDist_vec));


% --- indices of (global) redundant conf cell.
% -----------------------------------------------------------
redundantConfGlobalCell_ind = selectedConfGlobalCell_ind(find(minDist_vec));
[redundantConfGlobalCell_r,redundantConfGlobalCell_c] =...
    ind2sub(size(map_env),redundantConfGlobalCell_ind);


% --- plot (global) redundant conf cells
figure('name','redundant confs'); imshow(map_env); hold on;
plot(redundantConfGlobalCell_c,redundantConfGlobalCell_r,'*r');


% --- vector of (global) redundant conf.
% -----------------------------------------------------------
redundantConfGlobal_vec = zeros(size(dataHotspots.oV,2),1);
redundantConfGlobal_vec(redundantConfGlobal_ind) = 1;


% --- subject hotspots (vec) of redundant conf
% -----------------------------------------------------------
subjectHotspots_vec = zeros(size(hotspotSEQ));
for i = 1:size(selectedConfHotspots_vec,2)
    for j = 1:numel(redundantConfGlobal_ind)
        
        %find(find(confHotspots(:,i))==redundantConfGlobal_ind(j))
        
        % if any of selected conf for the hotspot is in the list of
        % redundant confs.           
        if find(find(selectedConfHotspots_vec(:,i))==redundantConfGlobal_ind(j))            
            subjectHotspots_vec(i) = 1;
        end
    end
end
subjectHotspots_vec


% --- subject hotspots (ind) of redundant conf
% -----------------------------------------------------------
subjectHotspots_ind = hotspotSEQ(find(subjectHotspots_vec))
[subjectHotspots_r,subjectHotspots_c] = ind2sub(size(map_env),subjectHotspots_ind);
subHotspotsCELL.ind = subjectHotspots_ind;
subHotspotsCELL.sub = [subjectHotspots_r,subjectHotspots_c];
% uncertainConfGlobal = confGlobal_ind(minDist_ConfGlobal_r)


% --- selected conf for the hotspots in question
% -----------------------------------------------------------
selectedConfSubjectHotspots_vec = selectedConfHotspots_vec(:,find(subjectHotspots_vec));

% NumberOfSelectedConf_sppHotspots = 3;
selectedConfSubjectHotspots_ind =...
    zeros(NumberOfSelectedConf_sppHotspots,size(selectedConfSubjectHotspots_vec,2));
for i = 1:size(selectedConfSubjectHotspots_vec,2)    
    selectedConfSubjectHotspots_ind(:,i) = find(selectedConfSubjectHotspots_vec(:,i));
end
selectedConfSubjectHotspots_ind


redundantConfGlobal_ind

redundantConfPairs

redundantConfSet = selectedConfGlobal_ind(find(redundantConfPairs(1,:)))


%% CONDENSED CONFIGURATIONS
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% condensed conf (bin) vector for each hotspot
condensedConfHotspots_vec = selectedConfHotspots_vec;
for i = 1:size(condensedConfHotspots_vec,2)    
    condensedConfHotspots_vec(find(selectedConfHotspots_vec(:,i).*redundantConfGlobal_vec),i) = 0;
    
    find(condensedConfHotspots_vec(:,i))
    
end

% condensed conf counts for each hotspot
condensedConfHotspots_num = zeros(1,size(condensedConfHotspots_vec,2));
for i = 1:numel(condensedConfHotspots_num)    
    find(condensedConfHotspots_vec(:,i))    
    condensedConfHotspots_num(i) = numel(find(condensedConfHotspots_vec(:,i)));
end
condensedConfHotspots_num

pause
%% CANDIDATE CONF MAP (INITIAL)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% map for the candidate sensing conf
map_candConfFusion_1st = zeros(size(map_env));
map_candConfFusion_1st(redundantConfGlobalCell_ind) = 1;

% plot
figure('name','cand nearest conf'); imshow(map_env); hold on;
for i = 1:size(map_candConfFusion_1st,1)
    for j = 1:size(map_candConfFusion_1st,2)        
        if map_candConfFusion_1st(i,j)
            plot(j,i,'xb','MarkerFaceColor','b');
        end
    end
end
plot(subjectHotspots_c,subjectHotspots_r,'or','MarkerFaceColor','r');

% grow area
radiusFreeCell_cell = 8 %paramPreprocess.radiusFreeCell_m/paramPreprocess.MapRes;

se = strel('disk',radiusFreeCell_cell,0); % outer limit
map_candConfFusion_1st = imdilate(map_candConfFusion_1st,se);

% plot
figure('name','cand nearest conf enlarged area'); imshow(map_env); hold on;
for i = 1:size(map_candConfFusion_1st,1)
    for j = 1:size(map_candConfFusion_1st,2)
        if map_candConfFusion_1st(i,j)
            plot(j,i,'xb','MarkerFaceColor','b');
        end
    end
end
plot(subjectHotspots_c,subjectHotspots_r,'or','MarkerFaceColor','r');



%% VISIBILITY MATRIX:
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

disp('V...')
[ V,oV,confKept,confRemoved,cellKept,cellRemoved,FoVfaces ] = ...
                fV( map_candConfFusion_1st,...
                    map_coverage,...
                    map_env,...
                    n,...
                    f,...
                    c,...
                    FoV,...
                    FoVfaces,...
                    SensorRange,...
                    OBS,...
                    hotspotSEQ,...
                    subjectHotspots_vec,...
                    paramPreprocess );

size(V)
size(oV)

cellKept

hotspotSEQ

subjectHotspots_num = find(subjectHotspots_vec)

hotspotSEQ(subjectHotspots_num)

find(ismember(cellKept,hotspotSEQ(subjectHotspots_num)))

cellKept(find(ismember(cellKept,hotspotSEQ(subjectHotspots_num))))

% V to be used for the contraint
constV = V(find(ismember(cellKept,hotspotSEQ(subjectHotspots_num))),:)

size(V)



%% UPDATED CANDIDATE CONF MAP :
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% cand conf that do not cover hotspot are removed from the second map.

% initialize final candidate conf map
map_candConfFusion = zeros(size(map_candConfFusion_1st));

% corresponding cell numbers for conf kept.
confCELL = fix(confKept/f)+(~~mod(confKept,f)); % cell num for each conf num.

% update the map
map_candConfFusion(confCELL) = 1;
% figure; imshow(map_candConf)
% numel(find(map_candConf))*4

% plot on environment map
figure('name','cand conf for fusion - final'); imshow(map_env); hold on;
% [candCell_r,candCell_c] = find(map_candConf);
for i = 1:size(map_candConfFusion,1)
    for j = 1:size(map_candConfFusion,2)        
        if map_candConfFusion(i,j)
            plot(j,i,'xb','MarkerFaceColor','b');
        end
    end
end
plot(subjectHotspots_c,subjectHotspots_r,'or','MarkerFaceColor','r');
numel(find(map_candConfFusion))

%% INFORMATION GAIN FOR THE CONDENSED CONFIGURATIONS.
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% subjectHotspots_num = find(subjectHotspots)

% get G, Uc and Ud of condensed conf

redundantConfGlobal_ind
selectedConfGlobal_ind

redundantConfGlobalCell_ind
selectedConfGlobalCell_ind


% initialize
condensedConfSubHotspot_infoGainG  = zeros(size(subjectHotspots_num));
condensedConfSubHotspot_infoGainUc = zeros(size(subjectHotspots_num));
condensedConfSubHotspot_infoGainUd = zeros(size(subjectHotspots_num));

% for each center of mass
subjectHotspots_num
for i = 1:numel(subjectHotspots_num)
    
    % condensed configurations
    condConf = find(condensedConfHotspots_vec(:,subjectHotspots_num(i)))
    
    % load preprocessing
    dataHotspots = load([SolutionLoc,'preprocess_coverage_crossangles_hotspot',num2str(subjectHotspots_num(i)),'.mat'],...
        'V','confKept','conf_crossAngles_G','D');
    
    % conf 
    confKeptCondensed = find(ismember(dataHotspots.confKept,condConf))
    
    C = zeros(size(dataHotspots.confKept,2),1);
    C(confKeptCondensed) = 1;
       
    % G
    G = alpha*dataHotspots.conf_crossAngles_G;
    condensedConfSubHotspot_infoGainG(i) = C'*(G*C);
    
    % Uc
    Uc = gamma*((ones(size(dataHotspots.V,1),1)'*dataHotspots.V)/max((ones(size(dataHotspots.V,1),1)'*dataHotspots.V)));
    
    condensedConfSubHotspot_infoGainUc(i) = sum(Uc(confKeptCondensed));
    
    % Ud
    Ud = (1-gamma)*(1-(dataHotspots.D/max(dataHotspots.D)));
    condensedConfSubHotspot_infoGainUd(i) = sum(Ud(confKeptCondensed));    
    
end



%% CROSS ANGLES:
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

disp('Cross angles...')
[ confTheta,...
  confCrossAngles,...
  confCrossAngles_G] =...
                 fCrossAngles( f,...
                               map_env,...
                               oV,...
                               confKept,...
                               redundantConfGlobal_ind,...
                               selectedConfSubjectHotspots_ind,....
                               condensedConfHotspots_vec,...
                               subHotspotsCELL,...
                               subjectHotspots_vec,...
                               selectedConfThetaHotspots_vec,...
                               hotspotSEQ,...
                               condensedConfHotspots_num,...
                               paramPreprocess );

                           
%% DISTANCE MATRIX:
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
disp('Distance matrix...')

% ---- D (distance) matrix is manhattan distance between cand conf position
% and reference position
confKeptCells_ind = fix(confKept/f)+(~~mod(confKept,f));
[confKeptCells_r,confKeptCells_c] = ind2sub(size(map_env),confKeptCells_ind);

% subjectHotspots_num = find(subjectHotspots)

D = zeros(numel(subjectHotspots_ind),numel(confKeptCells_ind));

for i = 1:numel(subjectHotspots_ind)
    
    for j = 1:numel(confKeptCells_ind)
        
        manhattan_distance = abs(confKeptCells_r(j)-referenceDistPointHotspots(subjectHotspots_num(i),1))+...
            abs(confKeptCells_c(j)-referenceDistPointHotspots(subjectHotspots_num(i),2));
        
        D(i,j) = manhattan_distance;
    end
end

                           
%% Z:
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% U = sum(V,1)
% U = beta*(U/max(U))
% 
% Z = alpha*confCrossAngles_G + [U;U]

sppHotspots_optVal_G
sppHotspots_optVal

condensedConfSubHotspot_infoGainG
condensedConfSubHotspot_infoGainUc
condensedConfSubHotspot_infoGainUd


% max 1Vs to narmalize 1Vs
norm1Vs = 1./max1Vs(subjectHotspots_num);

% max Ds to narmalize Ds
normDs = 1./maxDs(subjectHotspots_num);

% Uc = gamma * (ones(numel(subjectHotspotsCELL_ind),1)*(sum(V,1)/max(sum(V,1))))

Uc = gamma     * (norm1Vs*(sum(V,1))) + (condensedConfSubHotspot_infoGainUc*ones(1,size(D,2)));
Ud = (1-gamma) * (1-((normDs*ones(1,size(D,2))).*D)) + condensedConfSubHotspot_infoGainUd*ones(1,size(D,2));

% normD = [D(1,:)/maxDs(subjectHotspots_num(1));...
%          D(2,:)/maxDs(subjectHotspots_num(2))];
% Ud = (1-gamma) * (1-normD); %(D/max(D))

U = beta*(Uc+Ud);

% G = alpha*((2*confCrossAngles_G)+(2*0.9954));

G = alpha*(2*confCrossAngles_G) + condensedConfSubHotspot_infoGainG*ones(1,size(D,2));

Z = G+U;

% check the solution - 1
C = ismember(confKept,redundantConfGlobal_ind(1));
Z*C'
G*C'

% check the solution - 2
C = ismember(confKept,redundantConfGlobal_ind(2));
Z*C'
G*C'


% info gain required for the spp-fusion
infoGain_sppHotspots = sppHotspots_optVal(subjectHotspots_num)

% check the solution - XX
C = zeros(size(confKept));
C(5) = 1;
Z*C'
G*C'


end

function [confTheta,...
          confCrossAngles,...
          confCrossAngles_G] =...
                         fCrossAngles( f,...
                                       map_env,...
                                       oV,...
                                       confKept,...
                                       redundantConfSubHotspots_ind,...
                                       selectedConfSubjectHotspots_ind,....
                                       condensedConfHotspots_vec,...
                                       subHotspotsCELL,...
                                       subjectHotspots_vec,...
                                       selectedConfThetaHotspots_vec,...
                                       hotspotSEQ,...
                                       condensedConfHotspots_num,...
                                       paramPreprocess )

% conf_theta            := conf orientations [1,c]
% conf_crossAngles      := cross angles between conf [c,c]
% conf_crossAngles_G    := cross angles gain [c,c]

%{
center_of_mass = paramPreprocess.CoM; % reference cell
[center_of_mass_y,center_of_mass_x] = ind2sub(size(map_candConf),center_of_mass);

% ----- conf orientations
conf_theta = 1000*ones(1,size(oV,2)); % vector contains orientations of candidate conf.

for conf = 1:size(oV,2)
    conf_cell_num = fix(conf/f)+(~~mod(conf,f));   % cell num for each conf num.
    [conf_cell_y,conf_cell_x] = ind2sub(size(map_candConf),conf_cell_num); % row and col of cell num.
    
    if map_candConf(conf_cell_y,conf_cell_x)
        deltaX = conf_cell_x - center_of_mass_x;
        deltaY = conf_cell_y - center_of_mass_y;
        conf_theta(conf) = atan2(deltaY,deltaX);
    end
    
end
%}


% ----- conf orientations
% ---------------------------------------------------------------------------------
% disp('----- conf orientations old -------')
% 
% confTheta = 1000*ones(numel(subHotspotsCELL.ind),size(oV,2)); % matrix contains orientations of candidate conf.
% 
% for i = 1:numel(subHotspotsCELL.ind)
%     for j = 1:numel(confKept)
%         
%         confKeptCELL_ind = fix(confKept(j)/f)+(~~mod(confKept(j),f));   % cell num for each conf num.
%         [confKeptCELL_r,confKeptCELL_c] = ind2sub(size(map_env),confKeptCELL_ind); % row and col of cell num.
%         
%         deltaR = confKeptCELL_r - subHotspotsCELL.sub(i,1); %TODO make sure subHotspotsCELL row and col are correct
%         deltaC = confKeptCELL_c - subHotspotsCELL.sub(i,2);
%         
%         confTheta(i,confKept(j)) = atan2(deltaC,deltaR);
%         
%     end
% end   


% ----- conf orientations
% ---------------------------------------------------------------------------------
% Since, facing angle to a candidate conf is different for each hotspot.
disp('----- conf orientations -------')

confTheta = 1000*ones(numel(subHotspotsCELL.ind),numel(confKept)); % matrix contains orientations of candidate conf.

for i = 1:numel(subHotspotsCELL.ind)
    for j = 1:numel(confKept)
        
        confKeptCELL_ind = fix(confKept(j)/f)+(~~mod(confKept(j),f));   % cell num for each conf num.
        [confKeptCELL_r,confKeptCELL_c] = ind2sub(size(map_env),confKeptCELL_ind); % row and col of cell num.
        
        deltaR = confKeptCELL_r - subHotspotsCELL.sub(i,1); %TODO make sure subHotspotsCELL row and col are correct
        deltaC = confKeptCELL_c - subHotspotsCELL.sub(i,2);
        
        confTheta(i,j) = atan2(deltaC,deltaR);
        
    end
end

confTheta

% % ----- cross angles old

% disp('----- cross angles old ------')
% % ---------------------------------------------------------------------------------
% 
% % weakConfGlobal_ind
% % confSubjectHotspots_ind
% 
% % confHotspotsStrong
% 
% % 
% confHotspotsStrong_num = zeros(1,size(confHotspotsStrong,2));
% for i = 1:numel(confHotspotsStrong_num)
%     confHotspotsStrong_num(i) = numel(find(confHotspotsStrong(:,i)));
% end
% 
% confHotspotsStrong_num*subjectHotspots
% % confThetaHotspots
% 
% confCrossAngles = sparse(confHotspotsStrong_num*subjectHotspots,size(oV,2));
% 
% subjectHotspots_ind = find(subjectHotspots)
% 
% 
% % plot for debugging
% figure('name','cross angles debug'); imshow(map_env); hold on;
% 
% fixConfNum = 1;
% for i = 1:numel(subjectHotspots_ind) %1:numel(subjectHotspotsCELL_ind)*NumOfFixedConfPerSubjectHotspot
%     
%     fixedConf_ind = find(confHotspotsStrong(:,subjectHotspots_ind(i)))
%     
%     currentHotspot_ind = hotspotSEQ(subjectHotspots_ind(i))    
%     [currentHotspot_r,currentHotspot_c] = ind2sub(size(map_env),currentHotspot_ind);
%     currentHotspot_H = plot(currentHotspot_c,currentHotspot_r,'or','MarkerFaceColor','r');
%     
%     
%     for j = 1:confHotspotsStrong_num(subjectHotspots_ind(i)) %1:numel(confKept)
%         
%         fixConfAngle = confThetaHotspots(fixedConf_ind(j),subjectHotspots_ind(i))
%         
%         fixedConf_ind(j)
%         
%         currentFixedConfCELL_ind = fix(fixedConf_ind(j)/f)+(~~mod(fixedConf_ind(j),f));
%         [currentFixedConfCELL_r,currentFixedConfCELL_c] = ind2sub(size(map_env),currentFixedConfCELL_ind);
%         
%         currentFixedConfCELL_H = plot(currentFixedConfCELL_c,currentFixedConfCELL_r,'ob');
%         
%         for k = 1:numel(confKept)
%             
%             fixConfNum
%             confKept(k)
%             confCrossAngles(fixConfNum,confKept(k)) = angleAbsDiff(fixConfAngle,confTheta(i,confKept(k)));
%             %confCrossAngles(((i-1)* numel(subjectHotspots_ind))+j,confKept(k)) = 
%             
%             % -- for debugging            
%             %{
%             fixConfAngle_degree = rad2deg(fixConfAngle)
%             subConfAngle_degree = rad2deg(confTheta(i,confKept(k)))
%             cross_angles_degree = rad2deg(angleAbsDiff(fixConfAngle,confTheta(i,confKept(k))))
%             
%             subConfCell_ind = fix(confKept(k)/f)+(~~mod(confKept(k),f));
%             [subConfCell_r,subConfCell_c] = ind2sub(size(map_env),subConfCell_ind);
%             subConfCell_H = plot(subConfCell_c,subConfCell_r,'+b','MarkerFaceColor','b');
%             
%             pause;
%             delete(subConfCell_H);
%             clc;
%             %}
%             
%         end
%         
%         fixConfNum = fixConfNum+1;
%         delete(currentFixedConfCELL_H);
%     end
%     delete(currentHotspot_H);
% end
% 
% size(confCrossAngles)


% ----- cross angles
disp('----- cross angles ------')
% ---------------------------------------------------------------------------------


% % 
% condensedConfHotspots_num = zeros(1,size(condensedConfHotspots_vec,2));
% for i = 1:numel(condensedConfHotspots_num)
%     
%     find(condensedConfHotspots_vec(:,i))
%     
%     condensedConfHotspots_num(i) = numel(find(condensedConfHotspots_vec(:,i)));
% end




condensedConfHotspots_num

subjectHotspots_vec

condensedConfHotspots_num*subjectHotspots_vec



% confThetaHotspots

% confCrossAngles = sparse(confHotspotsStrong_num*subjectHotspots,numel(confKept));
confCrossAngles = 1000*ones(condensedConfHotspots_num*subjectHotspots_vec,numel(confKept));

size(confCrossAngles)


subjectHotspots_ind = find(subjectHotspots_vec)


% plot for debugging
figure('name','cross angles debug'); imshow(map_env); hold on;

fixConfNum = 1;


for i = 1:numel(subjectHotspots_ind) % for subject hotspots
    
    % indices of fixed conf.s
    fixedConf_ind = find(condensedConfHotspots_vec(:,subjectHotspots_ind(i)))
    
    % index of current hotspot
    currentHotspot_ind = hotspotSEQ(subjectHotspots_ind(i))    
    [currentHotspot_r,currentHotspot_c] = ind2sub(size(map_env),currentHotspot_ind);
    currentHotspot_H = plot(currentHotspot_c,currentHotspot_r,'or','MarkerFaceColor','r');
    
    
    for j = 1:condensedConfHotspots_num(subjectHotspots_ind(i)) % for number of fixed conf.s
        
        % orientation of fixed conf.
        fixConfAngle = selectedConfThetaHotspots_vec(fixedConf_ind(j),subjectHotspots_ind(i))
        
        fixedConf_ind(j)
        
        currentFixedConfCELL_ind = fix(fixedConf_ind(j)/f)+(~~mod(fixedConf_ind(j),f));
        [currentFixedConfCELL_r,currentFixedConfCELL_c] = ind2sub(size(map_env),currentFixedConfCELL_ind);
        
        currentFixedConfCELL_H = plot(currentFixedConfCELL_c,currentFixedConfCELL_r,'ob');
        
        for k = 1:numel(confKept)
            
            %fixConfNum
            %confKept(k)
            confCrossAngles(fixConfNum,k) = angleAbsDiff(fixConfAngle,confTheta(i,k));
            %confCrossAngles(((i-1)* numel(subjectHotspots_ind))+j,confKept(k)) = 
            
            % -- for debugging            
            %{
            fixConfAngle_degree = rad2deg(fixConfAngle)
            subConfAngle_degree = rad2deg(confTheta(i,k))
            cross_angles_degree = rad2deg(angleAbsDiff(fixConfAngle,confTheta(i,k)))
            
            subConfCell_ind = fix(confKept(k)/f)+(~~mod(confKept(k),f));
            [subConfCell_r,subConfCell_c] = ind2sub(size(map_env),subConfCell_ind);
            subConfCell_H = plot(subConfCell_c,subConfCell_r,'+b','MarkerFaceColor','b');
            
            pause;
            delete(subConfCell_H);
            clc;
            %}
            
        end
        
        fixConfNum = fixConfNum+1;
        delete(currentFixedConfCELL_H);
    end
    delete(currentHotspot_H);
end

size(confCrossAngles)



% % conf variables corresponding to candidate conf map only: -- OLD --
% % ------------------------------------------------------------------------------
% 
% % --- instead, use info from fV:
% % update conf orientation vector
% confTheta = confTheta(confKept);
% 
% % update conf cross angles matrix
% confCrossAngles = confCrossAngles(:,confKept);
% % confCrossAngles = confCrossAngles(confKept,:);
% confCrossAngles = full(confCrossAngles);

% ------ cross angles gain
disp('----- gain')
gain_sigma = deg2rad(10); % standard deviation
gain_meu   = deg2rad(60); % expected value or center
confCrossAngles_G1 = gaussmf(confCrossAngles,[gain_sigma,gain_meu]);


gain_sigma = deg2rad(10); % standard deviation
gain_meu   = deg2rad(120); % expected value or center
confCrossAngles_G2 = gaussmf(confCrossAngles,[gain_sigma,gain_meu]);

confCrossAngles_G = confCrossAngles_G1+confCrossAngles_G2;
 
size(confCrossAngles_G)

I = zeros(numel(subjectHotspots_ind),size(confCrossAngles_G,1))


condensedConfHotspots_num

confCrossAngles = sparse(condensedConfHotspots_num*subjectHotspots_vec,size(oV,2));

subjectHotspots_ind

confSubHotspotsStrong_num = condensedConfHotspots_num(subjectHotspots_ind)

for i = 1:numel(subjectHotspots_ind)
    
    %I(i:i+confSubHotspotsStrong_num(i)-1)    
    I(i,((i-1)*confSubHotspotsStrong_num(i))+1:((i-1)*confSubHotspotsStrong_num(i))+confSubHotspotsStrong_num(i)) = 1;
end
I

confCrossAngles_G = I*confCrossAngles_G

size(confCrossAngles_G)


pause

% ------ publish gussian gain function for test only
%{ 
angles = 0:1:180; % in degrees
gain_sigma = 10; % standard deviation
% gain-1
gain_meu   = 60; % expected value or center
guss_gain1 = gaussmf(angles,[gain_sigma,gain_meu]);
% gain-2
gain_meu   = 120;% expected value or center
guss_gain2 = gaussmf(angles,[gain_sigma,gain_meu]);
% final
guss_gain = guss_gain1+guss_gain2;
% plot
figure; plot(angles,guss_gain1,'-r')
figure; plot(angles,guss_gain2,'-g')
figure; plot(angles,guss_gain ,'-b')
%}


end




function [ T,E,eE,A,oA,D,P ] = fETAC( map_env,cnt,n,f,c,o,e,oO,eO,cell_length_m,OBS,startConf,confKept,paramPreprocess )
%fTt returns traversing time matrix Ct and predecessor matrix P (which is not
%used so far).
%

% CONSTANTS:
% ----------
[map_row,map_col] = size(map_env);

ID1oO = find(oO(1,:)); % ID 1
ID2oO = find(oO(2,:)); % ID 2
ID3oO = find(oO(3,:)); % ID 3
ID4oO = find(oO(4,:)); % ID 4
ID5oO = find(oO(5,:)); % ID 5
if ~all(size(ID1oO))
    ID1oO = 0;
end
if ~all(size(ID2oO))
    ID2oO = 0;
end
if ~all(size(ID3oO))
    ID3oO = 0;
end
if ~all(size(ID4oO))
    ID4oO = 0;
end
if ~all(size(ID5oO))
    ID5oO = 0;
end
ID1eO = find(eO(1,:)); % ID 1
ID2eO = find(eO(2,:)); % ID 2
ID3eO = find(eO(3,:)); % ID 3
ID4eO = find(eO(4,:)); % ID 4
ID5eO = find(eO(5,:)); % ID 5
if ~all(size(ID1eO))
    ID1eO = 0;
end
if ~all(size(ID2eO))
    ID2eO = 0;
end
if ~all(size(ID3eO))
    ID3eO = 0;
end
if ~all(size(ID4eO))
    ID4eO = 0;
end
if ~all(size(ID5eO))
    ID5eO = 0;
end

TstepS = paramPreprocess.TravCost*cell_length_m*1.000; % translational step - straight
TstepD = paramPreprocess.TravCost*cell_length_m*1.400; % Translational step - diagonal
Rstep1 = paramPreprocess.TravCost*cell_length_m*0.125; % rotational step - 1st deg
Rstep2 = paramPreprocess.TravCost*cell_length_m*0.250; % rotatioanl step - 2nd deg


%%                                    T
% --------------------------------------------------------------------------------------------------------------------------------
% Initialize T:
T = zeros(c,o);
size(T)
for i = 1:n
    if map_env(i)
        [rr,col] = ind2sub(size(map_env),i); % index number            
        
        if cnt == 4 && f~=8
            % -- List of all neighbors - straight
            N = [rr+0 col+1;...
                 rr-1 col+0;...
                 rr+0 col-1;...
                 rr+1 col+0];
        elseif cnt == 8 || f == 8
            N = [rr+0 col+1;...
                 rr-1 col+1;...
                 rr-1 col+0;...
                 rr-1 col-1;...
                 rr+0 col-1;...
                 rr+1 col-1;...
                 rr+1 col+0;...
                 rr+1 col+1];
        end
        
        for j = 1:f
            oCONFn = (i-1)*f+j; % conf number (from/outgoing)

            for k = 1:o
                l = (j*(f~=1))+(k*(f==1)); % selecting neighbor cell.
                   
                if (any(ID1oO==k)) || (any(ID2oO==k)) || (ID2eO==k)
                    
                    if ~(N(l,1) <= 0 || N(l,1) > map_row ||...
                           N(l,2) <= 0 || N(l,2) > map_col ||...
                           (~~ID3eO*(~mod(l,2))) ||...
                           map_env(N(l,1),N(l,2)) == 0)
                       
                       T(oCONFn,k) = TstepS;
                    else
                       T(oCONFn,k) = 1e300;
                    end
                end
                
                if ((any(ID2oO==k))*(mod(l,2)))||((any(ID2eO==k))*(~mod(l,2)))
                    
                    if ~(N(l,1) <= 0 || N(l,1) > map_row ||...
                           N(l,2) <= 0 || N(l,2) > map_col ||...
                           map_env(N(l,1),N(l,2)) == 0)
                       
                       T(oCONFn,k) = TstepD;
                    else
                       T(oCONFn,k) = 1e300;
                    end
                end
                
                if (any(ID5oO==k)||any(ID4oO==k)) && (k==2||k==3)
                    T(oCONFn,k) = (Rstep1*(any(ID4oO==k)))+(Rstep2*(any(ID5oO==k)));
                end
            end
        end
    end
end

oT = T; % keep T with obs
size(oT)

% OBSTACLE FREE:
% Removing obstacle indices:
for i = 1:numel(OBS.ind)
    T(((OBS.ind(i)-1)*f+1)-(f*(i-1)):((OBS.ind(i)-1)*f+f)-(f*(i-1)),:)  = [];
end

size(T)

%%                                    E
% --------------------------------------------------------------------------------------------------------------------------------

% Initialize E:
E = zeros(e,1);

for i = 1:n
    if map_env(i)
        [rr,col] = ind2sub(size(map_env),i); % index number            
        
        if cnt == 4 && f~=8
            % -- List of all neighbors - straight
            N = [rr+0 col+1;...
                 rr-1 col+0;...
                 rr+0 col-1;...
                 rr+1 col+0];
        elseif cnt == 8 || f == 8
            N = [rr+0 col+1;...
                 rr-1 col+1;...
                 rr-1 col+0;...
                 rr-1 col-1;...
                 rr+0 col-1;...
                 rr+1 col-1;...
                 rr+1 col+0;...
                 rr+1 col+1];
        end

        for j = 1:f            
            oCONFn = (i-1)*f+j; % conf number (from/outgoing)             

            for k = 1:o
                oEDGEn = (oCONFn-1)*o+k;   % edge num (from/outgoing)
                l = (j*(f~=1))+(k*(f==1)); % selecting neighbor cell.

                if (any(ID1oO==k)) || (any(ID2oO==k)) || (ID2eO==k)
                    
                    if ~(N(l,1) <= 0 || N(l,1) > map_row ||...
                           N(l,2) <= 0 || N(l,2) > map_col ||...
                           map_env(N(l,1),N(l,2)) == 0)
                       
                       Ni = sub2ind(size(map_env),N(l,1),N(l,2));
                       iCONFn = (Ni-1)*f+j;         % conf num (to/incoming)
                       iEDGEn = (iCONFn-1)*o+k;     % edge num (to/incoming)
                       E(oEDGEn) = iEDGEn;
                       
                    else
                       iCONFn = oCONFn+(((f/2)-((j/f>0.5)*f))*~mod(f,2));
                       iEDGEn = ((oEDGEn+(((o/2)-((k/o>0.5)*o))*~mod(o,2)))*(iCONFn==oCONFn))+...
                           (((iCONFn-1)*o+k)*(iCONFn~=oCONFn));
                       E(oEDGEn) = iEDGEn;
                    end                    
                end

                if (any(ID5oO==k)||any(ID4oO==k)) && (k==2)
                    % neighbor to the 2nd edge
                    iCONFn = oCONFn+((1-((~mod(j,f))*f))*~mod(f,2));  % incoming conf#                    
                    iEDGEn = (iCONFn-1)*o+k;
                    E(oEDGEn) = iEDGEn;
                end

                if (any(ID5oO==k)||any(ID4oO==k)) && (k==3)                    
                    % neighbor to the 3rd edge
                    iCONFn = oCONFn-1+((((mod(j,f)==1)*f))*~mod(f,2)); % incoming conf#
                    iEDGEn = (iCONFn-1)*o+k;
                    E(oEDGEn) = iEDGEn;                     
                end
            end
        end
    end
end


% Rearranging IDs of connected neighbours for E:
for i = numel(OBS.ind):-1:1
    E(E>(OBS.ind(i)*(f*o)))= E(E>(OBS.ind(i)*(f*o)))-(f*o);
end

% OBSTACLE FREE:
% Removing obstacle indices:
for i = 1:numel(OBS.ind)
    E(((OBS.ind(i)-1)*(f*o)+1)-((f*o)*(i-1)):((OBS.ind(i)-1)*(f*o)+(f*o))-((f*o)*(i-1)),:) = [];
end

eE = sub2ind([numel(find(map_env))*f o],fix(E/o)+(~~mod(E,o)),mod(E,o)+(~mod(E,o))*o);
% eE = sub2ind([c o],fix(E/o)+(~~mod(E,o)),mod(E,o)+(~mod(E,o))*o);

%%                                    A
% --------------------------------------------------------------------------------------------------------------------------------
A  = []; % initialize A.
oA = []; % initialize A with obstacles.

if paramPreprocess.A
    
    % Initialize A:
    A = inf(length(T(:,1))); % for all obstacle free conf
    for i = 1:length(T(:,1))
        a = find((T(i,:))<=TstepD);
        for j = 1:numel(a)        
            oEDGE = (i-1)*o+a(j);
            iEDGE = E(oEDGE);
            iCONF = fix(iEDGE/o)+(~~mod(iEDGE,o));
            A(i,iCONF) = T(i,a(j));
        end
    end
    oA = A; % A with obstacles.
    for i = 1:numel(OBS.ind)
        oA = [oA(:,1:(OBS.ind(i)-1)*f) inf(length(oA(:,1)),f) ...
            oA(:,(OBS.ind(i)-1)*f+1:end)];
        oA = [oA(1:(OBS.ind(i)-1)*f,:); inf(f,length(oA(1,:))); ...
            oA((OBS.ind(i)-1)*f+1:end,:)];
    end

end

%%                             Distance (D) and Predecessor (P) matrix
% --------------------------------------------------------------------------------------------------------------------------------
D = []; P = [];

if paramPreprocess.D

    % GENERATING MATRIX Ct and P:
    % ---------------------------
    %Ct = inf(c); % initiate Ct 
    %P  = inf(c); % initiate P

    originalA = A; % save A
    %[A.rp,A.ci,A.ai]=sparse_to_csr(obsA);
    [rp,ci,ai]=sparse_to_csr(oA);
    A = [];
    A.rp = rp; A.ci = ci; A.ai = ai;

    %for i = 1:c
    %    %[d, pred] = dijkstra(obsA,i);
    %    [d, pred] = dijkstra(A,i);
    %    Ct(i,:)   = d;
    %    P(i,:)    = pred;
    %end
    oA
    size(oA)
    [D,P] = dijkstra(oA,startConf);
    
    A = [];
    A = originalA; % recover original A.

    % Removing OBS
    % ------------
    %{
    for i = 1:numel(OBS.ind)    
        % Ct
        Ct(:,((OBS.ind(i)-1)*f+1)-(f*(i-1)):((OBS.ind(i)-1)*f+f)-(f*(i-1))) = []; % from
        Ct(((OBS.ind(i)-1)*f+1)-(f*(i-1)):((OBS.ind(i)-1)*f+f)-(f*(i-1)),:) = []; % to
        % P
        P(:,((OBS.ind(i)-1)*f+1)-(f*(i-1)):((OBS.ind(i)-1)*f+f)-(f*(i-1))) = []; % from
        P(((OBS.ind(i)-1)*f+1)-(f*(i-1)):((OBS.ind(i)-1)*f+f)-(f*(i-1)),:) = []; % to
    end
    %}
    
    % update distance and predecessor matrix for the candidate conf only.
    D = D(confKept);
    P = P(confKept);
    
    
end

% % Number of edges in the map
% nEDGES = numel(find(A~=inf & A~=0));


end



function [ V,oV,confKept,confRemoved,cellKept,cellRemoved,FoVfaces ] = ...
                fV( map_candConf,...
                    map_coverage,...
                    map_env,...
                    n,...
                    f,...
                    c,...
                    FoV,...
                    FoVfaces,...
                    SensorRange,...
                    OBS,...
                    hotspotSEQ,...
                    subjectHotspots,...
                    paramPreprocess )
%fV generates Visibility matrix V.

% --- ANGLES:
FoVfaces.lgt = [1/FoVfaces.num:1/FoVfaces.num:1-(1/FoVfaces.num) 0]'*360;
FoVfaces.alf = FoVfaces.lgt-(FoV/2);
FoVfaces.bay = FoVfaces.lgt+(FoV/2);
% compensate negative alf
FoVfaces.alf(FoVfaces.alf<0) = 360+FoVfaces.alf(FoVfaces.alf<0);
FoVfaces.ang = [FoVfaces.alf FoVfaces.bay];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SensorDomain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('.. Sensor Domain')

switch paramPreprocess.SenDom
    
    case 'ComputeNew'
        
        % For each cell of the map.
        sD = cell(1,f); % sensor domain (visible).
        oD = cell(1,f); % sensor's osbtacle domain (used to check the 
                        % visibility of "sD" for each cell in "oD"]

        % a single reference cell in the map.
        for i = 1
            
            % sensor position.
            sensor = [0,0];
            
            % For each sensing configuration.
            for j = 1:f

                % zawia := initial, middle and final (opening) angle.
                zawia = [FoVfaces.ang(j,1) FoVfaces.lgt(j) FoVfaces.ang(j,2)];

                % Finding sensor domain, i.e. visible cells for the
                % configuration. 
                [ visibleRange ] = fSenDomain( sensor,SensorRange,zawia );
                % visibleRange := visible range (row,col).

                % post-processing for sorting and removing duplicates
                visibleRange = unique(visibleRange,'rows'); % remove repeated numbers
                visibleRange = sortrows(visibleRange);      % sort rows
                
                % here we have jth sensor's domain.
                sD{j} = visibleRange;
                
                
                % since, to varify if a cell in sensor's domain is visible,
                % we need to check all the cells that effect the
                % visibility. Some of them are even not in the sensing
                % range, therefore, we need to compute obstacle domain that
                % contains the cells in the sensor domain and connecting.
                % "obsDomain" will be used to check visibility of the
                % sensor's domain.

                % obstacle domain
                [ obsRange ] = fSenDomainOBS( sensor,SensorRange,zawia );
                % obsRange := obstacle range (row,col).
                
                % post-processing for sorting and removing duplicates
                obsRange = unique(obsRange,'rows'); % remove repeated numbers
                obsRange = sortrows(obsRange);      % sort rows
                
                % here we have jth sensor's domain.
                oD{j} = obsRange;
                
            end
        end

        %--------------------------------------
        % EFFECTED VISIBILITY BY EACH OBSTACLE
        %--------------------------------------
        
        disp('.. Visibility vs Obstacles')
        VObs = cell(1,f); % Visibility vs Obstacle
        
        for i = 1
            sensor = [0,0];
            % For each configuration of a cell.
            for j = 1:f
                
                % let take sensor's mask
                sensorDom  = sD{j};
                % and obstacle domain
                obstacDom  = oD{j};
                
                VObs{j} = zeros(size(obstacDom,1),size(sensorDom,1)); % no more a square matrix
                % rows    := obstacle cell,
                % columns := visibility effect.
                
                %for k = 1:size(sensorDom,1)
                for k = 1:size(obstacDom,1)
                    
                    % kth cell is an obstacle
                    %obstacDomOBS(k,3)  = 0; 
                    obs = obstacDom(k,:);
                    
                    % visibility effect on all the cells
                    [ vDom ] = fVisSDvObs( sensor,sensorDom,obs );
                    
                    % update visibility against kth obstacle
                    VObs{j}(k,:) = vDom;
                   
                end
               
            end
        end
        
        % save the mask and visibility info to a file.
        save(['SenDom_Rng',num2str(SensorRange.cell),'FoV',num2str(FoV)],...
            'sD','oD','VObs');
        
    case 'UseEarlier'
        
        load(['SenDom_Rng',num2str(SensorRange.cell),'FoV',num2str(FoV)]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('.. Visibility Matrix')
V = sparse(n,c); % initializing
% V = zeros(n,c); % initializing

% For each cell of the map.
for i = 1:n
    if map_candConf(i) % if its not occupied
        [rr,cc] = ind2sub(size(map_candConf),i);
        sensor = [rr,cc];
        
        % For each configuration of a cell.
        for j = 1:f
            
            % sensorDom := sensor domain translation w.r.t sensor position.
            sensorDom = [sD{j}(:,1)+sensor(1) sD{j}(:,2)+sensor(2)];
            
            % obstacDom := obstalce domain translated w.r.t sensor's
            % position.
            obstacDom = [oD{j}(:,1)+sensor(1) oD{j}(:,2)+sensor(2)];
            
            % VOBS := Visibility vs Obstacles matrix for jth FoV.
            VOBS = VObs{j};
            
            % visibile cells 
            [ visDom ] = fVisibileDomain( map_coverage,map_env,sensorDom,obstacDom,VOBS );
            
            % Update Visibility matrix.
            V(visDom,((i-1)*f)+j) = 1;
        end
    end
end



% sensing config can not visualize own cell
for i = 1:numel(map_candConf)    
    V(i,(i-1)*f+1:(i-1)*f+f) = 0;
end

oV = V; % oV := visibility with obstacles

% V = sparse(V);

%{
% OBSTACLE FREE VISIBILITY:
% ------------------------------------------------------------------------------
% -- removing cells
V(OBS.ind,:) = [];

% -- removing conf.
confOBS = zeros(numel(OBS.ind)*f,1);
for i = 1:numel(OBS.ind)
    % removing sensing configurations.
    %V(:,((OBS.ind(i)-1)*f+1)-(f*(i-1)):((OBS.ind(i)-1)*f+f)-(f*(i-1))) = [];
    confOBS(((i-1)*f)+1:((i-1)*f)+f) = (((OBS.ind(i)-1)*f)+1):(((OBS.ind(i)-1)*f)+f);
end
V(:,confOBS) = [];
%}

% VISIBILITY MATRIX CORRESPOND TO THE HOTSPOT ONLY:
% ------------------------------------------------------------------------------

% ---- cells:

% cells covered by candidate conf.
% covered_ind = ( spotV*ones(size(spotV,2),1) ); % or
covered_ind = sum(oV,2);
covered_num = find(covered_ind);
covered_cell_vector = zeros(1,size(oV,1));
covered_cell_vector(covered_num) = 1;

% cells need to be covered -- according to the coverage map
coverage_ind = find(map_coverage);
coverage_cell_vector = zeros(1,size(oV,1));
% size(coverage_cell_vector)
coverage_cell_vector(coverage_ind) = 1;
% size(coverage_cell_vector)

% unoccupied cells
free_cell_vector = ones(1,size(oV,1));
free_cell_vector(OBS.ind) = 0;
find(free_cell_vector)

% size(coverage_cell_vector)

% pause
% --- selecting cells that are covered and need to be covered and unoccupied.
% size(covered_cell_vector)
% size(coverage_cell_vector)
% size(free_cell_vector)
cov_cell     = covered_cell_vector.*coverage_cell_vector.*free_cell_vector;
cellKept     = find(cov_cell);
cellRemoved  = find(cov_cell==0);
V            = V(cellKept,:);

% ---- configurations:

% --- selecting conf number corresponding to candidate conf map only
% cells for the candidate conf.
candConfCells = find(map_candConf); 
% conf numbers
% candConfNum = zeros(numel(map_candConf)*f,1);
candConfNum = zeros(1,size(oV,2));
for i = 1:numel(candConfCells)
    candConfNum( (candConfCells(i)-1)*f+1:(candConfCells(i)-1)*f+f ) = 1;
end

% --- conf that cover hotspot ---------------
% tmpV = oV;

subjectHotspots_ind = find(subjectHotspots)

numel(subjectHotspots_ind)

hotspotSEQ(subjectHotspots_ind)

tmpV = zeros(numel(subjectHotspots_ind),size(oV,2));

for i = 1:numel(subjectHotspots_ind)
    tmpV(i,:) = oV(hotspotSEQ(subjectHotspots_ind(i)),:);
end

sumV_hotspots = ones(1,numel(subjectHotspots_ind))*tmpV;
% conf_cover_all_hotspots_num = find(sumV_hotspots==numel(subjectHotspots_ind));
conf_cover_any_hotspots_num = find(sumV_hotspots>=1);

% tmpV = tmpV(hotspot,:);
% conf_hotspot_num = find(tmpV);
conf_hotspot = zeros(1,size(oV,2));
% conf_hotspot(conf_hotspot_num) = 1;
conf_hotspot(conf_cover_any_hotspots_num) = 1;
% ---------------------------------------------

% conf belong to unoccupied cells
conf_free = ones(1,size(oV,2));
for i = 1:numel(OBS.ind)
    conf_free( (OBS.ind(i)-1)*f+1:(OBS.ind(i)-1)*f+f ) = 0;
end

% --- keeping candidate conf in the list that are from candidate conf map
% AND cover the hotspot AND are placed over unoccupied cells 

size(candConfNum)
size(conf_hotspot)
size(conf_free)

conf_vector  = candConfNum.*conf_hotspot.*conf_free;
confKept     = find(conf_vector);
confRemoved  = find(conf_vector==0);
V            = V(:,confKept);
V            = full(V);


end

function [ visibleRange ] = fSenDomain( sensor,SensorRange,zawia )
%fSenDomain returns visible cells (domain/mask) for a given sensing
%configuration. 
%

% INTERPOLATING INTERMEDIATE ANGLES:
% -------------------------------------
DELTA_ANGLE = 1;
% --- Intermediate angles of visibile window. ---
% First half:
angl_1h = zawia(1):DELTA_ANGLE:zawia(2);
if zawia(1) > zawia(2) % if final angle is greater than 360 degrees.
    angl_1h = zawia(1):DELTA_ANGLE:360;
end
% Second half:
angl_2h = zawia(2):DELTA_ANGLE:zawia(3);
if zawia(2) > zawia(3) % if final angle is greater than 360 degrees.
    angl_2h = [zawia(2):DELTA_ANGLE:360 0:DELTA_ANGLE:zawia(3)];
end
angl = unique([angl_1h angl_2h]);


visibleRange      = [];        %
% subCellANGspecial = cell(1,4); % subject cells of special angles.

for i = 1:numel(angl)
    ALF = angl(i);
    % Find visible cells along one edge point.
    [ visCells ] = fRayTracSD(sensor,SensorRange,ALF);
    
    % Updating visible range.
    visibleRange = unique([visibleRange; visCells],'rows');   

end

end

function [ effective_visCells ] = fRayTracSD(sensor,SensorRange,ALF)
% fRayTracSD is ray tracing for sensor domain.

% CONSTANTS:
stepSize = 0.1; % step size of the ray length while sampling from origion to max. 

% initialize some high number of indices
effective_visCells = zeros(round(2*length(0.5:stepSize:SensorRange.cell)),2);
num = 0;
for d = [0.5:stepSize:SensorRange.cell SensorRange.cell]
    dr = 0.001*round((d*cosd(ALF))/0.001);  % row
    dc = 0.001*round((d*sind(ALF))/0.001);  % column
    
    vis_r = round(sensor(1)+round(dr)); % visible cell if validated
    vis_c = round(sensor(2)+round(dc)); % visible cell if validated
    
    if abs((dc-0.5)-fix(dc-0.5))>0.5 && abs((dr-0.5)-fix(dr-0.5))>0.5
        num = num+1;
        effective_visCells(num,:) = [vis_r vis_c];
    end
end
effective_visCells(num+1:end,:) = [];
effective_visCells = unique(effective_visCells,'rows');


end

function [ obsRange ] = fSenDomainOBS( sensor,SensorRange,zawia )
%fSenDomainOBS returns sensor's domain to check visibility against given
%obstacles. It contains more cells than the sensor's domain (mask).
%

% INTERPOLATING INTERMEDIATE ANGLES:
% ----------------------------------
DELTA_ANGLE = 1;
% --- Intermediate angles of visibile window. ---
% First half:
angl_1h = zawia(1):DELTA_ANGLE:zawia(2);
if zawia(1) > zawia(2) % if final angle is greater than 360 degrees.
    angl_1h = zawia(1):DELTA_ANGLE:360;
end
% Second half:
angl_2h = zawia(2):DELTA_ANGLE:zawia(3);
if zawia(2) > zawia(3) % if final angle is greater than 360 degrees.
    angl_2h = [zawia(2):DELTA_ANGLE:360 0:DELTA_ANGLE:zawia(3)];
end
angl = unique([angl_1h angl_2h]);

% initialize the obstacle range.
obsRange = []; 

% find the list of cells (for obstacles) in the sensing range.
for i = 1:numel(angl)
    
    % angle.
    ALF = angl(i);
    
    % Finding cells that are intersected along the ray.
    [ obsCells ] = fRayTracOD(sensor,SensorRange,ALF);
    
    % Updating the list of cells to be checked for obstacles.
    obsRange = unique([obsRange; obsCells],'rows');

end


end

function [ obsCells ] = fRayTracOD(sensor,SensorRange,ALF)
% fRayTracOD is ray tracing to find sensor domain for visibility check
% against listed obstacles

% point corresponding to the sensor range and ray angle.
subPoint_x = sensor(1) + (SensorRange.cell*cosd(ALF));
subPoint_y = sensor(2) + (SensorRange.cell*sind(ALF));

subPoint = [subPoint_x,subPoint_y];


% conditions.
cond_x1 = ALF>=0&&ALF<=45;
cond_x2 = ALF>=135&&ALF<=225;
cond_x3 = ALF>=315&&ALF<=360;

cond_y1 = ALF>=45&&ALF<=135;
cond_y2 = ALF>=225&&ALF<=315;

% sampling length.
delta = 0.01;
% Note: fine/coarse sampling changes the results for intersection points
% and hence changes the list of cells interseted. TODO: Find the optimal
% value for delta.


if cond_x1 || cond_x3
    
    % intermediate x and y points.
    x = [sensor(1):delta:subPoint(1) subPoint(1)];
    y = x*tand(ALF);
    
    % 4 decimal precision
    x = round(x*1e+4)/1e+4;
    y = round(y*1e+4)/1e+4;
    
    % corresponding cell numbers.
    cell_xy = round([x',y']);
    
    % cells that belong to the intersection point (x.5,y.5)
    ind_x  = find(abs(fix(x)-x)==0.5);
    ind_yy = abs(fix(y(ind_x))-y(ind_x))==0.5;
    ind_y  = ind_x(ind_yy);
    
    for i = 1:numel(ind_y)
        cell_xy = [cell_xy;...
                   ceil(x(ind_y(i))),  ceil(y(ind_y(i)))  ;...
                   ceil(x(ind_y(i))),  floor(y(ind_y(i))) ;...
                   floor(x(ind_y(i))), ceil(y(ind_y(i)))  ;...
                   floor(x(ind_y(i))), floor(y(ind_y(i)))];
    end    
    
elseif cond_x2
    
    % intermediate x and y points.    
    x = [sensor(1):-delta:subPoint(1) subPoint(1)];
    y = x*tand(ALF);
    
    % 4 decimal precision
    x = round(x*1e+4)/1e+4;
    y = round(y*1e+4)/1e+4;
    
    % corresponding cell numbers.
    cell_xy = round([x',y']);
    
    % cells that belong to the intersection point (x.5,y.5)
    ind_x  = find(abs(fix(x)-x)==0.5);
    ind_yy = find(abs(fix(y(ind_x))-y(ind_x))==0.5);
    ind_y  = ind_x(ind_yy);
    
    for i = 1:numel(ind_y)
        cell_xy = [cell_xy;...
                   ceil(x(ind_y(i))),  ceil(y(ind_y(i)))  ;...
                   ceil(x(ind_y(i))),  floor(y(ind_y(i))) ;...
                   floor(x(ind_y(i))), ceil(y(ind_y(i)))  ;...
                   floor(x(ind_y(i))), floor(y(ind_y(i)))];
    end
    
    
elseif cond_y1
    
    % intermediate x and y points.
    y = [sensor(2):delta:subPoint(2) subPoint(2)];
    x = y/tand(ALF);
    
    % 4 decimal precision
    x = round(x*1e+4)/1e+4;
    y = round(y*1e+4)/1e+4;
    
    % corresponding cell numbers.
    cell_xy = round([x',y']);
    
    % cells that belong to the intersection point (x.5,y.5)
    ind_x  = find(abs(fix(x)-x)==0.5);
    ind_yy = find(abs(fix(y(ind_x))-y(ind_x))==0.5);
    ind_y  = ind_x(ind_yy);
    
    for i = 1:numel(ind_y)
        cell_xy = [cell_xy;...
                   ceil(x(ind_y(i))),  ceil(y(ind_y(i)))  ;...
                   ceil(x(ind_y(i))),  floor(y(ind_y(i))) ;...
                   floor(x(ind_y(i))), ceil(y(ind_y(i)))  ;...
                   floor(x(ind_y(i))), floor(y(ind_y(i)))];
    end
    
    
elseif cond_y2
    
    % intermediate x and y points.
    y = [sensor(2):-delta:subPoint(2) subPoint(2)];
    x = y/tand(ALF);
    
    % 4 decimal precision
    x = round(x*1e+4)/1e+4;
    y = round(y*1e+4)/1e+4;
    
    % corresponding cell numbers.
    cell_xy = round([x',y']);
    
    % cells that belong to the intersection point (x.5,y.5)
    ind_x  = find(abs(fix(x)-x)==0.5);
    ind_yy = find(abs(fix(y(ind_x))-y(ind_x))==0.5);
    ind_y  = ind_x(ind_yy);
    
    for i = 1:numel(ind_y)
        cell_xy = [cell_xy;...
                   ceil(x(ind_y(i))),  ceil(y(ind_y(i)))  ;...
                   ceil(x(ind_y(i))),  floor(y(ind_y(i))) ;...
                   floor(x(ind_y(i))), ceil(y(ind_y(i)))  ;...
                   floor(x(ind_y(i))), floor(y(ind_y(i)))];
    end
    
end


 obsCells = unique(cell_xy,'rows');

end

function [ vDom ] = fVisSDvObs( sensor,sensorDom,obs )
%fVisibileArea returns visible cells for a given sensing configuration.
%

% initialize the list of visible cells.
vDom = zeros(size(sensorDom,1),1);

% for all cells in the sensor domain, check visibility against obs.
for i = 1:size(sensorDom,1)
    %i
    % cell to be validated for visibility
    subCell = [sensorDom(i,1) sensorDom(i,2)];

    % validatation for the visibility
    %[ cellValid ] = fCellValidSD( sensor,subCell,obstacDomOBS );
    [ validation ] = fCellValidSD( sensor,subCell,obs );

    % update visDom
    vDom(i) = validation;
    %pause
end

end


function [ validation ] = fCellValidSD( sensor,subCell,obs )
% 
% 

% angle to the center of the cell
ang = atan2(((subCell(2))-(sensor(2))),((subCell(1))-(sensor(1))))*(180/pi);
ang = ang+(360*(ang<0));

% conditions to project the line along (+/-x,+/-y).
cond_x1 = ang>=0&&ang<=45;
cond_x2 = ang>=135&&ang<=225;
cond_x3 = ang>=315&&ang<=360;

cond_y1 = ang>=45&&ang<=135;
cond_y2 = ang>=225&&ang<=315;

% sampling length.
delta = 0.01;
% Note: fine/coarse sampling changes the results for intersection points
% and hence changes the list of cells interseted. TODO: Find the optimal
% value for delta.

if cond_x1 || cond_x3
    
    % intermediate x and y points.
    x = sensor(1):delta:subCell(1);
    y = x*tand(ang);
    
    % 4 decimal precision
    x = round(x*1e+4)/1e+4;
    y = round(y*1e+4)/1e+4;
    
    % corresponding cell numbers.
    cell_xy = round([x',y']);
    
    % cells that belong to the intersection point (x.5,y.5)
    ind_x  = find(abs(fix(x)-x)==0.5);
    ind_yy = find(abs(fix(y(ind_x))-y(ind_x))==0.5);
    ind_y  = ind_x(ind_yy);
    
    for i = 1:numel(ind_y)
        cell_xy = [cell_xy;...
                   ceil(x(ind_y(i))),  ceil(y(ind_y(i)))  ;...
                   ceil(x(ind_y(i))),  floor(y(ind_y(i))) ;...
                   floor(x(ind_y(i))), ceil(y(ind_y(i)))  ;...
                   floor(x(ind_y(i))), floor(y(ind_y(i)))];
    end    
    
elseif cond_x2
    
    % intermediate x and y points.
    x = sensor(1):-delta:subCell(1);
    y = x*tand(ang);
    
    % 4 decimal precision
    x = round(x*1e+4)/1e+4;
    y = round(y*1e+4)/1e+4;
    
    % corresponding cell numbers.
    cell_xy = round([x',y']);
    
    % cells that belong to the intersection point (x.5,y.5)
    ind_x  = find(abs(fix(x)-x)==0.5);
    ind_yy = find(abs(fix(y(ind_x))-y(ind_x))==0.5);
    ind_y  = ind_x(ind_yy);
    
    for i = 1:numel(ind_y)
        cell_xy = [cell_xy;...
                   ceil(x(ind_y(i))),  ceil(y(ind_y(i)))  ;...
                   ceil(x(ind_y(i))),  floor(y(ind_y(i))) ;...
                   floor(x(ind_y(i))), ceil(y(ind_y(i)))  ;...
                   floor(x(ind_y(i))), floor(y(ind_y(i)))];
    end
    
    
elseif cond_y1
    
    % intermediate x and y points.
    y = sensor(2):delta:subCell(2);
    x = y/tand(ang);
    
    % 4 decimal precision
    x = round(x*1e+4)/1e+4;
    y = round(y*1e+4)/1e+4;
    
    % corresponding cell numbers.
    cell_xy = round([x',y']);
    
    % cells that belong to the intersection point (x.5,y.5)
    ind_x  = find(abs(fix(x)-x)==0.5);
    ind_yy = find(abs(fix(y(ind_x))-y(ind_x))==0.5);
    ind_y  = ind_x(ind_yy);
    
    for i = 1:numel(ind_y)
        cell_xy = [cell_xy;...
                   ceil(x(ind_y(i))),  ceil(y(ind_y(i)))  ;...
                   ceil(x(ind_y(i))),  floor(y(ind_y(i))) ;...
                   floor(x(ind_y(i))), ceil(y(ind_y(i)))  ;...
                   floor(x(ind_y(i))), floor(y(ind_y(i)))];
    end
    
    
elseif cond_y2
    
    % intermediate x and y points.
    y = sensor(2):-delta:subCell(2);
    x = y/tand(ang);
    
    % 4 decimal precision
    x = round(x*1e+4)/1e+4;
    y = round(y*1e+4)/1e+4;
    
    % corresponding cell numbers.
    cell_xy = round([x',y']);
    
    % cells that belong to the intersection point (x.5,y.5)
    ind_x  = find(abs(fix(x)-x)==0.5);
    ind_yy = find(abs(fix(y(ind_x))-y(ind_x))==0.5);
    ind_y  = ind_x(ind_yy);
    
    for i = 1:numel(ind_y)
        cell_xy = [cell_xy;...
                   ceil(x(ind_y(i))),  ceil(y(ind_y(i)))  ;...
                   ceil(x(ind_y(i))),  floor(y(ind_y(i))) ;...
                   floor(x(ind_y(i))), ceil(y(ind_y(i)))  ;...
                   floor(x(ind_y(i))), floor(y(ind_y(i)))];
    end
    
end

% remove repeated cells
cell_xy = unique(cell_xy,'rows');

% validation for the visibility
ind  = find(cell_xy(:,1)==obs(1)&cell_xy(:,2)==obs(2));
% if size(ind,1) > 0
%     validation = 0;
% end
validation = (size(ind,1) < 1);

end

function [ visDom ] = fVisibileDomain( map_coverage,map_env,sensorDom,obstacDom,VOBS )
%fVisibileDomain returns visible cells for a given sensing configuration.
%It takes the mask (sensor domain), finds the non-visible cells due to
%obstacles in the map and remove them from the mask.
%
% INPUTS:
% ------
% map       := [x,y] map
% sensorDom := [n,2] mask for the sensing configuration in question.
% obstacDom := [m,2] mask for the sensing configuration in question (used
%                    to check visibility against obstacles).
% VOBS      := [m,n] effect on visibility (columns) against each obstacle 
%                    (rows)
% OUTPUTS:
% ------
% visDom    := [k] visibile cells for the configuration in question 
%              (subscripts)

% rows and columns for the map
[map_row,map_col] = size(map_coverage);

% temporary domain/masks to handle intermediate computation quickly.
visDomX = sensorDom;
obsDomX = obstacDom;

% preRemoveList:= cells outside the map (sensor domain).
% along x-axis:
   i = find(visDomX(:,1)<=0);      % less than row number
   j = find(visDomX(:,1)>map_row); % greater than row number
   preRemoveList = [i;j];
% along y-axis:
   i = find(visDomX(:,2)<=0);      % less than col number
   j = find(visDomX(:,2)>map_col); % greater than col number
   preRemoveList = [preRemoveList;i;j];

% remove preRemoveList from temporary visible domain/mask.
% visDomX(preRemoveList,:) = [];

% preRemoveListOD:= cells outside the map (obstacle domain).
% along x-axis:
   i = find(obsDomX(:,1)<=0);      % less than row number
   j = find(obsDomX(:,1)>map_row); % greater than row number
   preRemoveListOD = [i;j];
% along y-axis:
   i = find(obsDomX(:,2)<=0);      % less than col number
   j = find(obsDomX(:,2)>map_col); % greater than col number
   preRemoveListOD = [preRemoveListOD;i;j];

% remove preRemoveList from temporary visible domain/mask.
obsDomX(preRemoveListOD,:) = [];


% convert visDomX from ind to sub.
% visDomX_ind = sub2ind(size(map),visDomX(:,1),visDomX(:,2));
obsDomX_ind = sub2ind(size(map_coverage),obsDomX(:,1),obsDomX(:,2));

% find the list of obstacles in visDomX;
% obs_list = find(map(visDomX_ind)==0);
obs_listOD = find(map_env(obsDomX_ind)==0);

% VOBSx:= VOBS with outside cells are removed.
% VOBSx = VOBS; 
% VOBSx(preRemoveList,:) = [];

VOBSx = VOBS; 
VOBSx(preRemoveListOD,:) = [];

% non-visible cells due to obstacles (logical values).
% nVisCells = (VOBSx(obs_list',:)==0);
nVisCells = (VOBSx(obs_listOD',:)==0);

% convert logical values to a list.
list = find(nVisCells'==1);

% translate the list to the column numbers of VOBS, which is remove list.
[removeList,~] = ind2sub(size(VOBSx'),list);

% update the RemoveList with preRemoveList and removeList.
RemoveList = [preRemoveList;removeList];

% finally, we know non-visible cells in the domain/mask.
visDom = sensorDom;
visDom(RemoveList,:) = [];

% convert indices into subscripts and sort them.
if visDom  
    visDom = sort(sub2ind(size(map_coverage),visDom(:,1),visDom(:,2)));
end

end

% ---------------------------------------------------- End of the document -------------------------------------------------------


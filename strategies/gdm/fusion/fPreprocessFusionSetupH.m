function [ selectedConf_sppH,...
           selectedConfGlobal_sppH,...
           selectedConfGlobalCell_sppH,...
           selectedConfTheta_sppH,...           
           infoGain_sppH,...
           infoGainG_sppH,...
           infoGainU_sppH,...
           maxDs_sppH,...
           max1Vs_sppH,...
           referenceDistPoint_sppH,...
           distMatrix_selectedConfGlobal_sppH,...
           DistForRedConfs_cell,...
           minDistSelectedConfGlobal,...           
           redundantConfGlobal_sppH,...
           redundantConfGlobalCell_sppH,...
           redundantConfSets_sppH,...
           redundantConfPairsAll_sppH,...
           dist_redundantConfPairsAll_sppH,...
           subH_RedundantConfPairsAll_sppH] = fPreprocessFusionSetupH( FoVfaces,...
                                                                      hotspotSEQ,...
                                                                      map_env,...
                                                                      alpha,...
                                                                      beta,...
                                                                      dir_,...
                                                                      para_,...
                                                                      visualize )
%fPreprocessFusionSetup identify the redundant pairs of selected conf from
%SPP-HOTSPOTS.



%% SENSING PARAMETERS AND CONSTANTS:
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
FoVfaces.num        = para_.numFoVs;   % number of conf per cell.
f                   = FoVfaces.num;       % for simplification.


%% HOTSPOTS

hotspotSEQ;
Hotspots.ind = hotspotSEQ;
[r,c] = ind2sub(size(map_env),Hotspots.ind);
Hotspots.sub = [r,c];


%% RETRIEVE INFO FROM SPP-HOTSPOTS
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Lets retrieve objective values (G+U, G, and U) of subproblems, reference
% distance points, maximum Ds and 1Vs.
disp('... retrieving data from SPP-Hotspots solutions')

% load preprocessing of first center of mass for initialization.
priorH = load([dir_.TomographyLocalSol,'preprocess_XVD_hotspot',num2str(1),'.mat'],...
    'oV','map_candConf','map_coverage');

% initialize.
selectedConfGlobal_sppH.vec  = zeros(size(priorH.oV,2),1); % global conf vector
selectedConf_sppH.vec        = zeros(size(priorH.oV,2),numel(hotspotSEQ)); % conf vector for each hotspot
selectedConfTheta_sppH.vec   = zeros(size(priorH.oV,2),numel(hotspotSEQ)); % orientations of conf in each hotspot
referenceDistPoint_sppH      = zeros(size(hotspotSEQ,1),2); % reference point used to calculate D of each hotspots [r,c]
infoGain_sppH                = zeros(size(hotspotSEQ)); % objective value (G+U) for each spp-hotspot subproblem
infoGainG_sppH               = zeros(size(hotspotSEQ)); % objective value (G) for each spp-hotspot subproblem
infoGainU_sppH               = zeros(size(hotspotSEQ)); % objective value (U) for each spp-hotspot subproblem
maxDs_sppH                   = zeros(size(hotspotSEQ)); % maximum distance value of each hotspot to be used to normalize D vector
max1Vs_sppH                  = zeros(size(hotspotSEQ)); % maximum 1V value of each hotspot to be used to normalize 1V vector

if visualize == 1
figure('name','redundant'); 
imshow(map_env','InitialMagnification',700); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1]);
end

col_lines = lines(numel(hotspotSEQ));

% for each hotspot
for i = 1:numel(hotspotSEQ)
    
    % load preprocessing
    priorH = load([dir_.TomographyLocalSol,'preprocess_XVD_hotspot',num2str(i),'.mat'],...
        'oV','V','confKept','conf_theta','conf_crossAngles_G','refDistPoint','D');
    
    % load SPP solution for the center of mass       
    selectedConf = load([dir_.TomographyLocalSol,'selectedConf_hotspot',num2str(i),'.dat']);
    selectedConf = round(selectedConf);
    
    % conf vector with all cells in the map.
    conf = zeros(size(priorH.oV,2),1);
    conf(priorH.confKept) = selectedConf;
    
    % update overall conf vector.
    selectedConfGlobal_sppH.vec = selectedConfGlobal_sppH.vec+conf;
    
    % update conf vector for all hotspots
    selectedConf_sppH.vec(:,i) = conf;
    
    % orientations of the selected conf in each spp-hotspots
    selectedConfTheta_sppH.vec(priorH.confKept(find(selectedConf)),i) =...
        priorH.conf_theta(find(selectedConf));
    
    % debug
    %{    
    conf_debug_ind = dataHotspots.confKept(find(selectedConf))
    confcell_debug_ind = fix(conf_debug_ind/f)+(~~mod(conf_debug_ind,f)); 
    [confcell_debug_r,confcell_debug_c] = ind2sub(size(map_env),confcell_debug_ind)
    H_conf_debug = plot(confcell_debug_c,confcell_debug_r,'xr')
    hot_debug_ind = hotspotSEQ(i)
    [hot_debug_r,hot_debug_c] = ind2sub(size(map_env),hot_debug_ind)
    H_hotspot = plot(hot_debug_c,hot_debug_r,'ob')  
    rad2deg(dataHotspots.conf_theta(find(selectedConf)))
    pause;
    delete(H_conf_debug)
    delete(H_hotspot)
    %}
    
    
    % debug
    if visualize == 1
    conf_debug_ind = priorH.confKept(find(selectedConf));
    confcell_debug_ind = fix(conf_debug_ind/f)+(~~mod(conf_debug_ind,f)); 
        
    [confcell_debug_r,confcell_debug_c] = ind2sub(size(map_env),confcell_debug_ind);
    plot(confcell_debug_r,confcell_debug_c,'x','Color',col_lines(i,:));
    
    hot_debug_ind = hotspotSEQ(i);
    [hot_debug_r,hot_debug_c] = ind2sub(size(map_env),hot_debug_ind);
    plot(hot_debug_r,hot_debug_c,'o','Color',col_lines(i,:));
    end
    %}
    
    
    
    % spp-hotspots solution/objective values
    C = selectedConf;
    G = alpha*priorH.conf_crossAngles_G;
    
    Uc = beta*((ones(size(priorH.V,1),1)'*priorH.V)/max((ones(size(priorH.V,1),1)'*priorH.V)));
    Ud = (1-beta)*(1-(priorH.D/max(priorH.D)));
    U  = (1-alpha)*(Uc'+Ud);
    
    infoGain_sppH(i) = C'*(G*C+U);
    %C'*(alpha*conf_crossAngles_G*C+beta*(  (gamma*((ones(size(V,1),1)'*V)/max((ones(size(V,1),1)'*V))))' +  (1-gamma)*(1-(D/max(D)))  ))
    
    infoGainG_sppH(i) = C'*(G*C);
    
    infoGainU_sppH(i) = C'*U;
    
    % max distances in Ds of each hotpsot
    maxDs_sppH(i) = max(priorH.D);
    
    % max visibilites in 1Vs of each hotpsot
    max1Vs_sppH(i) = max(sum(priorH.V,1));
    
    % reference distance point [r,c] in each hotspot
    referenceDistPoint_sppH(i,:) = priorH.refDistPoint;
    
end



% selected conf indices
% -----------------------------------------------------------
selectedConf_sppH.ind = zeros(para_.NumberOfConf_SPPH,size(selectedConf_sppH.vec,2));
for i = 1:size(selectedConf_sppH.vec,2)
    selectedConf_sppH.ind(:,i) = find(selectedConf_sppH.vec(:,i));
end

% selected golbal conf indices
% -----------------------------------------------------------
selectedConfGlobal_sppH.ind = find(selectedConfGlobal_sppH.vec); 

% selected global conf cells (indices/subscripts)
% -----------------------------------------------------------
% indices of global conf cells
selectedConfGlobalCell_sppH.ind =...
    fix(selectedConfGlobal_sppH.ind/f)+(~~mod(selectedConfGlobal_sppH.ind,f)); 
 %subscripts of global conf indices
[r,c] = ind2sub(size(map_env),selectedConfGlobalCell_sppH.ind);
selectedConfGlobalCell_sppH.sub = [r,c];

% --- plot global conf cells
% -----------------------------------------------------------
% figure('name','all selected conf from subproblems'); imshow(map_env); hold on;
% plot(selectedConfGlobalCell_sppH.sub(:,2),selectedConfGlobalCell_sppH.sub(:,1),'*b');

%% REDUNDANT CONFIGURATIONS
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
disp('... finding redundant configurations')

% % --- Binomial coefficient or all combinations of hotspots
% hotspot_num = 1:numel(hotspotSEQ);
% % [p,q] = meshgrid(hotspot_num,hotspot_num);
% % hotspotsPairs = [p(:) q(:)];
% % a = unique(sort(hotspotsPairs,2),'rows')
% hotspotsAllPairs = nchoosek(hotspot_num,2);


% distance between all selected conf
% -----------------------------------------------------------
distMatrix_selectedConfGlobal_sppH = zeros(numel(selectedConfGlobalCell_sppH.ind));
for i = 1:numel(selectedConfGlobalCell_sppH.ind)
    for j = 1:numel(selectedConfGlobalCell_sppH.ind)
        distMatrix_selectedConfGlobal_sppH(i,j) =...
            sqrt(((selectedConfGlobalCell_sppH.sub(j,2)-selectedConfGlobalCell_sppH.sub(i,2))^2)+...
                 ((selectedConfGlobalCell_sppH.sub(j,1)-selectedConfGlobalCell_sppH.sub(i,1))^2));
    end
end
% distance to itself is inf
for i = 1:numel(selectedConfGlobalCell_sppH.ind)
    distMatrix_selectedConfGlobal_sppH(i,i) = inf;
end


% check for redundant conf
% -----------------------------------------------------------
% min distance (in cells) to check redundant conf.
DistForRedConfs_cell = para_.DistanceForRedundantConfs_m/para_.MapRes;
% --- finding redundant conf from the distance matrix
minDistSelectedConfGlobal.vec = zeros(size(distMatrix_selectedConfGlobal_sppH,1),1);

for i = 1:size(distMatrix_selectedConfGlobal_sppH,1)
    %[i,find(distMatrix_selectedConfGlobal_sppH(i,:) <= DistForRedConfs_cell)]
    if any(distMatrix_selectedConfGlobal_sppH(i,:) <= DistForRedConfs_cell)
        minDistSelectedConfGlobal.vec(i) = 1;
    end
end
% minDistSelectedConfGlobal.vec


% indices/vector of global redundant conf.
% ------------------------------------------------------------------------------------
% indices
redundantConfGlobal_sppH.ind = selectedConfGlobal_sppH.ind(find(minDistSelectedConfGlobal.vec));
% vec
redundantConfGlobal_sppH.vec = zeros(size(priorH.oV,2),1);
redundantConfGlobal_sppH.vec(redundantConfGlobal_sppH.ind) = 1;


% indices/subscripts of global redundant conf cells.
% ------------------------------------------------------------------------------------
redundantConfGlobalCell_sppH.ind = selectedConfGlobalCell_sppH.ind(find(minDistSelectedConfGlobal.vec));
[r,c] = ind2sub(size(map_env),redundantConfGlobalCell_sppH.ind);
redundantConfGlobalCell_sppH.sub = [r,c];

% sets of redundant confs.
% -----------------------------------------------------------
redundantConfSets_sppH = zeros(numel(selectedConfGlobalCell_sppH.ind));
for i = 1:numel(selectedConfGlobalCell_sppH.ind)
    redundantConfSets_sppH(i,i) = 1; % include ith conf vs ith conf
    redundantConfSets_sppH(i,find(distMatrix_selectedConfGlobal_sppH(i,:)<= DistForRedConfs_cell)) = 1;    
end
% remove duplicate rows
redundantConfSets_sppH = unique(redundantConfSets_sppH,'rows');

% QUESTIONABLE
while sum(sum(redundantConfSets_sppH,1),2)~=size(redundantConfSets_sppH,2)
    B = [];
    for i = 1:size(redundantConfSets_sppH,2)
        B = [B;sum(redundantConfSets_sppH(find(redundantConfSets_sppH(:,i)),:),1)];        
    end
    C = unique(B,'rows');
    D = zeros(size(C));
    D(find(C)) = 1;
    redundantConfSets_sppH = D;
end

% finding single conf rows (to be removed from the list)
toRemove_ind = zeros(size(redundantConfSets_sppH,1),1);
for i = 1:size(redundantConfSets_sppH,1)
    if numel(find(redundantConfSets_sppH(i,:)))<=1
        toRemove_ind(i) = 1;
    end
end
% now remove them
redundantConfSets_sppH(find(toRemove_ind),:) = [];
% redundantConfSets_sppH


% all pairs of redundant conf
% ------------------------------------------------------------------------------------
% conf num
redundantConfPairsAll_sppH.num = [];
for i = 1:size(redundantConfSets_sppH,1)
    ind = find(redundantConfSets_sppH(i,:));
    redundantConfPairsAll_sppH.num = [redundantConfPairsAll_sppH.num; nchoosek(ind,2)];
end

% conf indices
% conf indices
redundantConfPairsAll_sppH.ind = zeros(size(redundantConfPairsAll_sppH.num));
for i = 1:size(redundantConfPairsAll_sppH.ind,1)
    redundantConfPairsAll_sppH.ind(i,:) = selectedConfGlobal_sppH.ind(redundantConfPairsAll_sppH.num(i,:));
end

% lets arrange them w.r.t distance values
dist_redundantConfPairsAll_sppH = zeros(size(redundantConfPairsAll_sppH.ind,1),1);
for i = 1:size(redundantConfPairsAll_sppH.ind,1)
    dist_redundantConfPairsAll_sppH(i) = distMatrix_selectedConfGlobal_sppH(redundantConfPairsAll_sppH.num(i,1),redundantConfPairsAll_sppH.num(i,2));
end
%dist_redundantConfPairsAll_sppH
[~,ind] = sort(dist_redundantConfPairsAll_sppH);
% now arrange them.
redundantConfPairsAll_sppH.num  = redundantConfPairsAll_sppH.num(ind,:);
redundantConfPairsAll_sppH.ind  = redundantConfPairsAll_sppH.ind(ind,:);
dist_redundantConfPairsAll_sppH = dist_redundantConfPairsAll_sppH(ind);

% pairs to test in the first place
% ------------------------------------------------------------------------------------
%{
% conf nums
% redundantConfPairs1stSPP_sppH.ind = sortrows(redundantConfPairs_sppH.ind)
pairs = [];
for i = 1:size(distMatrix_selectedConfGlobal_sppH,1)
    if any(distMatrix_selectedConfGlobal_sppH(i,:) <= DistForRedConfs_cell)
        if ~ismember(i,pairs)
            ind = find(distMatrix_selectedConfGlobal_sppH(i,:) <= DistForRedConfs_cell);
            ind(find(ismember(ind,pairs))) = [];
            if ind
                %nchoosek([i,ind],2)
                pairs = [pairs; nchoosek([i,ind],2)];
            end
        end
    end
end
redundantConfPairs1stSPP_sppH.num = pairs;

% conf indices
redundantConfPairs1stSPP_sppH.ind = zeros(size(redundantConfPairs1stSPP_sppH.num));
for i = 1:size(redundantConfPairs1stSPP_sppH.ind,1)
    redundantConfPairs1stSPP_sppH.ind(i,:) = selectedConfGlobal_sppH.ind(redundantConfPairs1stSPP_sppH.num(i,:));
end
%}


% --- plot (global) redundant conf cells
% figure('name','redundant confs'); imshow(map_env); hold on;
% plot(redundantConfGlobalCell_sppH_c,redundantConfGlobalCell_sppH_r,'*r');


% --- subject hotspots (vec) of redundant conf
% -----------------------------------------------------------
% subjectH_vec = zeros(size(hotspotSEQ));
% for i = 1:size(selectedConf_sppH_vec,2)
%     for j = 1:numel(redundantConfGlobal_sppH_ind)
%         
%         %find(find(confHotspots(:,i))==redundantConfGlobal_ind(j))
%         
%         % if any of selected conf for the hotspot is in the list of
%         % redundant confs.           
%         if find(find(selectedConf_sppH_vec(:,i))==redundantConfGlobal_sppH_ind(j))            
%             subjectH_vec(i) = 1;
%         end
%     end
% end
% subjectH_vec

%% SUBJECT HOTSPOTS

% % subject hotspots for all pairs of redundant conf
% % ------------------------------------------------------------------------------------
% subH_RedundantConfPairsAll_sppH.vec = zeros(size(redundantConfSets_sppH,1),size(hotspotSEQ,1));
% for i = 1:size(subH_RedundantConfPairsAll_sppH.vec,1)
%     conf_ind = find(redundantConfSets_sppH(i,:));
%     for j = 1:numel(conf_ind)
%          [~,hotspot_num] = find(selectedConf_sppH.ind==selectedConfGlobal_sppH.ind(conf_ind(j)));
%          subH_RedundantConfPairsAll_sppH.vec(i,hotspot_num) = 1;
%     end
% end


% subject hotspots for all pairs of redundant conf
% ------------------------------------------------------------------------------------
subH_RedundantConfPairsAll_sppH.vec = zeros(size(redundantConfPairsAll_sppH.num,1),size(hotspotSEQ,1));
for i = 1:size(subH_RedundantConfPairsAll_sppH.vec,1)
    %conf_ind = find(redundantConfSets_sppH(i,:));
    conf_num = redundantConfPairsAll_sppH.num(i,:);
    for j = 1:numel(conf_num)
         [~,hotspot_num] = find(selectedConf_sppH.ind==selectedConfGlobal_sppH.ind(conf_num(j)));
         subH_RedundantConfPairsAll_sppH.vec(i,hotspot_num) = 1;
    end
end

% check redundant conf pairs for cross angles.
% ------------------------------------------------------------------------------------

col_lines = lines(size(redundantConfPairsAll_sppH.num,1));

validity_vec = zeros(size(redundantConfPairsAll_sppH.num,1),1);
for i = 1:size(redundantConfPairsAll_sppH.num,1)
    
    conf_sub = selectedConfGlobalCell_sppH.sub(redundantConfPairsAll_sppH.num(i,:),:);
    
    redConf_num     = redundantConfPairsAll_sppH.num(i,:);
    %redConf_ind     = selectedConfGlobal_sppH.ind(redConf_num);
    redConfCell_sub = selectedConfGlobalCell_sppH.sub(redConf_num,:);
    
    subHotspots_num = find(subH_RedundantConfPairsAll_sppH.vec(i,:));
    %subHotspots_ind = Hotspots.ind(subHotspots_num);
    subHotspots_sub = Hotspots.sub(subHotspots_num,:);
    
    redConfTheta = zeros(numel(redConf_num),numel(subHotspots_num));    
    for j = 1:numel(redConf_num)        
        for k = 1:numel(subHotspots_num)            
            deltaR = redConfCell_sub(j,1)-subHotspots_sub(k,1);
            deltaC = redConfCell_sub(j,2)-subHotspots_sub(k,2);            
            redConfTheta(j,k) = atan2(deltaC,deltaR);            
        end
    end    
    %rad2deg(redConfTheta)
    %rad2deg(angleAbsDiff(redConfTheta(1,:),redConfTheta(2,:)))
    redConfCrossAnglesForSubHotspots = ...
        angleAbsDiff(redConfTheta(1,:),redConfTheta(2,:));
    
    if any(rad2deg(redConfCrossAnglesForSubHotspots)<=para_.XAnglesForRedundantConfs_deg)
        validity_vec(i) = 1;
    end
    
    if visualize == 1
    h_plot = plot(conf_sub(:,1),conf_sub(:,2),'s','Color',col_lines(i,:),'MarkerSize',5);
        
    disp('---- validity')
    %validity_vec(i)
    %pause
    
    delete(h_plot)
    end
    
end


% validity_vec

% size(redundantConfPairsAll_sppH.num)

valid_ind = find(validity_vec);
% update 
redundantConfPairsAll_sppH.num  = redundantConfPairsAll_sppH.num(valid_ind,:);
redundantConfPairsAll_sppH.ind  = redundantConfPairsAll_sppH.ind(valid_ind,:);
dist_redundantConfPairsAll_sppH = dist_redundantConfPairsAll_sppH(valid_ind);
subH_RedundantConfPairsAll_sppH.vec = subH_RedundantConfPairsAll_sppH.vec(valid_ind,:);

% size(redundantConfPairsAll_sppH.num)


% plot
% col_lines = colormap(lines(size(redundantConfPairsAll_sppH.num,1)));
% colormap(gray)
% for i = 1:size(redundantConfPairsAll_sppH.num,1)
%     
%     conf_sub = selectedConfGlobalCell_sppH.sub(redundantConfPairsAll_sppH.num(i,:),:);
%     plot(conf_sub(:,2),conf_sub(:,1),'s','Color',col_lines(i,:),'MarkerSize',5+(i*0.25));
%     
%     disp('---- validity')
%     validity_vec(i)
%     %hot_sub = Hotspots.sub(find(subjectH_redundantConfPairs1stSPP_sppH.vec(i,:)),:)
%     %plot(hot_sub(:,2),hot_sub(:,1),'+','Color',col_lines(i,:));
%     pause
% end


% plot for debugging
% {

if visualize == 1
col_lines = lines(size(redundantConfPairsAll_sppH.num,1));

for i = 1:size(redundantConfPairsAll_sppH.num,1)
    
    conf_sub = selectedConfGlobalCell_sppH.sub(redundantConfPairsAll_sppH.num(i,:),:);
    h_conf = plot(conf_sub(:,1),conf_sub(:,2),'s','Color',col_lines(i,:),'MarkerSize',8);
    
    
    subHotspots_num = find(subH_RedundantConfPairsAll_sppH.vec(i,:));
    %subHotspots_ind = Hotspots.ind(subHotspots_num);
    subHotspots_sub = Hotspots.sub(subHotspots_num,:);
    h_hot = plot(subHotspots_sub(:,1),subHotspots_sub(:,2),'og','MarkerSize',8);
    
    
    %disp('---- validity')
    %validity_vec(i)
    %hot_sub = Hotspots.sub(find(subjectH_redundantConfPairs1stSPP_sppH.vec(i,:)),:)
    %plot(hot_sub(:,2),hot_sub(:,1),'+','Color',col_lines(i,:));
    %pause
    delete(h_conf)
    delete(h_hot)
end
end
%}



end

% ----------------------------- End of the document --------------------------------------


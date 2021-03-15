function [ selectedConf_sppHF,...
           selectedConfGlobal_sppHF,...
           selectedConfGlobalCell_sppHF,...
           distMatrix_selectedConfGlobal_sppHF,...
           DistForRedConfs_cell,...
           minDistSelectedConfGlobal,...           
           redundantConfGlobal_sppHF,...
           redundantConfGlobalCell_sppHF,...
           redundantConfSets_sppHF,...
           redundantConfPairsAll_sppHF,...
           redundantConfPairs1stSPP_sppHF,...
           subH_RedundantConfPairs1stSPP_sppHF,...
           subH_RedundantConfPairsAll_sppHF] = fPreprocessFusionSetupHF( FoVfaces,...
                                                                      hotspotSEQ,...
                                                                      map_env,...
                                                                      SolutionLoc,...
                                                                      selectedConfGlobal_sppHF,...                                                                      
                                                                      selectedConf_sppHF,...
                                                                      paramPreprocess )
%fPreprocessFusionSetup identify the redundant pairs of selected conf from
%SPP-HOTSPOTS.



%% SENSING PARAMETERS AND CONSTANTS:
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
FoVfaces.num        = paramPreprocess.numFoVs;   % number of conf per cell.
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
priorH = load([SolutionLoc,'preprocess_coverage_crossangles_hotspot',num2str(1),'.mat'],...
    'oV','map_coverage');

% priorHF = load([SolutionLoc,'fusion_integration_1st_spp_pairs.mat']);


% initialize.
% selectedConfGlobal_sppHF.vec  = priorHF.selectedConfGlobal_sppHF.vec; %zeros(size(priorH.oV,2),1); % global conf vector
% selectedConf_sppHF.vec        = priorHF.selectedConf_sppHF.vec; %zeros(size(priorH.oV,2),numel(hotspotSEQ)); % conf vector for each hotspot


% figure('name','debug'); imshow(map_env); hold on;
       
% col_lines =[0         0    1.0000;...
%            0    0.5000         0;...
%            1.0000         0         0;...
%            0    0.7500    0.7500;...
%            0.7500         0    0.7500;...
%            0.7500    0.7500         0;...
%            0.2500    0.2500    0.2500;...
%            0         0    1.0000;...
%            0    0.5000         0;...
%            1.0000         0         0;...
%            0    0.7500    0.7500;...
%            0.7500         0    0.7500];

col_lines =[0         0    1.0000;...
           0    0.5000         0;...
           1.0000         0         0;...
           0    0.7500    0.7500;...
           0.7500         0    0.7500;...
           0.7500    0.7500         0;...
           0.2500    0.2500    0.2500;...
           0         0    1.0000;...
           0    0.5000         0;...
           1.0000         0         0;...
           0    0.7500    0.7500;...
           0         0    1.0000;...
           0    0.5000         0;...
           1.0000         0         0;...
           0    0.7500    0.7500;...

           0.7500         0    0.7500];    
% figure('name','redundant'); imshow(map_env); hold on;
figure('name','redundant'); imshow(map_env,'InitialMagnification',700); hold on;


% for each hotspot
%{
for i = 1:numel(hotspotSEQ)
    
    % load preprocessing
    priorH = load([SolutionLoc,'preprocess_coverage_crossangles_hotspot',num2str(i),'.mat'],...
        'oV','V','confKept','conf_theta','conf_crossAngles_G','refDistPoint','D');
    
    % load SPP solution for the center of mass
    selectedConf = load([SolutionLoc,'/spp_selectedconf_hotspot',num2str(i),'.dat']);    
    selectedConf = round(selectedConf);
    
    % conf vector with all cells in the map.
    conf = zeros(size(priorH.oV,2),1);
    conf(priorH.confKept) = selectedConf;
    
    % update overall conf vector.
    selectedConfGlobal_sppHF.vec = selectedConfGlobal_sppHF.vec+conf;
    
    % update conf vector for all hotspots
    selectedConf_sppHF.vec(:,i) = conf;
    
    % orientations of the selected conf in each spp-hotspots
    selectedConfTheta_sppHF.vec(priorH.confKept(find(selectedConf)),i) =...
        priorH.conf_theta(find(selectedConf));
    
    
    conf_debug_ind = priorH.confKept(find(selectedConf));
    confcell_debug_ind = fix(conf_debug_ind/f)+(~~mod(conf_debug_ind,f)); 
        
    [confcell_debug_r,confcell_debug_c] = ind2sub(size(map_env),confcell_debug_ind);
    plot(confcell_debug_c,confcell_debug_r,'x','Color',col_lines(i,:));
    
    hot_debug_ind = hotspotSEQ(i);
    [hot_debug_r,hot_debug_c] = ind2sub(size(map_env),hot_debug_ind);
    plot(hot_debug_c,hot_debug_r,'o','Color',col_lines(i,:));
    
    
    
    
    
    % spp-hotspots solution/objective values
    C = selectedConf;
    G = alpha*priorH.conf_crossAngles_G;
    
    Uc = gamma*((ones(size(priorH.V,1),1)'*priorH.V)/max((ones(size(priorH.V,1),1)'*priorH.V)));
    Ud = (1-gamma)*(1-(priorH.D/max(priorH.D)));
    U  = beta*(Uc'+Ud);
    
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
%}

% {
for i = 1:numel(hotspotSEQ)
    
    conf_debug_ind = find(selectedConf_sppHF.vec(:,i));
    confcell_debug_ind = fix(conf_debug_ind/f)+(~~mod(conf_debug_ind,f)); 
    
    [confcell_debug_r,confcell_debug_c] = ind2sub(size(map_env),confcell_debug_ind);
    plot(confcell_debug_c,confcell_debug_r,'x','Color',col_lines(i,:));
    
    hot_debug_ind = hotspotSEQ(i);
    [hot_debug_r,hot_debug_c] = ind2sub(size(map_env),hot_debug_ind);
    plot(hot_debug_c,hot_debug_r,'o','Color',col_lines(i,:));
    %pause
end
%}


% selected conf indices
% -----------------------------------------------------------
% selectedConf_sppHF.ind = zeros(NumberOfSelectedConf_sppH,size(selectedConf_sppHF.vec,2));
% for i = 1:size(selectedConf_sppHF.vec,2)
%     selectedConf_sppHF.ind(:,i) = find(selectedConf_sppHF.vec(:,i));
% end

% selected golbal conf indices
% -----------------------------------------------------------
% selectedConfGlobal_sppHF.ind = find(selectedConfGlobal_sppHF.vec); 
% selectedConfGlobal_sppHF

% selected global conf cells (indices/subscripts)
% -----------------------------------------------------------
% indices of global conf cells
selectedConfGlobalCell_sppHF.ind = fix(selectedConfGlobal_sppHF.ind/f)+(~~mod(selectedConfGlobal_sppHF.ind,f)); 
 %subscripts of global conf indices
[r,c] = ind2sub(size(map_env),selectedConfGlobalCell_sppHF.ind);
selectedConfGlobalCell_sppHF.sub = [r,c];

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
distMatrix_selectedConfGlobal_sppHF = zeros(numel(selectedConfGlobalCell_sppHF.ind));
for i = 1:numel(selectedConfGlobalCell_sppHF.ind)
    for j = 1:numel(selectedConfGlobalCell_sppHF.ind)
        distMatrix_selectedConfGlobal_sppHF(i,j) =...
            sqrt(((selectedConfGlobalCell_sppHF.sub(j,2)-selectedConfGlobalCell_sppHF.sub(i,2))^2)+...
                 ((selectedConfGlobalCell_sppHF.sub(j,1)-selectedConfGlobalCell_sppHF.sub(i,1))^2));
    end
end
% distance to itself is inf
for i = 1:numel(selectedConfGlobalCell_sppHF.ind)
    distMatrix_selectedConfGlobal_sppHF(i,i) = inf;
end


% check for redundant conf
% -----------------------------------------------------------
% min distance (in cells) to check redundant conf.
DistForRedConfs_cell = paramPreprocess.DistanceForRedundantConfs_m/paramPreprocess.MapRes;
% --- finding redundant conf from the distance matrix
minDistSelectedConfGlobal.vec = zeros(size(distMatrix_selectedConfGlobal_sppHF,1),1);

for i = 1:size(distMatrix_selectedConfGlobal_sppHF,1)
    %[i,find(distMatrix_selectedConfGlobal_sppH(i,:) <= DistForRedConfs_cell)]
    if any(distMatrix_selectedConfGlobal_sppHF(i,:) <= DistForRedConfs_cell)
        minDistSelectedConfGlobal.vec(i) = 1;
    end
end
% minDistSelectedConfGlobal.vec


% indices/vector of global redundant conf.
% ------------------------------------------------------------------------------------
% indices
redundantConfGlobal_sppHF.ind = selectedConfGlobal_sppHF.ind(find(minDistSelectedConfGlobal.vec));
% vec
redundantConfGlobal_sppHF.vec = zeros(size(priorH.oV,2),1);
redundantConfGlobal_sppHF.vec(redundantConfGlobal_sppHF.ind) = 1;


% indices/subscripts of global redundant conf cells.
% ------------------------------------------------------------------------------------
redundantConfGlobalCell_sppHF.ind = selectedConfGlobalCell_sppHF.ind(find(minDistSelectedConfGlobal.vec));
[r,c] = ind2sub(size(map_env),redundantConfGlobalCell_sppHF.ind);
redundantConfGlobalCell_sppHF.sub = [r,c];

% sets of redundant confs.
% -----------------------------------------------------------
redundantConfSets_sppHF = zeros(numel(selectedConfGlobalCell_sppHF.ind));
for i = 1:numel(selectedConfGlobalCell_sppHF.ind)
    redundantConfSets_sppHF(i,i) = 1; % include ith conf vs ith conf
    redundantConfSets_sppHF(i,find(distMatrix_selectedConfGlobal_sppHF(i,:)<= DistForRedConfs_cell)) = 1;    
end
% remove duplicate rows
redundantConfSets_sppHF = unique(redundantConfSets_sppHF,'rows');

% QUESTIONABLE
while sum(sum(redundantConfSets_sppHF,1),2)~=size(redundantConfSets_sppHF,2)
    B = [];
    for i = 1:size(redundantConfSets_sppHF,2)
        B = [B;sum(redundantConfSets_sppHF(find(redundantConfSets_sppHF(:,i)),:),1)];        
    end
    C = unique(B,'rows');
    D = zeros(size(C));
    D(find(C)) = 1;
    redundantConfSets_sppHF = D;
end

% finding single conf rows (to be removed from the list)
toRemove_ind = zeros(size(redundantConfSets_sppHF,1),1);
for i = 1:size(redundantConfSets_sppHF,1)
    if numel(find(redundantConfSets_sppHF(i,:)))<=1
        toRemove_ind(i) = 1;
    end
end
% now remove them
redundantConfSets_sppHF(find(toRemove_ind),:) = [];
% redundantConfSets_sppH


% all pairs of redundant conf
% ------------------------------------------------------------------------------------
% conf num
redundantConfPairsAll_sppHF.num = [];
for i = 1:size(redundantConfSets_sppHF,1)
    ind = find(redundantConfSets_sppHF(i,:));
    redundantConfPairsAll_sppHF.num = [redundantConfPairsAll_sppHF.num; nchoosek(ind,2)];
end

% conf indices
% conf indices
redundantConfPairsAll_sppHF.ind = zeros(size(redundantConfPairsAll_sppHF.num));
for i = 1:size(redundantConfPairsAll_sppHF.ind,1)
    redundantConfPairsAll_sppHF.ind(i,:) = selectedConfGlobal_sppHF.ind(redundantConfPairsAll_sppHF.num(i,:));
end


% pairs to test in the first place
% ------------------------------------------------------------------------------------
% conf nums
% redundantConfPairs1stSPP_sppH.ind = sortrows(redundantConfPairs_sppH.ind)
pairs = [];
for i = 1:size(distMatrix_selectedConfGlobal_sppHF,1)
    if any(distMatrix_selectedConfGlobal_sppHF(i,:) <= DistForRedConfs_cell)
        if ~ismember(i,pairs)
            ind = find(distMatrix_selectedConfGlobal_sppHF(i,:) <= DistForRedConfs_cell);
            ind(find(ismember(ind,pairs))) = [];
            if ind
                %nchoosek([i,ind],2)
                pairs = [pairs; nchoosek([i,ind],2)];
            end
        end
    end
end
redundantConfPairs1stSPP_sppHF.num = pairs;

% conf indices
redundantConfPairs1stSPP_sppHF.ind = zeros(size(redundantConfPairs1stSPP_sppHF.num));
for i = 1:size(redundantConfPairs1stSPP_sppHF.ind,1)
    redundantConfPairs1stSPP_sppHF.ind(i,:) = selectedConfGlobal_sppHF.ind(redundantConfPairs1stSPP_sppHF.num(i,:));
end

% to add prior variables
priorF = load([SolutionLoc,'preprocess_fusion_setup_H.mat'],...
                            'redundantConfPairs1stSPP_sppH',...
                            'subH_RedundantConfPairs1stSPP_sppH',...
                            'subH_RedundantConfPairsAll_sppH');

redundantConfPairs1stSPP_sppHF.ind = [priorF.redundantConfPairs1stSPP_sppH.ind;...
                                             redundantConfPairs1stSPP_sppHF.ind];

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


% subject hotspots for all pairs of redundant conf
% ------------------------------------------------------------------------------------
subH_RedundantConfPairsAll_sppHF.vec = zeros(size(redundantConfSets_sppHF,1),size(hotspotSEQ,1));
for i = 1:size(subH_RedundantConfPairsAll_sppHF.vec,1)
    conf_num = find(redundantConfSets_sppHF(i,:));    
    conf_ind = selectedConfGlobal_sppHF.ind(conf_num);    
    for j = 1:numel(conf_ind)        
        hotspot_num = find(selectedConf_sppHF.vec(conf_ind(j),:));
        %[~,hotspot_num] = find(selectedConf_sppHF.ind==selectedConfGlobal_sppHF.ind(conf_ind(j)));
        subH_RedundantConfPairsAll_sppHF.vec(i,hotspot_num) = 1;
    end
end

% add previous list
subH_RedundantConfPairsAll_sppHF.vec = [priorF.subH_RedundantConfPairsAll_sppH.vec;...
                                        subH_RedundantConfPairsAll_sppHF.vec];

% subject hotspots for 1st pairs of redundant conf
% ------------------------------------------------------------------------------------
subH_RedundantConfPairs1stSPP_sppHF.vec = zeros(size(redundantConfPairs1stSPP_sppHF.num,1),size(hotspotSEQ,1));
for i = 1:size(subH_RedundantConfPairs1stSPP_sppHF.vec,1)
    %conf_ind = find(redundantConfSets_sppH(i,:));
    for j = 1:2%numel(conf_ind)
        %[~,hotspot_num] = find(selectedConf_sppHF.ind==selectedConfGlobal_sppHF.ind(redundantConfPairs1stSPP_sppH.num(i,j)));
        hotspot_num = find(selectedConf_sppHF.vec(redundantConfPairs1stSPP_sppHF.ind(j),:));
        subH_RedundantConfPairs1stSPP_sppHF.vec(i,hotspot_num) = 1;
    end
end
% add previous list
subH_RedundantConfPairs1stSPP_sppHF.vec = [priorF.subH_RedundantConfPairs1stSPP_sppH.vec;...
                                           subH_RedundantConfPairs1stSPP_sppHF.vec];


% plot
for i = 1:size(redundantConfPairs1stSPP_sppHF.num,1)
    
    conf_sub = selectedConfGlobalCell_sppHF.sub(redundantConfPairs1stSPP_sppHF.num(i,:),:);
    plot(conf_sub(:,2),conf_sub(:,1),'s','Color',col_lines(i,:));
    
    %hot_sub = Hotspots.sub(find(subjectH_redundantConfPairs1stSPP_sppH.vec(i,:)),:)
    %plot(hot_sub(:,2),hot_sub(:,1),'+','Color',col_lines(i,:));
    
end




%% CONDENSED CONFIGURATIONS
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%{
% condensed conf (bin) vector for each hotspot
condensedConf_sppH_vec = selectedConf_sppH_vec;
for i = 1:size(condensedConf_sppH_vec,2)    
    condensedConf_sppH_vec(find(selectedConf_sppH_vec(:,i).*redundantConfGlobal_sppH_vec),i) = 0;
    
    find(condensedConf_sppH_vec(:,i))
    
end

% condensed conf counts for each hotspot
condensedConf_sppH_num = zeros(1,size(condensedConf_sppH_vec,2));
for i = 1:numel(condensedConf_sppH_num)    
    find(condensedConf_sppH_vec(:,i))    
    condensedConf_sppH_num(i) = numel(find(condensedConf_sppH_vec(:,i)));
end
condensedConf_sppH_num



% ---

selectedConf_sppH.vec = selectedConf_sppH_vec
selectedConf_sppH.ind = selectedConf_sppH_ind

selectedConfGlobal_sppH.ind = selectedConfGlobal_sppH_ind
selectedConfGlobal_sppH.vec = selectedConfGlobal_sppH_vec

selectedConfGlobalCell_sppH.ind = selectedConfGlobalCell_sppH_ind
selectedConfGlobalCell_sppH.sub = [selectedConfGlobalCell_sppH_r,selectedConfGlobalCell_sppH_c]

selectedConfTheta_sppH.vec = selectedConfTheta_sppH_vec

redundantConfGlobal_sppH.ind = redundantConfGlobal_sppH_ind
redundantConfGlobal_sppH.vec = redundantConfGlobal_sppH_vec
redundantConfGlobalCell_sppH.ind = redundantConfGlobalCell_sppH_ind
redundantConfGlobalCell_sppH.sub =[redundantConfGlobalCell_sppH_r,redundantConfGlobalCell_sppH_c]

condensedConf_sppH.vec = condensedConf_sppH_vec
condensedConf_sppH.num = condensedConf_sppH_num

infoGain_sppH
infoGainG_sppH
maxDs_sppH
max1Vs_sppH
referenceDistPoint_sppH

distMatrix_selectedConfGlobal_sppH
DistForRedConfs_cell
minDistSelectedConfGlobal_vec

redundantConfSets_sppH
subjectH_RedundantConfPairs_sppH_vec



% plot
% figure('name','pair wise redundant confs'); imshow(map_env); hold on;

for i = 1:size(redundantConfSets_sppH,1)
    conf_ind_sub = selectedConfGlobalCell_sppH.sub(find(redundantConfSets_sppH(i,:)),:);
    H_pairwise_conf = plot(conf_ind_sub(:,2),conf_ind_sub(:,1),'s','Color',col_lines(i,:));
    %pause(5);
    %delete(H_pairwise_conf)
end


pause

% redundant conf pairs with hotspot numbers.
redundantConfPairs_hotNum_sppH = zeros(size(redundantConfSets_sppH));
for i = 1:size(redundantConfSets_sppH,1)
    conf_num = find(redundantConfSets_sppH(i,:));    
    conf_ind = selectedConfGlobal_sppH_ind(conf_num);    
    for j = 1:numel(conf_ind)        
        [~,hot_num] = find(selectedConf_sppH_ind==conf_ind(j));        
        redundantConfPairs_hotNum_sppH(i,conf_num(j)) = hot_num;
    end
end


redundantConfPairs_hotNum_sppH

redundantConfSets_sppH



selectedConfGlobal_sppH_ind
selectedConf_sppH_ind

%}


end

% ---------------------------------------------------- End of the document -------------------------------------------------------


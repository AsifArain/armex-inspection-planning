function [ ConfPair_ioh ] = fIterativeFusionSetup1t( map_env,...
                                                     PairNumberToSelect,...
                                                     confs_executed_hc,...
                                                     currentFusedConfs,...
                                                     cellsize_env,...
                                                     para_ )
%fFusionSetup identify the redundant pairs of selected conf from the local
%solutions.



%**********************************************************************
% -- list of conf with index and orientation
%**********************************************************************
% conf_io = [selectedFusedConf_ind(:),selectedFusedConf_orn(:)];
% conf_io = sortrows(conf_io,1);

conf_io = [(cell2mat(currentFusedConfs.ind))',(cell2mat(currentFusedConfs.orn))'];
% conf_io = sortrows(conf_io,1);
conf_io = unique(conf_io,'rows');

%-- remove executed conf
for i = 1:size(confs_executed_hc,1)
    conf_io(ismember(conf_io(:,1),confs_executed_hc(i,1)),:) = [];
end
  
%-- corresponding (local) hc number
conf_ioh = [conf_io,zeros(size(conf_io,1),numel(currentFusedConfs.ind))];
for i = 1:size(conf_ioh,1)
    for j = 1:numel(currentFusedConfs.ind)        
        if any(ismember(currentFusedConfs.ind{j},conf_ioh(i,1)) &...
               ismember(currentFusedConfs.orn{j},conf_ioh(i,2)))
           %conf_ioh(i,3) = j;
           conf_ioh(i,2+j) = 1;
        end
    end
end


%**********************************************************************
% -- distance matrix between all the confs.
%**********************************************************************
[conf_r,conf_c] = ind2sub(size(map_env),conf_ioh(:,1));
dist = zeros(size(conf_ioh,1));
for i = 1:size(dist,1)
    % -- dist from this conf to all
    dist(i,:) = pdist2([conf_r(i),conf_c(i)],[conf_r,conf_c],'euclidean');
    
    % -- dist from this conf to this conf
    %dist(i,i) = inf;
    
    % -- dist from this conf to the other confs in the same hotspots    
    % - hotspots number
    %{
    [hot_num,~] = find(ismember(selectedFusedConf_ind,conf_io(i,1)) &...
                       ismember(selectedFusedConf_orn,conf_io(i,2)));
                   
    [hot_num,~] = find(ismember(localConf.ind{:},conf_io(i,1)) &...
                       ismember(localConf.orn{:},conf_io(i,2)));
                   
    [hot_num,~] = find(ismember(localConf.ind{:},5) &...
                       ismember(localConf.orn{:},15));
                   
    [max_size, max_index] = max(cellfun('size', localConf.ind, 5))
    %}
    
    % - conf numbers in the hotspots
    %ind_shared_hot = ismember(conf_io(:,1),selectedFusedConf_ind(hot_num,:));
        
    this_conf_hc_num = find(conf_ioh(i,3:end));
    for j = 1:numel(this_conf_hc_num)
        ind_shared_hot = find(conf_ioh(:,2+this_conf_hc_num(j)));
        dist(i,ind_shared_hot) = inf;
    end
       
    %{
    ind_shared_hot = ismember(conf_ioh(:,3),conf_ioh(i,3));
    
    % - dist to the confs is inf
    dist(i,ind_shared_hot) = inf;
    %}
end


%**************
confABcomb_sub = nchoosek(1:size(conf_io,1),2);
confABcomb_ind = sub2ind(size(dist),confABcomb_sub(:,1),confABcomb_sub(:,2));
confAB_dist = sortrows([confABcomb_sub,dist(confABcomb_ind)],3);

% -- distance to look for the pairs
% DistForRedConfs_cell = para_.DistanceForRedundantConfs_m/dataYAML.resolution;
DistForRedConfs_cell = para_.DistanceForRedundantConfs_m/cellsize_env;

confAB_dist(confAB_dist(:,3)>DistForRedConfs_cell,:) = [];


if PairNumberToSelect>size(confAB_dist,1)
    ConfPair_ioh = [];
else
    
    conf_num_ab = confAB_dist(PairNumberToSelect,1:2);
    ConfPair_ioh = conf_ioh(conf_num_ab,:);
    
end

%{
%**********************************************************************
% -- conf pair with minimum distance
%**********************************************************************
[val,ind] = sort(dist(:));


% -- distance to look for the pairs
DistForRedConfs_cell = para_.DistanceForRedundantConfs_m/dataYAML.resolution;

if all(val(PairNumberToSelect:end)>DistForRedConfs_cell)
    ConfPair_ioh = [];
    return
end


[conf_num_a,conf_num_b] = ind2sub(size(dist),ind(PairNumberToSelect));
conf_num_ab = sort([conf_num_a;conf_num_b]);


%**********************************************************************
% -- minimum dist conf pair with index, orinetation and hotspots
%**********************************************************************
%{
% - initialize the pair
ConfPair_ioh = zeros(2,2+numel(localConf.ind));

% - index and orinetations
ConfPair_ioh(:,1:2) = conf_ioh(conf_num_ab,1:2);

% - respective hotspots
%{
ind_a = find(selectedFusedConf_ind(:)==ConfPair_ioh(1,1) &...
    selectedFusedConf_orn(:)==ConfPair_ioh(1,2));

ind_b = find(selectedFusedConf_ind(:)==ConfPair_ioh(2,1) &...
    selectedFusedConf_orn(:)==ConfPair_ioh(2,2));

[hots_a,~] = ind2sub(size(selectedFusedConf_ind),ind_a);
[hots_b,~] = ind2sub(size(selectedFusedConf_ind),ind_b);


ConfPair_ioh(1,2+hots_a) = 1;
ConfPair_ioh(2,2+hots_b) = 1;
%}

ConfPair_ioh(1,2+hots_a) = 1;
ConfPair_ioh(2,2+hots_b) = 1;

%}

ConfPair_ioh = conf_ioh(conf_num_ab,:);

%}

end

%------------------------------- End of the document -------------------------------------


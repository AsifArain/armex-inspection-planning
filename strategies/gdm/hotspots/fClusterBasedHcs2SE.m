function fClusterBasedHcs2SE
%


map_hotspot = zeros(size(map_recon));
% ind = find(gd_map_hotspot(:))>500;
% gd_map_hotspot(ind) = 1;
map_hotspot((map_recon(:)>para_.highConcentrationThreshold_PPM)) = 1;
% test_map_resized2_obs = find(test_map_resized2(:)==0);
% gd_map_coverage(test_map_resized2_obs) = 0;


% -- hotspots based on clustering

% - high concentration cells
[high_r,high_c] = find(map_hotspot);
Hcells = [high_r,high_c];


% plot
% visualize = 1
if visualize==1
    figure('name','Positive concentrations'); 
    imshow(map_env','InitialMagnification',350); hold on;
    set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
    color_positive_con = [255,140,0]/255;
    pause(1)
    for i = 1:size(Hcells,1)
        plot(Hcells(i,1),Hcells(i,2),...
             's',...
             'color',color_positive_con,...
             'MarkerFaceColor',color_positive_con,...
             'MarkerSize',3.75); %2.75
    end
    pause(2)
    filename = 'hotspots-positive-concentrations';
    % print('-painters','-dpdf','-r500',filename); % pdf    
    export_fig([dir_.TomographyMaps,filename], '-pdf','-r500')
    pause(1)
end


%*******************************************************************
% 
%               CLUSTERING
% 
%*******************************************************************

%-- linkage cutoff for clusters
cut_off = 3.5; %median([Z(end-2,3) Z(end-1,3)])


%--- filter outliers
%---------------------------------------------
if size(Hcells,1)>1
    
    
    %-- cluster linkage tree
    Z = linkage(Hcells,'single');
    %figure; dendrogram(Z,'ColorThreshold',cut_off)

    % - high concentration clusters
    clusters2Filter = cluster(Z,'cutoff',cut_off,'criterion','distance');
        
    % - number of clusters
    num_of_clusters = numel(unique(clusters2Filter));
    
    %-- weighted map
    Wmap = map_recon./sum(map_recon(:)); % total weight is 1    
    
    %-- filter
    filter_ind = [];
    for i = 1:num_of_clusters
        
        %-- this cluster
        thisCluster = Hcells(clusters2Filter==i,:);
        
        %-- if number of observation are less than a number
        if size(thisCluster,1) < 5
            
            %-- concentration weight of this cluster
            thisWeight = 0;
            for j = 1:size(thisCluster,1)
                thisWeight = thisWeight + Wmap(thisCluster(j,1),thisCluster(j,2));
            end
            %thisWeight
            %-- if concentration weight of this cluster is less than a
            %number
            if thisWeight < 0.01
                %-- filter out the observations (cells)
                %Hcells_filtered(clusters2Filter==i,:) = [];
                
                filter_ind = [filter_ind;find(clusters2Filter==i)];
                
            end
        end
    end
    
    Hcells(filter_ind,:) = [];
    
    
    
    % plot
    if visualize==1
        figure('name','Filtered observations (cells)'); 
        imshow(map_env','InitialMagnification',350); hold on;
        set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
        pause(1)
        for i = 1:size(Hcells,1)
            plot(Hcells(i,1),Hcells(i,2),...
                 's',...
                 'color','k',...%color_positive_con,...
                 'MarkerFaceColor',color_positive_con,...
                 'MarkerSize',3.75); %2.75
        end
        pause(2)
        filename = 'hotspots-filtered-cells';
        % print('-painters','-dpdf','-r500',filename); % pdf    
        export_fig([dir_.TomographyMaps,filename], '-pdf','-r500')
        pause(1)
    end
        
end




if size(Hcells,1)>1
    
    
    %*******************************************
    %               FIRST LAYER
    %*******************************************

    % - cluster linkage tree
    Z = linkage(Hcells,'single');
    figure; dendrogram(Z);
    
    %{
    Z = linkage(Hcells,'average');
    figure; dendrogram(Z);
    
    Z = linkage(Hcells,'centroid');
    figure; dendrogram(Z);
    
    Z = linkage(Hcells,'complete');
    figure; dendrogram(Z);
    
    Z = linkage(Hcells,'median');
    figure; dendrogram(Z);
    
    Z = linkage(Hcells,'ward');
    figure; dendrogram(Z);
    
    Z = linkage(Hcells,'weighted');
    figure; dendrogram(Z);
    %}
        
    %figure; dendrogram(Z,'ColorThreshold',cut_off); hold on;
    %plot([0,50],[cut_off,cut_off],'-r');
    
    % - high concentration clusters
    Hclusters1 = cluster(Z,'cutoff',cut_off,'criterion','distance');
    
    
    % Hclusters = cluster(Z,'cutoff',1);
    % Hclusters = cluster(Z,'cutoff',2.5,'Depth',2);
    % Hclusters = cluster(Z,'cutoff',1,'MaxClust',5);

    % - number of clusters
    num_of_clusters1 = numel(unique(Hclusters1));


    % plot
    % visualize = 1
    if visualize==1
    figure('name','clusters - 1st layer'); 
    imshow(map_env','InitialMagnification',350); hold on;
    set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
    pause(1)
    % colors = lines(num_of_clusters1);
    %colors = lines(9);
    % colors = linspecer(20);
    %colors = linspecer(30);
    colors = linspecer(num_of_clusters1);
    for i = 1:num_of_clusters1
        plot(Hcells(Hclusters1==i,1),Hcells(Hclusters1==i,2),...
             's',...
             'color','k',...%colors(i,:),...
             'MarkerFaceColor',colors(i,:),...
             'MarkerSize',3.75);
    end    
    pause(2)
    filename = 'hotspots-clusters1-agglomerative';
    % print('-painters','-dpdf','-r500',filename); % pdf    
    export_fig([dir_.TomographyMaps,filename], '-pdf','-r500')
    pause(1)    
    end

    %*******************************************
    %           SECOND LAYER
    %*******************************************


    hot_clusters = zeros(size(Hcells,1),1);
    clust_num = 1;

    for i = 1:num_of_clusters1

        cluster1_this = [Hcells(Hclusters1==i,1),Hcells(Hclusters1==i,2)];

        dist_this = zeros(size(cluster1_this,1));
        for j = 1:size(dist_this,1)
            dist_this(j,:) = pdist2(cluster1_this(j,:),cluster1_this,'euclidean');
        end

        num_of_desired_clusters2_this = ...
            round(max(dist_this(:))/(SensorRange.cell/2));
        if num_of_desired_clusters2_this<1
            num_of_desired_clusters2_this = 1;
        end


        if num_of_desired_clusters2_this>1

            %if num_of_clusters2_this<1
            %    num_of_clusters2_this = 1;
            %end

            %cluster_tree2_this = linkage(cluster1_this,'single');    
            %Hclusters2_this = cluster(cluster_tree2_this,'MaxClust',num_of_clusters2_this)
            [Hclusters2_this,C] = kmeans(cluster1_this,num_of_desired_clusters2_this,'Distance','cityblock');
            %[Hclusters2_this,C] = kmeans(cluster1_this,num_of_clusters2_this,'Distance','cosine');
            %[Hclusters2_this,C] = kmeans(cluster1_this,num_of_clusters2_this,'Distance','correlation');
            %[Hclusters2_this,C] = kmeans(cluster1_this,num_of_clusters2_this,'Distance','hamming'); % beykar


            % - number of clusters
            num_of_desired_clusters2_this = numel(unique(Hclusters2_this));
            hot_clusters(Hclusters1==i) = Hclusters2_this+clust_num-1;
            clust_num = clust_num+num_of_desired_clusters2_this;

        else

            hot_clusters(Hclusters1==i) = clust_num;
            clust_num = clust_num+1;
        end

        %pause

    end


    % - number of clusters
    num_of_clusters = numel(unique(hot_clusters));



    % plot
    if visualize==1
    figure('name','clusters - 2nd layer'); 
    imshow(map_env','InitialMagnification',350); hold on;
    set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
    pause(1)
    % colors = lines(num_of_clusters);
    colors = linspecer(num_of_clusters1+num_of_clusters);
    for i = 1:num_of_clusters
        plot(Hcells(hot_clusters==i,1),Hcells(hot_clusters==i,2),...
             's',...
             'color',colors(num_of_clusters1+i,:),...
             'MarkerFaceColor',colors(num_of_clusters1+i,:),...
             'MarkerSize',3.75);
    end
    pause(2)
    filename = 'hotspots-clusters2-kmeans';
    % print('-painters','-dpdf','-r500',filename); % pdf    
    export_fig([dir_.TomographyMaps,filename], '-pdf','-r500')
    pause(1)
    end 

    % -- hotspots are weighted mean
    Wmap = map_recon./sum(map_recon(:)); % total weight is 1
    hotspots = zeros(num_of_clusters,2);
    hcs_weights = zeros(num_of_clusters,1);
    for i = 1:num_of_clusters

        Hcluster = [Hcells(hot_clusters==i,1),Hcells(hot_clusters==i,2)];
        
        Hcluster_wmean_row = 0;
        Hcluster_wmean_col = 0;
        Hcluster_wmean_wgt = 0;
        for j = 1:size(Hcluster,1)
            %
            Hcluster_wmean_row = Hcluster_wmean_row + ...
                Wmap(Hcluster(j,1),Hcluster(j,2))*Hcluster(j,1);
            
            Hcluster_wmean_col = Hcluster_wmean_col + ...
                Wmap(Hcluster(j,1),Hcluster(j,2))*Hcluster(j,2);
            
            Hcluster_wmean_wgt = Hcluster_wmean_wgt + ...
                Wmap(Hcluster(j,1),Hcluster(j,2));
        end

        Hcluster_wmean = [Hcluster_wmean_row/Hcluster_wmean_wgt,...
                          Hcluster_wmean_col/Hcluster_wmean_wgt];

        hotspots(i,:) = round(Hcluster_wmean);
        hcs_weights(i) = Hcluster_wmean_wgt;
        
    end
    
    hcw_before = [hotspots,hcs_weights];
    
    %-- moving hcs from occupied cells to the nearest unoccupied cells
    [ hotspots ] = fMovingOccupiedHcsToFreeSpace( hotspots,...
                                                  OBS,...
                                                  map_env,...
                                                  map_hotspot,...
                                                  dataYAML,...
                                                  para_ );
    hcw_after = [hotspots,hcs_weights];
    
    find(hcw_after(:)~=hcw_before(:))
    %pause
    
    if visualize==1
        
        figure('name','clusters-2 with hotspots'); 
        imshow(map_env','InitialMagnification',350); hold on;
        set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
        pause(1)
        
        colors = lines(num_of_clusters);
        for i = 1:num_of_clusters
            plot(Hcells(hot_clusters==i,1),Hcells(hot_clusters==i,2),...
                 's',...
                 'color',colors(i,:),...
                 'MarkerFaceColor',colors(i,:),...
                 'MarkerSize',4); %5
        end        
        % -- hotspots
        for i = 1:size(hotspots,1)
            plot(hotspots(i,1),hotspots(i,2),...
                 'o','color','k','MarkerFaceColor','r','MarkerSize',5); %7
        end
        pause(2)
        filename = 'hotspots-clusters2-kmeans-hc';
        % print('-painters','-dpdf','-r500',filename); % pdf    
        export_fig([dir_.TomographyMaps,filename], '-pdf','-r500')
        pause(1)
    end
    
    
    
    %*******************************************
    %           nth LAYER
    %*******************************************
    
    hcs_prevCluster = hotspots;
    ws_prevCluster  = hcs_weights;
    
    cond_nClustering = size(hcs_prevCluster,1)>1;
    
    
    while cond_nClustering
    
        % - cluster linkage tree
        Z = linkage(hcs_prevCluster,'single');
        %figure; dendrogram(Z);

        % - high concentration clusters
        clusters_this = cluster(Z,'cutoff',8.5,'criterion','distance');

        % - number of clusters
        num_of_clusters_this = numel(unique(clusters_this));
        
        % -- hotspots are weighted mean
        hcs_thisCluster = zeros(num_of_clusters_this,2);
        ws_thisCluster  = zeros(num_of_clusters_this,1);
        
        for i = 1:num_of_clusters_this

            candidateCells = [hcs_prevCluster(clusters_this==i,1),...
                              hcs_prevCluster(clusters_this==i,2)];
            candiateWeights = ws_prevCluster(clusters_this==i);
            
            thisHc_r = 0;
            thisHc_c = 0;
            thisHc_w = 0;
            for j = 1:size(candidateCells,1)
                %
                thisHc_r = thisHc_r + candiateWeights(j)*candidateCells(j,1);
                thisHc_c = thisHc_c + candiateWeights(j)*candidateCells(j,2);
                thisHc_w = thisHc_w + candiateWeights(j);
            end
            
            thisHc = [thisHc_r/thisHc_w,thisHc_c/thisHc_w];

            hcs_thisCluster(i,:) = round(thisHc);
            ws_thisCluster(i) = thisHc_w;

        end
        
        if visualize==1
        
            figure('name','clusters-n with hotspots'); 
            imshow(map_env','InitialMagnification',350); hold on;
            set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
            pause(1)

            colors = lines(num_of_clusters_this);
            for i = 1:num_of_clusters_this
                plot(hcs_prevCluster(clusters_this==i,1),...
                     hcs_prevCluster(clusters_this==i,2),...
                     's',...
                     'color','k',...
                     'MarkerFaceColor',colors(i,:),...
                     'MarkerSize',4); %5
            end

            % -- hotspots
            for i = 1:size(hcs_thisCluster,1)
                plot(hcs_thisCluster(i,1),hcs_thisCluster(i,2),...
                     'o','color','k','MarkerFaceColor','r','MarkerSize',5); %7
            end
        end
        
        cond_nClustering = size(hcs_thisCluster,1)>1 & ...
                           (size(hcs_prevCluster,1)-size(hcs_thisCluster,1))>0;
        
        hcs_prevCluster = hcs_thisCluster;
        ws_prevCluster = ws_thisCluster;
        
        %pause
        
    end
    
    if visualize==1
    pause(2)
    filename = 'hotspots-clusters-n-hcs';
    % print('-painters','-dpdf','-r500',filename); % pdf    
    export_fig([dir_.TomographyMaps,filename], '-pdf','-r500')
    pause(1)
    end
    
    
    hotspots = round(hcs_prevCluster);
    
    %-- moving hcs from occupied cells to the nearest unoccupied cells
    [ hotspots ] = fMovingOccupiedHcsToFreeSpace( hotspots,...
                                                  OBS,...
                                                  map_env,...
                                                  map_hotspot,...
                                                  dataYAML,...
                                                  para_ );
    hotspots = unique(hotspots,'rows');
        
    %{
    disp('Im out .sössösö')
    pause
        % plot
        % visualize = 1
        if visualize==1
            figure('name','clusters - 3rd layer'); 
            imshow(map_env','InitialMagnification',350); hold on;
            set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
            pause(1)
            % colors = lines(num_of_clusters1);
            %colors = lines(9);
            % colors = linspecer(20);
            %colors = linspecer(30);
            colors = linspecer(num_of_clusters_this);
            for i = 1:num_of_clusters_this
                plot(Hc_nCluster(clusters_this==i,1),Hc_nCluster(clusters_this==i,2),...
                     's',...
                     'color',colors(i,:),...
                     'MarkerFaceColor',colors(i,:),...
                     'MarkerSize',5);
            end

        end
    
    pause
    
    %}
    
    %{
    
    %-- find repeated hcs
    %[u,I,J] = unique(hotspots, 'rows', 'first');
    %hasDuplicates = size(u,1) < size(hotspots,1)
    %ixDupHcs = setdiff(1:size(hotspots,1), I)
    %dupHcs = hotspots(ixDupHcs,:)

    
    %-- merg all hcs close enough
    dist_ind_kaamkay = 1;    
    testind = 1;
    while numel(dist_ind_kaamkay)>=1 
        
        testind
        testind = testind+1;
        
        [hotspots,hcs_weights]
        
        size(hotspots)
        
        %hcs_dist = fPathDist2(hotspots,hotspots,map_env);
        
        
        hcs_comb_ind = nchoosek(1:size(hotspots,1),2);
        hcs_comb_dist = zeros(size(hcs_comb_ind,1),1);
        
        for i = 1:size(hcs_comb_dist,1)
            %i
            hc_from = hotspots(hcs_comb_ind(i,1),:);
            hc_to   = hotspots(hcs_comb_ind(i,2),:);
            
            dist_this = fPathDist2(hc_from,hc_to,map_env);
            hcs_comb_dist(i) = dist_this;
            
            %pause
        end
        
        hcs_comb_dist
        
        [dist_val,dist_ind] = sort(hcs_comb_dist);
        
        dist_ind_kaamkay = dist_ind(dist_val<20);
        
        for i = 1:numel(dist_ind_kaamkay)
            i
            [hotspots,hcs_weights]
            
            hc_first = hotspots(hcs_comb_ind(dist_ind_kaamkay(i),1),:)
            w_first = hcs_weights(hcs_comb_ind(dist_ind_kaamkay(i),1))
            
            hc_second = hotspots(hcs_comb_ind(dist_ind_kaamkay(i),2),:)
            w_second = hcs_weights(hcs_comb_ind(dist_ind_kaamkay(i),2))
            
            w_mirged = w_first + w_second
            
            hc_mirged_r = ((w_first*hc_first(1,1)) + (w_second*hc_second(1,1))) / w_mirged;
            hc_mirged_c = ((w_first*hc_first(1,2)) + (w_second*hc_second(1,2))) / w_mirged;
            
            hc_mirged = [hc_mirged_r,hc_mirged_c]
            
            hotspots(hcs_comb_ind(dist_ind_kaamkay(i),1),:)
            hotspots(hcs_comb_ind(dist_ind_kaamkay(i),1),:) = hc_mirged;
            hotspots(hcs_comb_ind(dist_ind_kaamkay(i),2),:)            
            hotspots(hcs_comb_ind(dist_ind_kaamkay(i),2),:) = hc_mirged;
            
            hcs_weights(hcs_comb_ind(dist_ind_kaamkay(i),1))
            hcs_weights(hcs_comb_ind(dist_ind_kaamkay(i),1)) = w_mirged;
            hcs_weights(hcs_comb_ind(dist_ind_kaamkay(i),2))
            hcs_weights(hcs_comb_ind(dist_ind_kaamkay(i),2)) = w_mirged;
            
            [hotspots,hcs_weights]
            
            
            pause
            clc
        end
        pause
        hotspots
        hotspots = round(hotspots)
        
        hcs_ws = [hotspots,hcs_weights]
        hcs_ws = unique(hcs_ws,'rows');
                       
        hotspots = hcs_ws(:,1:2)
        hcs_weights = hcs_ws(:,3)
        
        pause
    end
    
    
    disp('Im out from the loop......')
    pause
    
    
    %-- merg all hcs close enough
    testind = 1;
    toMergHcs = [0,0];    
    while numel(toMergHcs)>=2 
        
        testind
        testind = testind+1;
        size(hotspots)
        
        hcs_dist = fPathDist2(hotspots,hotspots,map_env);

        for i = 1:size(hotspots,1)
            hcs_dist(i,i) = inf;
        end
        [r,c] = ind2sub(size(hcs_dist), find(hcs_dist<10));

        toMergHcs = unique(sort([r,c],2),'rows');

        for i = 1:size(toMergHcs,1)
            hotspots(toMergHcs(i,1),1) = (hotspots(toMergHcs(i,1),1)+hotspots(toMergHcs(i,2),1))/2;
            hotspots(toMergHcs(i,1),2) = (hotspots(toMergHcs(i,1),2)+hotspots(toMergHcs(i,2),2))/2;

            hotspots(toMergHcs(i,2),1) = hotspots(toMergHcs(i,1),1);
            hotspots(toMergHcs(i,2),2) = hotspots(toMergHcs(i,1),2);

        end
        hotspots = round (hotspots);
        hotspots = unique(hotspots,'rows');
    end
    
    %}

else
    hotspots = Hcells;
    
end


 
%///////////////////////////////////////////// 
% HOTSPOTS OVER CLUSTERS
%/////////////////////////////////////////////
if visualize==1
figure('name','clusters with hotspots'); 
imshow(map_env','InitialMagnification',350); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
pause(1)
colors = lines(num_of_clusters1+num_of_clusters);
for i = 1:num_of_clusters
    plot(Hcells(hot_clusters==i,1),Hcells(hot_clusters==i,2),...
         's',...
         'color',colors(num_of_clusters1+i,:),...
         'MarkerFaceColor',colors(num_of_clusters1+i,:),...
         'MarkerSize',3.75); %5
end
% -- hotspots
for i = 1:size(hotspots,1)
    plot(hotspots(i,1),hotspots(i,2),...
         'o','color','k','MarkerFaceColor','r','MarkerSize',3); %7
end
pause(2)
filename = 'hotspots-clusters';
% print('-painters','-dpdf','-r500',filename); % pdf    
export_fig([dir_.TomographyMaps,filename], '-pdf','-r500')
pause(1)
end



%*******************************************************************
% HOTSPOTS OVER COARSE MAP
%*******************************************************************
if visualize==1
    
figure('name','coarse map with hotspots'); 
imshow(map_env','InitialMagnification',350); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])

% ------------- select color map for gdm --------------
gdm_colormap = flipud(hot(512)); % for main fig. 
max_concentration = max(map_recon(:)); % main fig.
delta = max_concentration/(size(gdm_colormap,1)-1);

[gd_row,gd_col] = size(map_recon);

for i = 1:gd_row
    for j = 1:gd_col
        % if positive concentration

        if map_recon(i,j)>=0
            
            if map_env(i,j) > 0

            % ---------------------- useful trash -------------------
            %col = gdm_colormap(round(map_recon(i,j)/delta)+1,:);

            % --- if concentration is greater than 1000 ppm, its still 1000 ppm
            %col = gdm_colormap( round(min(map_recon(i,j),3000)/delta)+1, :);
            %col = gdm_colormap( round(min(map_recon(i,j),1.78e+04)/delta)+1, :);
            %col = gdm_colormap( round(min(map_recon(i,j),inf)/delta)+1, :);
            % -------------------------------------------------------

            % ----------------------------------------------------
            % linear color code
            % ----------------------------------------------------
            %if strcmp(color_scale,'linear')
            % {
            linear_color_number = round( min(map_recon(i,j),max_concentration)/delta)+1;
            col = gdm_colormap( linear_color_number,:);                
            %}
            %end
            % ----------------------------------------------------

            %----- to avoid black borders of earch cell -----
            % {
            if sum(col) == 3
                col = col-(1e-10);
            end
            %}
            % ----- plot ------------
            %plot(j,i,'s','Color',col,'MarkerFaceColor',col,'MarkerSize',4);
            %plot(i,j,'s','Color',col,'MarkerFaceColor',col,'MarkerSize',4);
            %plot(cell_coords_x(i),cell_coords_y(j),'s','Color',col,'MarkerFaceColor',col,'MarkerSize',4);

            %plot(cell_coords_y(i)+0.5,cell_coords_x(j)+0.5,...
            %    's','Color',col,'MarkerFaceColor',col,'MarkerSize',5);

            %plot(cell_coords_x(j)+0.5,cell_coords_y(i)+0.5,...
            %    's','Color',col,'MarkerFaceColor',col,'MarkerSize',4.75);
            plot(i,j,...
                's',...
                'MarkerEdgeColor',col,...
                'MarkerFaceColor',col,...
                'MarkerSize',2.75); %4.85
            
            end
        end
    end
end   

pause(1)

for i = 1:size(hotspots,1)    
    plot(hotspots(i,1),hotspots(i,2),...
         'o','color','k','MarkerFaceColor','r','MarkerSize',3); %7
end
% pause(2)
% filename = 'hotspots-coarse-map';
% % print('-painters','-dpdf','-r500',filename); % pdf    
% export_fig([dir_.TomographyMaps,filename], '-pdf','-r500')
% pause(1)
% 
% pause
end
 

% -----------------------------------------------------------------------
% hotspots on occupied cells
% -----------------------------------------------------------------------

% {

[ hotspots ] = fMovingOccupiedHcsToFreeSpace( hotspots,...
                                              OBS,...
                                              map_env,...
                                              map_hotspot,...
                                              dataYAML,...
                                              para_ );

% hotspots_indices = sub2ind(size(map_env),hotspots(:,1),hotspots(:,2));
% occupiedHotspots = OBS.ind(ismember(OBS.ind,hotspots_indices));
% [occupiedHotspots_r,occupiedHotspots_c] = ind2sub(size(map_env),occupiedHotspots);
% 
% 
% occupiedHotspots_relocated_ind = zeros(numel(occupiedHotspots),1);
% 
% 
% for m = 1:numel(occupiedHotspots)
%     
%     % map to look for free cells
%     map_FreeCell = zeros(size(map_env));
%     map_FreeCell(occupiedHotspots(m)) = 1;
%         
%     % max radius to look for unoccupied cell
%     %radiusFreeCell_m = 2 % thickness of candidate conf circle (meters)
%     radiusFreeCell_cell = para_.radiusFreeCell_m/dataYAML.resolution;
% 
%     % grow area 
%     se = strel('disk',radiusFreeCell_cell,0); % outer limit
%     map_FreeCell = imdilate(map_FreeCell,se);
%     
%     % remove occupied cells from the region
%     map_FreeCell(OBS.ind(ismember(OBS.ind,find(map_FreeCell)))) = 0;
%     
%     % remove other high concentration cells
%     highConcentrationCells = find(map_hotspot);
%     map_FreeCell(highConcentrationCells(ismember(highConcentrationCells,find(map_FreeCell)))) = 0;
%     
%     % indices/subscripts of unoccupied (candidate) cells around
%     map_FreeCell_ind = find(map_FreeCell);
%     [map_FreeCell_r,map_FreeCell_c] = ind2sub(size(map_env),map_FreeCell_ind);
%         
%     % distance vector between occupied hotspot and unoccupied cells around
%     dist_hotspot2FreeCells = zeros(size(map_FreeCell_ind));
%     
%     for n = 1:numel(map_FreeCell_ind)
%         dist_hotspot2FreeCells(n) = sqrt(((map_FreeCell_c(n)-occupiedHotspots_c(m))^2)+...
%                                          ((map_FreeCell_r(n)-occupiedHotspots_r(m))^2));
%     end
%     
%     % nearest cell
%     [~,cell_ind] = min(dist_hotspot2FreeCells);
%     
%     % nearnest cell is new centroid hotspot.
%     map_hotspot(map_FreeCell_ind(cell_ind)) = 1;
%     
%     
%     occupiedHotspots_relocated_ind(m) = map_FreeCell_ind(cell_ind);
% end
% 
% 
% hotspots_indices(ismember(hotspots_indices,occupiedHotspots)) =...
%     occupiedHotspots_relocated_ind;
% [hotspots_r,hotspots_c] = ind2sub(size(map_env),hotspots_indices);
% 
% hotspots = [hotspots_r,hotspots_c];


if visualize==1
pause(1)
for i = 1:size(hotspots,1)
    plot(hotspots(i,1),hotspots(i,2),...
         'o','color','b','MarkerSize',7);
end
pause(2)
filename = 'hotspots-coarse-map';
% print('-painters','-dpdf','-r500',filename); % pdf    
export_fig([dir_.TomographyMaps,filename], '-pdf','-r500')
pause(1)

pause
end





end


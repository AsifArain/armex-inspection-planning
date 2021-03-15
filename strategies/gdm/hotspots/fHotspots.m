function [ hotspots ] = fHotspots( map_env,...
                                   map_recon,...
                                   SensorRange,...
                                   OBS,...
                                   para_,...
                                   cellsize_env,...
                                   dir_,...
                                   visualize ) 
%
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
         'MarkerSize',2.75);
end
pause(2)
filename = sprintf('%s%s-hotspots-positive-concentrations',...
    dir_.TomographyMaps,para_.ExperimentTitle)
% print('-painters','-dpdf','-r500',filename); % pdf
export_fig(filename, '-pdf','-r500')

pause(1)
end


if size(Hcells,1)>1
    %*******************************************************************
    % 
    %               CLUSTERING -- FIRST LAYER
    % 
    %*******************************************************************

    % - cluster linkage tree
    Z = linkage(Hcells,'single');
    % figure; dendrogram(Z);

    % - high concentration clusters
    Hclusters1 = cluster(Z,'cutoff',3.5,'criterion','distance');

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
    colors = linspecer(30);
    for i = 1:num_of_clusters1
        plot(Hcells(Hclusters1==i,1),Hcells(Hclusters1==i,2),...
             's',...
             'color',colors(i,:),...
             'MarkerFaceColor',colors(i,:),...
             'MarkerSize',2.75);
    end
    pause(2)
    filename = sprintf('%s%s-hotspots-clusters1-agglomerative',...
        dir_.TomographyMaps,para_.ExperimentTitle)

    %print('-painters','-dpdf','-r500',filename); % pdf
    export_fig(filename, '-pdf','-r500')

    pause(1)

    end

    %*******************************************************************
    % 
    %               CLUSTERING -- SECOND LAYER
    % 
    %*******************************************************************


    hot_clusters = zeros(size(Hcells,1),1);
    clust_num = 1;

    for i = 1:num_of_clusters1

        cluster1_this = [Hcells(Hclusters1==i,1),Hcells(Hclusters1==i,2)];

        dist_this = zeros(size(cluster1_this,1));
        for j = 1:size(dist_this,1)
            dist_this(j,:) = pdist2(cluster1_this(j,:),cluster1_this,'euclidean');
        end

        num_of_desired_clusters2_this = round(max(dist_this(:))/(SensorRange.cell/2));
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
    for i = 1:num_of_clusters
        plot(Hcells(hot_clusters==i,1),Hcells(hot_clusters==i,2),...
             's',...
             'color',colors(num_of_clusters1+i,:),...
             'MarkerFaceColor',colors(num_of_clusters1+i,:),...
             'MarkerSize',2.75);
    end
    pause(2)
    filename = sprintf('%s%s-hotspots-clusters2-kmeans',...
        dir_.TomographyMaps,para_.ExperimentTitle)
    %print('-painters','-dpdf','-r500',filename); % pdf
    export_fig(filename, '-pdf','-r500')

    pause(1)

    end 

    % -- hotspots are weighted mean
    Wmap = map_recon./max(map_recon(:));
    hotspots = zeros(num_of_clusters,2);
    for i = 1:num_of_clusters

        Hcluster = [Hcells(hot_clusters==i,1),Hcells(hot_clusters==i,2)];

        Hcluster_wmean_row = 0;
        Hcluster_wmean_col = 0;
        Hcluster_wmean_wgt = 0;
        for j = 1:size(Hcluster,1)
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

        
    end
    
    %-- merg all hcs close enough
    toMergHcs = [0,0];    
    while numel(toMergHcs)>=2 
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
% color_clusters = lines(num_of_clusters);
for i = 1:num_of_clusters
    plot(Hcells(hot_clusters==i,1),Hcells(hot_clusters==i,2),...
         's',...
         'color',colors(num_of_clusters1+i,:),...
         'MarkerFaceColor',colors(num_of_clusters1+i,:),...
         'MarkerSize',2.75); %5
end
% -- hotspots
for i = 1:size(hotspots,1)
    plot(hotspots(i,1),hotspots(i,2),...
         'o','color','k','MarkerFaceColor','r','MarkerSize',3); %7
end

pause(2)
filename = sprintf('%s%s-hotspots-clusters',...
    dir_.TomographyMaps,para_.ExperimentTitle)
%print('-painters','-dpdf','-r500',filename); % pdf
export_fig(filename, '-pdf','-r500')

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

pause(2)
filename = sprintf('%s%s-hotspots-coarse-map',...
    dir_.TomographyMaps,para_.ExperimentTitle)
%print('-painters','-dpdf','-r500',filename); % pdf
export_fig(filename, '-pdf','-r500')

pause
end
 

% -----------------------------------------------------------------------
% hotspots on occupied cells
% -----------------------------------------------------------------------

% {

hotspots_indices = sub2ind(size(map_env),hotspots(:,1),hotspots(:,2));
occupiedHotspots = OBS.ind(ismember(OBS.ind,hotspots_indices));
[occupiedHotspots_r,occupiedHotspots_c] = ind2sub(size(map_env),occupiedHotspots);


occupiedHotspots_relocated_ind = zeros(numel(occupiedHotspots),1);


for m = 1:numel(occupiedHotspots)
    
    % map to look for free cells
    map_FreeCell = zeros(size(map_env));
    map_FreeCell(occupiedHotspots(m)) = 1;
        
    % max radius to look for unoccupied cell
    %radiusFreeCell_m = 2 % thickness of candidate conf circle (meters)
    %radiusFreeCell_cell = para_.radiusFreeCell_m/dataYAML.resolution;
    radiusFreeCell_cell = para_.radiusFreeCell_m/cellsize_env;
    
    % grow area 
    se = strel('disk',radiusFreeCell_cell,0); % outer limit
    map_FreeCell = imdilate(map_FreeCell,se);
    
    % remove occupied cells from the region
    map_FreeCell(OBS.ind(ismember(OBS.ind,find(map_FreeCell)))) = 0;
    
    % remove other high concentration cells
    highConcentrationCells = find(map_hotspot);
    map_FreeCell(highConcentrationCells(ismember(highConcentrationCells,find(map_FreeCell)))) = 0;
    
    % indices/subscripts of unoccupied (candidate) cells around
    map_FreeCell_ind = find(map_FreeCell);
    [map_FreeCell_r,map_FreeCell_c] = ind2sub(size(map_env),map_FreeCell_ind);
        
    % distance vector between occupied hotspot and unoccupied cells around
    dist_hotspot2FreeCells = zeros(size(map_FreeCell_ind));
    
    for n = 1:numel(map_FreeCell_ind)
        dist_hotspot2FreeCells(n) = sqrt(((map_FreeCell_c(n)-occupiedHotspots_c(m))^2)+...
                                         ((map_FreeCell_r(n)-occupiedHotspots_r(m))^2));
    end
    
    % nearest cell
    [~,cell_ind] = min(dist_hotspot2FreeCells);
    
    % nearnest cell is new centroid hotspot.
    map_hotspot(map_FreeCell_ind(cell_ind)) = 1;
    
    
    occupiedHotspots_relocated_ind(m) = map_FreeCell_ind(cell_ind);
end


hotspots_indices(ismember(hotspots_indices,occupiedHotspots)) =...
    occupiedHotspots_relocated_ind;
[hotspots_r,hotspots_c] = ind2sub(size(map_env),hotspots_indices);

hotspots = [hotspots_r,hotspots_c];


if visualize==1
pause(1)
for i = 1:size(hotspots,1)
    plot(hotspots(i,1),hotspots(i,2),...
         'o','color','b','MarkerSize',7);
end
pause
end




end


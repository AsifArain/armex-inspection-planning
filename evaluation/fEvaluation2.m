function [ nearness,sc,jsd ] = fEvaluation2( recon_filename,...
                                             map_env,...
                                             para_,...
                                             dir_,...
                                             visualize )


% fEvaluation performs evaluation of reconstructions.
%
% Author:  Asif Arain
% Project: Simulator for the Evaluation of Sensing Geometries (JFR-2016)


% fprintf('\n\nEvaluation:\n');
% ______________________________________________________________________ %
%                   upload maps and reconstructions                      %
% ______________________________________________________________________ %


% mean_map = dlmread([results_dir,'mean_map.dat']);
% mean_map = dlmread([dir_.recon_mean,gdm_filename]);

% --- cell_coords_x
% cell_coords_x_file_name = strrep(gdm_filename,'mean_map','cell_coords_x');
% cell_coords_y_file_name = strrep(gdm_filename,'mean_map','cell_coords_y');

% cell_coords_x = dlmread([dir_.recon_cellx,cell_coords_x_file_name]);
% cell_coords_y = dlmread([dir_.recon_celly,cell_coords_y_file_name]);


load([dir_.recon,recon_filename],'mean_map','cell_coords_x','cell_coords_y')
map_recon = mean_map';

% map_conc_______gt = dlmread([dir_.maps,'map_conc_______gt.dat']);
% map_conc_______gt = dlmread([dir_.maps,gt_filename]);

% gt_filename = sprintf('map_conc_%d.dat',para_.GroundTruthTimeStamp);
gt_filename = sprintf('map_con.dat');
map_gt = dlmread([dir_.gt,gt_filename]);

% map_env = []
% SensorRange = []
% OBS = []
% para_ = []
% dataYAML = []
% visualize = 1
% 
% [ HighCells_rc,...
%   HighCells_ind,...
%   HighCells_clusterNums ] = fClusters( map_env,...
%                                        map_conc_______gt,...
%                                        SensorRange,...
%                                        OBS,...
%                                        para_,...
%                                        dataYAML,...
%                                        visualize );

map_coverage1 = map_gt;
map_coverage1(map_coverage1(:)>0) = 1;

% figure; imshow(map_coverage1','InitialMagnification',500); hold on;
% set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1]);

map_coverage2 = map_recon;
map_coverage2(map_coverage2(:)>0) = 1;

% figure; imshow(map_coverage2','InitialMagnification',500); hold on;
% set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1]);

map_coverage = map_coverage1+map_coverage2;
map_coverage(map_coverage(:)>0) = 1;

% figure; imshow(map_coverage','InitialMagnification',500); hold on;
% set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1]);



radius_candCircleOut = 1; %ceil(SensorRange.cell);

% grow hotspot to max radius
se = strel('disk',radius_candCircleOut,0); % outer limit
% se = strel('line',radius_candCircleOut,3); % outer limit
map_coverage = imdilate(map_coverage,se);
map_coverage(map_env(:)==0) = 0;

% figure; imshow(map_coverage','InitialMagnification',500); hold on; pause(2)
% set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1]);




% - high concentration cells
% [high_r,high_c] = find(map_conc_______gt);
% HighCells_rc = [high_r,high_c];
% HighCells_ind = find(map_conc_______gt);
coverage_ind = find(map_coverage);


% ______________________________________________________________________ %
%                           ground truth                                 %
% ______________________________________________________________________ %

% grndtruth_m = map_conc_______gt(round(cell_coords_x)+1,round(cell_coords_y));
% grndtruth_m = map_conc_______gt(round(cell_coords_x),round(cell_coords_y));
grndtruth_m = map_gt(coverage_ind);
map_recon   = map_recon(coverage_ind);

% mean_map = grndtruth_m;
% mean_map = mean_map/(max(mean_map(:)));
% mean_map = 1*ones(size(mean_map));
% mean_map(10) = 10;

% HighCells_ind1 = HighCells_ind(HighCells_clusterNums==3);
% grndtruth_m = map_conc_______gt(HighCells_ind1);
% mean_map = mean_map(HighCells_ind1);


% if visualize == 1
%     figure; 
%     subplot(1,2,1); pcolor(grndtruth_m); axis equal
%     subplot(1,2,2); pcolor(mean_map);    axis equal
% end
% 
% if visualize == 1
%     figure; 
%     subplot(1,2,1);
%     image(cell_coords_x,cell_coords_y,grndtruth_m,'CDataMapping','scaled');
%     %image(1:100,1:100,gt,'CDataMapping','scaled');hold on;        
%     % title(sprintf('Mean Map, file %s', A_FILE));
%     colormap hot;
%     %colormap gray; 
%     colormap(flipud(colormap))
%     colorbar;
%     %caxis([0 5000]);
%     %xlim([0,100]); ylim([0,100]);
%     %}
%     hold on;
%     axis equal;
% 
%     subplot(1,2,2);
%     image(cell_coords_x,cell_coords_y,mean_map,'CDataMapping','scaled');
%     %image(1:100,1:100,gt,'CDataMapping','scaled');hold on;        
%     % title(sprintf('Mean Map, file %s', A_FILE));
%     colormap hot; 
%     %colormap gray; 
%     colormap(flipud(colormap))
%     %title([num2str(table(i,2:n_c+1)),' <',...
%     %    num2str(rad2deg(normalizeAngle(deg2rad(table(i,3)-table(i,2))))),'>']);
%     colorbar;
%     %caxis([0 5000]);
%     %xlim([0,100]); ylim([0,100]);
%     %}
%     hold on;
%     axis equal;
% end


% size(grndtruth_m)
% size(mean_map)

% Display
% figure;
% subplot(2,3,1);
% pcolor(environment_model.cell_coords_x(1:end-1),environment_model.cell_coords_y(1:end-1),grndtruth_m); 
% hold on; 
% title('Ground Truth (Original Resolution)');
% set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
% 
% groundTruth.meanMap.originalResolution = grndtruth_m;


% ______________________________________________________________________ %
%                           MSE                                          %
% ______________________________________________________________________ %

mse = mean( (grndtruth_m(:)-map_recon(:)).^2 );
% fprintf('- MSE = %f\n',mse);
% new_MSE = MeanSquareError(grndtruth_m, mean_map)

% ______________________________________________________________________ %
%                           Nearness                                     %
% ______________________________________________________________________ %

nom = sum( (grndtruth_m(:)-map_recon(:)).^2 );
den = sum( ( grndtruth_m(:)- (sum(grndtruth_m(:))/numel(grndtruth_m)) ).^2 );

nearness = nom/den;
fprintf('- Nearness = %f\n',nearness);

% ______________________________________________________________________ %
%                         peak location error                            %
% ______________________________________________________________________ %

[~,num_gt] = max(grndtruth_m(:));
[gt_row,gt_col] = ind2sub(size(grndtruth_m),num_gt);

[~,num_mm] = max(map_recon(:)); 
[mm_row,mm_col] = ind2sub(size(map_recon),num_mm);


ple = sqrt((mm_row-gt_row)^2 + (mm_col-gt_col)^2);
% fprintf('- Peak location error = %f\n',ple);

% ______________________________________________________________________ %
%                         maximum difference                             %
% ______________________________________________________________________ %

md = MaximumDifference(grndtruth_m, map_recon);
% fprintf('- Maximum difference = %f\n',md);


% ______________________________________________________________________ %
%                     normalized absolute error                          %
% ______________________________________________________________________ %

nae = NormalizedAbsoluteError(grndtruth_m,map_recon);
% fprintf('- Normalized absolute error = %f\n',nae);


% ______________________________________________________________________ %
%                     normalized cross correlation                       %
% ______________________________________________________________________ %

nk = NormalizedCrossCorrelation(grndtruth_m,map_recon);
% fprintf('- Normalized cross correlation = %f\n',nk);


% ______________________________________________________________________ %
%                     Peak signal to noice ratio                         %
% ______________________________________________________________________ %

psnr = PeakSignaltoNoiseRatio(grndtruth_m,map_recon);
% fprintf('- Peak signal to noise ration = %f\n',psnr);


% ______________________________________________________________________ %
%                          Structural content                            %
% ______________________________________________________________________ %

sc = StructuralContent(grndtruth_m,map_recon);
fprintf('- Structural content = %f\n',sc);


% ______________________________________________________________________ %
%                   Structural SIMilarity (SSIM) index                   %
% ______________________________________________________________________ %

% % [mssim,~] = ssim_index(grndtruth_m,map_recon);
% 
% K = [0.05 0.05];
% % window = ones(8);
% window = ones(4);
% L = max(grndtruth_m(:));%100;
% [mssim ssim_map] = ssim_index(grndtruth_m,map_recon,K,window,L);
% fprintf('- Structural SIMilarity (SSIM) index = %f\n',mssim);
% 
% % figure; imshow(max(0, ssim_map).^4)

mssim = [];

% ______________________________________________________________________ %
%                   Histogram (probability distributions)                %
% ______________________________________________________________________ %
edges = linspace(0,max(grndtruth_m(:)),5000);
% edges = linspace(0,max([grndtruth_m(:);map_recon(:)]),5000);
% max([grndtruth_m(:);map_recon(:)])


[N_gt,edges] = histcounts(grndtruth_m,edges);
pdist_gt     = N_gt/numel(grndtruth_m);

[N_recon,edges] = histcounts(map_recon,edges);
pdist_recon     = N_recon/numel(map_recon);


% ______________________________________________________________________ %
%                        Kullback-Leibler Divergence                     %
% ______________________________________________________________________ %

[ kld ] = hcompare_KL( pdist_gt',pdist_recon');
% fprintf('- Kullback-Leibler Divergence = %f\n',kld);

% ______________________________________________________________________ %
%                        Jensen-Shannon divergence                       %
% ______________________________________________________________________ %

jsd = JSDiv( pdist_gt,pdist_recon );
fprintf('- Jensen-Shannon divergence = %f\n',jsd);

% ______________________________________________________________________ %
%         Pairwise distance between two observations - euclidean         %
% ______________________________________________________________________ %

d_eucl = pdist2(pdist_gt,pdist_recon,'euclidean');
% fprintf('- Pairwise euclidean dist = %f\n',d_eucl);


% ______________________________________________________________________ %
%         Pairwise distance between two observations - correlation       %
% ______________________________________________________________________ %

d_corr = pdist2(pdist_gt,pdist_recon,'correlation');
% fprintf('- Pairwise correlation dist = %f\n',d_corr);





% pause
% close all

end

function [ HighCells_rc,...
           HighCells_ind,...
           HighCells_clusterNums ] = fClusters( map_env,...
                                                map_conc_______gt,...
                                                SensorRange,...
                                                OBS,...
                                                para_,...
                                                dataYAML,...
                                                visualize )
%
%


map_conc_______gt(map_conc_______gt(:)>0) = 1;

% - high concentration cells
[high_r,high_c] = find(map_conc_______gt);

HighCells_rc = [high_r,high_c];
HighCells_ind = find(map_conc_______gt);

%*******************************************************************
% 
%               CLUSTERING -- FIRST LAYER
% 
%*******************************************************************

% - cluster linkage tree
Z = linkage(HighCells_rc,'single');
% figure; dendrogram(Z);

% - high concentration clusters
HighCells_clusterNums = cluster(Z,'cutoff',2,'criterion','distance');

% Hclusters = cluster(Z,'cutoff',1);
% Hclusters = cluster(Z,'cutoff',2.5,'Depth',2);
% Hclusters = cluster(Z,'cutoff',1,'MaxClust',5);

% - number of clusters
num_of_clusters1 = numel(unique(HighCells_clusterNums));



% plot
if visualize==1
figure('name','clusters - 1st layer'); 
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
pause(1)
colors = lines(num_of_clusters1);
for i = 1:num_of_clusters1
    plot(HighCells_rc(HighCells_clusterNums==i,1),HighCells_rc(HighCells_clusterNums==i,2),...
         's',...
         'color',colors(i,:),...
         'MarkerFaceColor',colors(i,:),...
         'MarkerSize',5);
end
end



%*******************************************************************
% 
%               CLUSTERING -- SECOND LAYER
% 
%*******************************************************************
% 
% 
% hot_clusters = zeros(size(HighCells_rc,1),1);
% clust_num = 0;
% 
% for i = 1:num_of_clusters1
%     
%     cluster1_this = [HighCells_rc(HighCells_clusterNums==i,1),HighCells_rc(HighCells_clusterNums==i,2)];
%     
%     dist_this = zeros(size(cluster1_this,1));
%     for j = 1:size(dist_this,1)
%         dist_this(j,:) = pdist2(cluster1_this(j,:),cluster1_this,'euclidean');
%     end
%     
%     num_of_clusters2_this = round(max(dist_this(:))/(SensorRange.cell/2));
%     
%     if num_of_clusters2_this<1
%         num_of_clusters2_this = 1;
%     end
%     
%     %cluster_tree2_this = linkage(cluster1_this,'single');    
%     %Hclusters2_this = cluster(cluster_tree2_this,'MaxClust',num_of_clusters2_this)
%     
%     [Hclusters2_this,C] = kmeans(cluster1_this,num_of_clusters2_this,'Distance','cityblock');
%     %[Hclusters2_this,C] = kmeans(cluster1_this,num_of_clusters2_this,'Distance','cosine');
%     %[Hclusters2_this,C] = kmeans(cluster1_this,num_of_clusters2_this,'Distance','correlation');
%     %[Hclusters2_this,C] = kmeans(cluster1_this,num_of_clusters2_this,'Distance','hamming'); % beykar
%     
%     
%     % - number of clusters
%     num_of_clusters2_this = numel(unique(Hclusters2_this));
%     
%     hot_clusters(HighCells_clusterNums==i) = Hclusters2_this+clust_num;
%     
%     clust_num = clust_num+num_of_clusters2_this;
%         
% end
% 
% 
% % - number of clusters
% num_of_clusters = numel(unique(hot_clusters));
% 
% 
% 
% % plot
% if visualize==1
% figure('name','clusters - 2nd layer'); 
% imshow(map_env','InitialMagnification',800); hold on;
% set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
% pause(1)
% colors = lines(num_of_clusters);
% for i = 1:num_of_clusters
%     plot(HighCells_rc(hot_clusters==i,1),HighCells_rc(hot_clusters==i,2),...
%          's',...
%          'color',colors(i,:),...
%          'MarkerFaceColor',colors(i,:),...
%          'MarkerSize',5);
% end
% end 
% 
% % -- hotspots are weighted mean
% Wmap = map_recon./max(map_recon(:));
% hotspots = zeros(num_of_clusters,2);
% for i = 1:num_of_clusters
%     
%     Hcluster = [HighCells_rc(hot_clusters==i,1),HighCells_rc(hot_clusters==i,2)];
%     
%     Hcluster_wmean_row = 0;
%     Hcluster_wmean_col = 0;
%     Hcluster_wmean_wgt = 0;
%     for j = 1:size(Hcluster,1)
%         Hcluster_wmean_row = Hcluster_wmean_row + ...
%             Wmap(Hcluster(j,1),Hcluster(j,2))*Hcluster(j,1);
%         Hcluster_wmean_col = Hcluster_wmean_col + ...
%             Wmap(Hcluster(j,1),Hcluster(j,2))*Hcluster(j,2);
%         Hcluster_wmean_wgt = Hcluster_wmean_wgt + ...
%             Wmap(Hcluster(j,1),Hcluster(j,2));
%     end
%     
%     
%     Hcluster_wmean = [Hcluster_wmean_row/Hcluster_wmean_wgt,...
%                       Hcluster_wmean_col/Hcluster_wmean_wgt];
%     
%     hotspots(i,:) = round(Hcluster_wmean);
%     
% end
% 
% if visualize==1
% pause(1)
% for i = 1:size(hotspots,1)
%     plot(hotspots(i,1),hotspots(i,2),...
%          'o','color','k','MarkerFaceColor','r','MarkerSize',7);
% end
% end
% 
% 
% % plot
% if visualize==1
% figure('name','coarse map with hotspots'); 
% imshow(map_env','InitialMagnification',800); hold on;
% set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
% 
% pause(1)
% colors = lines(num_of_clusters);
% for i = 1:num_of_clusters
%     plot(HighCells_rc(hot_clusters==i,1),HighCells_rc(hot_clusters==i,2),...
%          's',...
%          'color',colors(i,:),...
%          'MarkerFaceColor',colors(i,:),...
%          'MarkerSize',5);
% end
% end 
% 
% % plot
% if visualize==1
%     
% figure('name','coarse map with hotspots'); 
% imshow(map_env','InitialMagnification',800); hold on;
% set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
% 
% % ------------- select color map for gdm --------------
% gdm_colormap = flipud(colormap('hot')); % for main fig. 
% max_concentration = max(map_recon(:)); % main fig.
% delta = max_concentration/(size(gdm_colormap,1)-1);
% 
% [gd_row,gd_col] = size(map_recon);
% 
% for i = 1:gd_row
%     for j = 1:gd_col
%         % if positive concentration
% 
%         if map_recon(i,j)>=0
% 
%             % ---------------------- useful trash -------------------
%             %col = gdm_colormap(round(map_recon(i,j)/delta)+1,:);
% 
%             % --- if concentration is greater than 1000 ppm, its still 1000 ppm
%             %col = gdm_colormap( round(min(map_recon(i,j),3000)/delta)+1, :);
%             %col = gdm_colormap( round(min(map_recon(i,j),1.78e+04)/delta)+1, :);
%             %col = gdm_colormap( round(min(map_recon(i,j),inf)/delta)+1, :);
%             % -------------------------------------------------------
% 
%             % ----------------------------------------------------
%             % linear color code
%             % ----------------------------------------------------
%             %if strcmp(color_scale,'linear')
%             % {
%             linear_color_number = round( min(map_recon(i,j),max_concentration)/delta)+1;
%             col = gdm_colormap( linear_color_number,:);                
%             %}
%             %end
%             % ----------------------------------------------------
% 
%             %----- to avoid black borders of earch cell -----
%             % {
%             if sum(col) == 3
%                 col = col-(1e-10);
%             end
%             %}
%             % ----- plot ------------
%             %plot(j,i,'s','Color',col,'MarkerFaceColor',col,'MarkerSize',4);
%             %plot(i,j,'s','Color',col,'MarkerFaceColor',col,'MarkerSize',4);
%             %plot(cell_coords_x(i),cell_coords_y(j),'s','Color',col,'MarkerFaceColor',col,'MarkerSize',4);
% 
%             %plot(cell_coords_y(i)+0.5,cell_coords_x(j)+0.5,...
%             %    's','Color',col,'MarkerFaceColor',col,'MarkerSize',5);
% 
%             %plot(cell_coords_x(j)+0.5,cell_coords_y(i)+0.5,...
%             %    's','Color',col,'MarkerFaceColor',col,'MarkerSize',4.75);
%             plot(i,j,...
%                 's',...
%                 'MarkerEdgeColor',col,...
%                 'MarkerFaceColor',col,...
%                 'MarkerSize',4.85);
%         end
%     end
% end   
% 
% pause(1)
% % colors = lines(num_of_clusters);
% % for i = 1:num_of_clusters
% %     plot(Hcells(hot_clusters==i,1),Hcells(hot_clusters==i,2),...
% %          's',...
% %          'color',colors(i,:),...
% %          'MarkerFaceColor',colors(i,:),...
% %          'MarkerSize',5);
% %      %pause
% % end
% 
% for i = 1:size(hotspots,1)    
%     plot(hotspots(i,1),hotspots(i,2),...
%          'o','color','k','MarkerFaceColor','r','MarkerSize',7);
% end
% end
%  
% % -----------------------------------------------------------------------
% % hotspots on occupied cells
% % -----------------------------------------------------------------------
% 
% % {
% 
% hotspots_indices = sub2ind(size(map_env),hotspots(:,1),hotspots(:,2));
% occupiedHotspots = OBS.ind(ismember(OBS.ind,hotspots_indices));
% [occupiedHotspots_r,occupiedHotspots_c] = ind2sub(size(map_env),occupiedHotspots);
% 
% 
% occupiedHotspots_relocated_ind = zeros(numel(occupiedHotspots),1);
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
% 
% 
% if visualize==1
% pause(1)
% for i = 1:size(hotspots,1)
%     plot(hotspots(i,1),hotspots(i,2),...
%          'o','color','b','MarkerSize',7);
% end
end




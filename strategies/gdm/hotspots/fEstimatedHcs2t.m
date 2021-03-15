function [ hotspots ] = fEstimatedHcs2t( map_env,...
                                         map_recon,...
                                         SensorRange,...
                                         OBS,...
                                         cellsize_env,...
                                         para_,...
                                         dir_,...
                                         visualize ) 
%
%


% %------------------------------------
% %    PUBLISH PARAMETERS (PAPER)
% %------------------------------------
% paraPub.ZoomIn = 750;
% paraPub.CellMarkerSize = 10;%7.0;
% paraPub.BoundryLineWidth = 1;
% paraPub.BoundryLineColor = [0.3,0.3,0.3];
% paraPub.HcMarkerSize = 10;
% paraPub.HcLineWidth = 2;
% paraPub.UpperPPM = 20000; %30000
% %--------------------------------

%------------------------------------
%    PUBLISH PARAMETERS (THESIS)
%------------------------------------
paraPub.ZoomIn = 350;
paraPub.CellMarkerSize = 2.95 %10;%7.0;
paraPub.BoundryLineWidth = 1;
paraPub.BoundryLineColor = [0.3,0.3,0.3];
paraPub.HcMarkerSize = 4 %10;
paraPub.HcLineWidth = 1;
paraPub.UpperPPM = 20000; %30000

%--------------------------------

map_hotspot = zeros(size(map_recon));
% ind = find(gd_map_hotspot(:))>500;
% gd_map_hotspot(ind) = 1;
map_hotspot((map_recon(:)>para_.highConcentrationThreshold_PPM)) = 1;
% test_map_resized2_obs = find(test_map_resized2(:)==0);
% gd_map_coverage(test_map_resized2_obs) = 0;



%------------------------------------





% warning("Hardcore change for sample-env-thesis");
% map_recon(33,16) = 0;
% para_.highConcentrationThreshold_PPM




%------------------------------------



% -- hotspots based on clustering

% - high concentration cells
[high_r,high_c] = find(map_hotspot);
Hcells = [high_r,high_c];


%-------------------------------------
% plot: positive concentration
%-------------------------------------
% visualize = 1
if visualize==1
    figure('name','step 1: Non-zero concentrations'); 
    imshow(map_env','InitialMagnification',paraPub.ZoomIn); hold on; %350
    set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
    %color_positive_con = [255,140,0]/255;
    color_positive_con = [218,165,32]/255;
    pause(1)
    for i = 1:size(Hcells,1)
        plot(Hcells(i,1)+0.0,Hcells(i,2)+0.0,...
             's',...
             'color',color_positive_con,...
             'MarkerFaceColor',color_positive_con,...
             'MarkerSize',paraPub.CellMarkerSize); %2.75 %3.75
    end
    pause(2)
    %filename = 'hotspots-positive-concentrations';
    % print('-painters','-dpdf','-r500',filename); % pdf    
    %export_fig([dir_.TomographyMaps,filename], '-pdf','-r500')
    filename = sprintf('%s%s-hotspots-step1-positive-concentrations',...
        dir_.gdmplan_hotspots,para_.FilePrefix);
    export_fig(filename, '-pdf','-r500');
    pause(1)
    %pause
    
    %--- backup ---
    save([dir_.gdmplan_hotspots,'plot_step1_positive_concentration.mat'],...
        'Hcells');
end


%*******************************************************************
% 
%               CLUSTERING
% 
%*******************************************************************

%-- linkage cutoff for clusters
cut_off = 2.5 %3.5; %median([Z(end-2,3) Z(end-1,3)])



%--- filter outliers
%---------------------------------------------
% if size(Hcells,1)>1
%     
%     
%     %-- cluster linkage tree
%     Z = linkage(Hcells,'single');
%     %figure; dendrogram(Z,'ColorThreshold',cut_off)
% 
%     % - high concentration clusters
%     clusters2Filter = cluster(Z,'cutoff',cut_off,'criterion','distance');
%         
%     % - number of clusters
%     num_of_clusters = numel(unique(clusters2Filter));
%     
%     %-- weighted map
%     Wmap = map_recon./sum(map_recon(:)); % total weight is 1    
%     
%     %-- filter
%     filter_ind = [];
%     for i = 1:num_of_clusters
%         
%         %-- this cluster
%         thisCluster = Hcells(clusters2Filter==i,:);
%         
%         %-- if number of observation are less than a number
%         if size(thisCluster,1) < 5
%             
%             %-- concentration weight of this cluster
%             thisWeight = 0;
%             for j = 1:size(thisCluster,1)
%                 thisWeight = thisWeight + Wmap(thisCluster(j,1),thisCluster(j,2));
%             end
%             %thisWeight
%             %-- if concentration weight of this cluster is less than a
%             %number
%             if thisWeight < 0.005 %0.01
%                 %-- filter out the observations (cells)
%                 %Hcells_filtered(clusters2Filter==i,:) = [];
%                 
%                 filter_ind = [filter_ind;find(clusters2Filter==i)];
%                 
%             end
%         end
%     end
%     
%     Hcells(filter_ind,:) = [];
%     
%     
%     
%     % plot
%     if visualize==1
%         figure('name','Filtered observations (cells)'); 
%         imshow(map_env','InitialMagnification',350); hold on;
%         set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
%         pause(1)
%         for i = 1:size(Hcells,1)
%             plot(Hcells(i,1),Hcells(i,2),...
%                  's',...
%                  'color','k',...%color_positive_con,...
%                  'MarkerFaceColor',color_positive_con,...
%                  'MarkerSize',3.75); %2.75
%         end
%         pause(2)
%         filename = 'hotspots-filtered-cells';
%         % print('-painters','-dpdf','-r500',filename); % pdf    
%         export_fig([dir_.TomographyMaps,filename], '-pdf','-r500')
%         pause(1)
%     end
%         
%     pause
% end




if size(Hcells,1)>1
    
    
    %*******************************************
    %               FIRST LAYER
    %*******************************************
    
    %*******************************************
    % Step-2: CLUSTERING
    % clusters based on distance
    %*******************************************

    % - cluster linkage tree
    Z = linkage(Hcells,'single');
    %figure; dendrogram(Z);
    
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


    %-------------------------------------
    % plot: clusters 1st layer
    %-------------------------------------
    % visualize = 1
    
    
    
    
    if visualize==1
        figure('name','clusters - 1st layer'); 
        imshow(map_env','InitialMagnification',paraPub.ZoomIn); hold on; %350
        set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
        pause(1)
        % colors = lines(num_of_clusters1);
        %colors = lines(9);
        % colors = linspecer(20);
        %colors = linspecer(30);
        colors = linspecer(num_of_clusters1);
        for i = 1:num_of_clusters1
            plot(Hcells(Hclusters1==i,1)+0.0,Hcells(Hclusters1==i,2)+0.0,...
                 's',...
                 'color',colors(i,:),...
                 'MarkerFaceColor',colors(i,:),...
                 'MarkerSize',paraPub.CellMarkerSize); %3.75
        end    
        for i = 1:num_of_clusters1
             %-- boundry
             boundycandidates = Hcells(Hclusters1==i,:);
             boundypoints = boundary(boundycandidates,1);
             plot(boundycandidates(boundypoints,1),...
                  boundycandidates(boundypoints,2),...
                  '-',...
                  'color',paraPub.BoundryLineColor,...
                  'LineWidth',paraPub.BoundryLineWidth);
        end    
        pause(2)
        %filename = 'hotspots-clusters1-agglomerative';
        % print('-painters','-dpdf','-r500',filename); % pdf    
        %export_fig([dir_.TomographyMaps,filename], '-pdf','-r500')        
        filename = sprintf('%s%s-hotspots-step2-clusters1-agglomerative',...
            dir_.gdmplan_hotspots,para_.FilePrefix);
        export_fig(filename, '-pdf','-r500');    
        pause(1)
        %pause
        
        %--- backup ---
        save([dir_.gdmplan_hotspots,'plot_step2_cluster1.mat'],...
            'Hcells','Hclusters1');
    end

    %*******************************************
    %           SECOND LAYER
    %*******************************************
    
    %*******************************************
    % Step-3: CLUSTERING
    % sub-divide lengthy clusters
    %*******************************************

    hot_clusters = zeros(size(Hcells,1),1);
    clust_num = 1;

    for i = 1:num_of_clusters1

        cluster1_this = [Hcells(Hclusters1==i,1),Hcells(Hclusters1==i,2)];

        dist_this = zeros(size(cluster1_this,1));
        for j = 1:size(dist_this,1)
            dist_this(j,:) = pdist2(cluster1_this(j,:),cluster1_this,'euclidean');
        end

        %num_of_desired_clusters2_this = ...
        %    round(max(dist_this(:))/(SensorRange.cell/2));
        num_of_desired_clusters2_this = ...
            ceil(max(dist_this(:))/(SensorRange.cell/2));
        if num_of_desired_clusters2_this<1
            num_of_desired_clusters2_this = 1;
        end


        if num_of_desired_clusters2_this>1

            %if num_of_clusters2_this<1
            %    num_of_clusters2_this = 1;
            %end
            
            distance_type = 'cityblock';
            %distance_type = 'cosine';
            %distance_type = 'correlation';
            %distance_type = 'hamming'; % beykar            
            [Hclusters2_this,~] = kmeans( cluster1_this,...
                                          num_of_desired_clusters2_this,...
                                          'Distance',distance_type );

            % - number of clusters
            num_of_desired_clusters2_this = numel(unique(Hclusters2_this));
            hot_clusters(Hclusters1==i) = Hclusters2_this+clust_num-1;
            clust_num = clust_num+num_of_desired_clusters2_this;

        else
            hot_clusters(Hclusters1==i) = clust_num;
            clust_num = clust_num+1;
        end
    end


    % - number of clusters
    num_of_clusters = numel(unique(hot_clusters));
    % - number of clusters
    num_of_kmean_clusters = numel(unique(hot_clusters));



    % plot
    if visualize==1
        figure('name','clusters - 2nd layer'); 
        imshow(map_env','InitialMagnification',paraPub.ZoomIn); hold on;
        set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
        pause(1)
        %-- clusters
        % colors = lines(num_of_clusters);
        colors = linspecer(num_of_clusters1+num_of_clusters);
        for i = 1:num_of_clusters
            plot(Hcells(hot_clusters==i,1)+0.0,Hcells(hot_clusters==i,2)+0.0,...
                 's',...
                 'color',colors(num_of_clusters1+i,:),...
                 'MarkerFaceColor',colors(num_of_clusters1+i,:),...
                 'MarkerSize',paraPub.CellMarkerSize);
        end
        %-- boundries
        for i = 1:num_of_clusters
             %-- boundry
             boundycandidates = Hcells(hot_clusters==i,:);
             boundypoints = boundary(boundycandidates,1);
             plot(boundycandidates(boundypoints,1),...
                  boundycandidates(boundypoints,2),...
                  '-',...
                  'color',paraPub.BoundryLineColor,...
                  'LineWidth',paraPub.BoundryLineWidth);
        end        
        pause(2)
        %filename = 'hotspots-clusters2-kmeans';
        %%% print('-painters','-dpdf','-r500',filename); % pdf    
        %export_fig([dir_.TomographyMaps,filename], '-pdf','-r500')        
        filename = sprintf('%s%s-hotspots-step3-clusters2-kmeans',...
            dir_.gdmplan_hotspots,para_.FilePrefix);
        export_fig(filename, '-pdf','-r500');
        pause(1)
        %pause
        
        %--- backup ---
        save([dir_.gdmplan_hotspots,'plot_step3_cluster2.mat'],...
            'Hcells','hot_clusters');
    end 

    
    %************************************************
    % Step-4: ESTIMATE HCs
    % hotspot centers are weighted mean in each cluster
    %************************************************
    
    % -- hotspots are weighted mean
    Wmap = map_recon./sum(map_recon(:)); % total weight is 1
    hotspots = zeros(num_of_clusters,2);
    hcs_weights = zeros(num_of_clusters,1);
    hcs_ppmm = zeros(num_of_clusters,1);
    for i = 1:num_of_clusters

        Hcluster = [Hcells(hot_clusters==i,1),Hcells(hot_clusters==i,2)];
        
        Hcluster_wmean_row = 0;
        Hcluster_wmean_col = 0;
        Hcluster_wmean_wgt = 0;
        Hcluster_ppmm      = 0;
        for j = 1:size(Hcluster,1)
            %
            Hcluster_wmean_row = Hcluster_wmean_row + ...
                Wmap(Hcluster(j,1),Hcluster(j,2))*Hcluster(j,1);
            
            Hcluster_wmean_col = Hcluster_wmean_col + ...
                Wmap(Hcluster(j,1),Hcluster(j,2))*Hcluster(j,2);
            
            Hcluster_wmean_wgt = Hcluster_wmean_wgt + ...
                Wmap(Hcluster(j,1),Hcluster(j,2));
            
            Hcluster_ppmm = Hcluster_ppmm +...
                map_recon(Hcluster(j,1),Hcluster(j,2));
        end

        Hcluster_wmean = [Hcluster_wmean_row/Hcluster_wmean_wgt,...
                          Hcluster_wmean_col/Hcluster_wmean_wgt];

        hotspots(i,:) = round(Hcluster_wmean);
        hcs_weights(i) = Hcluster_wmean_wgt;
        hcs_ppmm(i) = Hcluster_ppmm;
        
    end
    
    %hcw_before = [hotspots,hcs_weights];
    
    %************************************************    
    % move estimated centers to unoccupied space
    %************************************************        
    [ hotspots ] = fMovingOccupiedHcsToFreeSpace( hotspots,...
                                                  OBS,...
                                                  map_env,...
                                                  cellsize_env,...
                                                  para_ );
    %hcw_after = [hotspots,hcs_weights];
    
    %find(hcw_after(:)~=hcw_before(:))
    %pause
    
    if visualize==1
        
        figure('name','step 4: clusters-2 with hotspots'); 
        imshow(map_env','InitialMagnification',paraPub.ZoomIn); hold on;
        set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
        pause(1)
        
        %-- coarse map
        fPublishRecon(map_recon,map_env,paraPub);
        pause(1)

        colors = linspecer(num_of_clusters1+num_of_clusters);        
        %colors = lines(num_of_clusters);
        for i = 1:num_of_clusters
            %plot(Hcells(hot_clusters==i,1)+0.0,Hcells(hot_clusters==i,2)+0.0,...
            %     's',...
            %     'color',colors(i,:),...
            %     'MarkerFaceColor',colors(i,:),...
            %     'MarkerSize',paraPub.CellMarkerSize); %5
             
             %-- boundry
             boundycandidates = Hcells(hot_clusters==i,:);
             boundypoints = boundary(boundycandidates,1);
             plot(boundycandidates(boundypoints,1),...
                  boundycandidates(boundypoints,2),...
                  '-',...
                  'color',paraPub.BoundryLineColor,...
                  'LineWidth',paraPub.BoundryLineWidth);
        end        
        % -- hotspots
        for i = 1:size(hotspots,1)
            plot(hotspots(i,1),hotspots(i,2),...
                 'o','color','k',...
                 'MarkerFaceColor','r',...
                 'LineWidth',paraPub.HcLineWidth,...
                 'MarkerSize',paraPub.HcMarkerSize); %7
        end
        pause(2)
        %filename = 'hotspots-clusters2-kmeans-hc';
        %%% print('-painters','-dpdf','-r500',filename); % pdf    
        %export_fig([dir_.TomographyMaps,filename], '-pdf','-r500')        
        filename = sprintf('%s%s-hotspots-step4-clusters2-kmeans-hc',...
            dir_.gdmplan_hotspots,para_.FilePrefix);
        export_fig(filename, '-pdf','-r500');        
        pause(1)
        %pause
        
        %--- backup ---
        save([dir_.gdmplan_hotspots,'plot_step4_cluster2.mat'],...
            'Hcells','hot_clusters');
    end
    
    %pause
    %************************************************
    % Step-5: ESTIMATE CENTERS
    % nearby centers are combined
    %************************************************   
    
    warning('with or without hotspot fusion? select or deselect here');
    %{
    warning('lets go with hotspot fusion');
    
    %-- merg all hcs close enough
    %dist_ind_kaamkay = 1;
    % -- temp fix    
    warning('temporary fix.');
    if size(hotspots,1) < 2
        dist_ind_kaamkay = []
    else 
        dist_ind_kaamkay = 1
    end
    
    testind = 0;
    while numel(dist_ind_kaamkay)>=1

        testind = testind+1;
        %hcs_dist = fPathDist2(hotspots,hotspots,map_env);
        
        %-- combinations and pairwise distance
        %-----------------------------------------------
        scs_comb_ind = nchoosek(1:size(hotspots,1),2);
        scs_comb_dist = inf(size(scs_comb_ind,1),1);
        for i = 1:size(scs_comb_dist,1)
            %-- from -- to
            sc_from = hotspots(scs_comb_ind(i,1),:);
            sc_to   = hotspots(scs_comb_ind(i,2),:);
            %-- traveling dist
            %dist_this = fPathDist2(sc_from,sc_to,map_env);
            dist_this = pdist2(sc_from,sc_to,'euclidean');
            scs_comb_dist(i) = dist_this;
            %if dist_this<=para_.SCsCombiningDist_cells
            %    break;
            %end
        end
        
        [dist_val,dist_ind] = sort(scs_comb_dist);
        dist_ind_kaamkay = dist_ind(dist_val<=para_.SCsCombiningDist_cells);
        
        if numel(dist_ind_kaamkay)>=1
            
            sc_first = hotspots(scs_comb_ind(dist_ind_kaamkay(1),1),:);
            w_first  = hcs_weights(scs_comb_ind(dist_ind_kaamkay(1),1));
            ppm_first = hcs_ppmm(scs_comb_ind(dist_ind_kaamkay(1),1));
            
            sc_second = hotspots(scs_comb_ind(dist_ind_kaamkay(1),2),:);
            w_second  = hcs_weights(scs_comb_ind(dist_ind_kaamkay(1),2));
            ppm_second = hcs_ppmm(scs_comb_ind(dist_ind_kaamkay(1),2));
            
            w_mirged = w_first + w_second;
            
            sc_mirged_r = ((w_first*sc_first(1,1)) + (w_second*sc_second(1,1))) / w_mirged;
            sc_mirged_c = ((w_first*sc_first(1,2)) + (w_second*sc_second(1,2))) / w_mirged;
            
            sc_mirged = [sc_mirged_r,sc_mirged_c];
            
            ppm_mirged = ppm_first + ppm_second;
            
            hotspots(scs_comb_ind(dist_ind_kaamkay(1),1),:) = sc_mirged;
            hotspots(scs_comb_ind(dist_ind_kaamkay(1),2),:) = [];
            
            hcs_weights(scs_comb_ind(dist_ind_kaamkay(1),1)) = w_mirged;
            hcs_weights(scs_comb_ind(dist_ind_kaamkay(1),2)) = [];
            
            hcs_ppmm(scs_comb_ind(dist_ind_kaamkay(1),1)) = ppm_mirged;
            hcs_ppmm(scs_comb_ind(dist_ind_kaamkay(1),2)) = [];
            
        end
        
        %-- tmp fix
        warning('temporary fix.');
        if size(hotspots,1) < 2
            dist_ind_kaamkay = [];
        end
    end
    
    
    %----------
    % plot
    %----------
    if visualize==1
        figure('name','Estimated hotspots - step5 combined centers');
        imshow(map_env','InitialMagnification',paraPub.ZoomIn); hold on;
        set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
        pause(1)
        
        %-- coarse map
        fPublishRecon(map_recon,map_env,paraPub);
        pause(1)

        %-- estimated sources
        sc_estimated_plot = sortrows(hotspots,1);
        for i = 1:size(hotspots,1)
            plot(sc_estimated_plot(i,1),sc_estimated_plot(i,2),...
                 'o',...
                 'color','k',...
                 'MarkerFaceColor','r',...
                 'LineWidth',paraPub.HcLineWidth,...
                 'MarkerSize',paraPub.HcMarkerSize); %7 %5
             %- estimated source num
             %text(sc_estimated_plot(i,1)+1.0,...
             %     sc_estimated_plot(i,2),num2str(i),...
             %     'FontSize',6)
        end
        pause(2)
        %-- save
        %filename = sprintf('%s%s-estimated-sources-clustering-step5',...
        %    dir_.Evaluation,para_.FilePrefix);
        %filename = sprintf('%s%s-estimated-hcs-clustering-step5',...
        %    dir_.TomographyMaps,para_.FilePrefix);
        filename = sprintf('%s%s-hotspots-step5-combined-centers',...
            dir_.gdmplan_hotspots,para_.FilePrefix);
        %print('-painters','-dpdf','-r500',filename); % pdf
                %filename = 'hotspots-clusters2-kmeans-hc';
        % print('-painters','-dpdf','-r500',filename); % pdf
        export_fig(filename, '-pdf','-r500');        
        pause(1)
        
        %--- backup ---
        save([dir_.gdmplan_hotspots,'plot_step5_combine.mat'],...
            'hotspots');
    end
    
    %}
    
    
    %*******************************************
    %           nth LAYER
    %*******************************************
    
    hcs_prevCluster = hotspots;
    ws_prevCluster  = hcs_weights;
    ppmm_prevCluster  = hcs_ppmm;
    
    %{
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
        ppmm_thisCluster  = zeros(num_of_clusters_this,1);
        
        for i = 1:num_of_clusters_this

            candidateCells = [hcs_prevCluster(clusters_this==i,1),...
                              hcs_prevCluster(clusters_this==i,2)];
            candiateWeights = ws_prevCluster(clusters_this==i);            
            candiatePPMM = ppmm_prevCluster(clusters_this==i);
            
            thisHc_r = 0;
            thisHc_c = 0;
            thisHc_w = 0;
            thisHc_ppmm = 0;
            for j = 1:size(candidateCells,1)
                %
                thisHc_r = thisHc_r + candiateWeights(j)*candidateCells(j,1);
                thisHc_c = thisHc_c + candiateWeights(j)*candidateCells(j,2);
                thisHc_w = thisHc_w + candiateWeights(j);
                thisHc_ppmm = thisHc_ppmm + candiatePPMM(j);
            end
            
            thisHc = [thisHc_r/thisHc_w,thisHc_c/thisHc_w];

            hcs_thisCluster(i,:) = round(thisHc);
            ws_thisCluster(i) = thisHc_w;
            ppmm_thisCluster(i) = thisHc_ppmm;

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
            pause
        end
        
        cond_nClustering = size(hcs_thisCluster,1)>1 & ...
                           (size(hcs_prevCluster,1)-size(hcs_thisCluster,1))>0;
        
        hcs_prevCluster = hcs_thisCluster;
        ws_prevCluster = ws_thisCluster;
        ppmm_prevCluster = ppmm_thisCluster;
        
        %pause
        
    end
    
    if visualize==1
        pause(2)
        filename = 'hotspots-clusters-n-hcs';
        % print('-painters','-dpdf','-r500',filename); % pdf    
        export_fig([dir_.TomographyMaps,filename], '-pdf','-r500')
        pause(1)
        pause
    end
    
    %}
    
    %hotspots = round(hcs_prevCluster);
    hotspots = round(hcs_prevCluster*100)/100;
    
    
    
    %*************************************************
    % Step-6: WEIGHT CONDITION
    %-- high pass filtering of hcs based on weights
    %*************************************************
    total_weight = sum(ws_prevCluster);
    hcs_weights = ws_prevCluster;
    hcs_ppmm = ppmm_prevCluster;
    num_of_hcs = size(hotspots,1);
    
    %qualification_threshould = (total_weight/num_of_hcs*0.5)
    qualification_w = (total_weight/num_of_hcs*para_.HcsCutoffWeightFactor);    
        
    %hotspots = hotspots(hcs_weights>qualification_threshould,:)
    %weight_condition = find(hcs_weights>qualification_w);
    %concentration_condition = find(hcs_ppmm>para_.HcsCutoffPPMM);
    
    hcs_weights
    hcs_ppmm
    
    qualification_w
    
    hotspots = hotspots(hcs_weights>qualification_w | hcs_ppmm>para_.HcsCutoffPPMM,:);
        
    %-- moving hcs from occupied cells to the nearest unoccupied cells
    [ hotspots ] = fMovingOccupiedHcsToFreeSpace( hotspots,...
                                                  OBS,...
                                                  map_env,...
                                                  cellsize_env,...
                                                  para_ );
    
    hotspots = unique(hotspots,'rows');
    
    %----------
    % plot
    %----------
    if visualize==1
        figure('name','Estimated hotspots - step6 - weight condition');
        imshow(map_env','InitialMagnification',paraPub.ZoomIn); hold on;
        set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
        pause(1)
        %-- coarse map
        fPublishRecon(map_recon,map_env,paraPub);
        pause(1)
        %-- estimated sources
        sc_estimated_plot = sortrows(hotspots,1);
        for i = 1:size(hotspots,1)
            plot(sc_estimated_plot(i,1),sc_estimated_plot(i,2),...
                 'o',...
                 'color','k',...
                 'MarkerFaceColor','r',...
                 'LineWidth',paraPub.HcLineWidth,...
                 'MarkerSize',paraPub.HcMarkerSize);
        end
        pause(2)
        %-- save
        filename = sprintf('%s%s-hotspots-step6-weight-condition',...
            dir_.gdmplan_hotspots,para_.FilePrefix);
        export_fig(filename, '-pdf','-r500');
        pause(1)
        
        %--- backup
        save([dir_.gdmplan_hotspots,'plot_step6_wcondition.mat'],...
            'hotspots');
    end
    
  
else
    hotspots = Hcells;
    
end


%*******************************************************************
% HOTSPOTS OVER COARSE MAP
%*******************************************************************
if visualize==1
    
    figure('name','coarse map with hotspots'); 
    imshow(map_env','InitialMagnification',paraPub.ZoomIn); hold on;
    set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])

    %-- coarse map
    fPublishRecon(map_recon,map_env,paraPub);
    pause(1)

    for i = 1:size(hotspots,1)    
        plot(hotspots(i,1),hotspots(i,2),...
             'o','color','k',...
             'MarkerFaceColor','r',...
             'LineWidth',paraPub.HcLineWidth,...
             'MarkerSize',paraPub.HcMarkerSize); %7
    end
    pause(2)
    % filename = 'hotspots-coarse-map';
    % % print('-painters','-dpdf','-r500',filename); % pdf    
    % export_fig([dir_.TomographyMaps,filename], '-pdf','-r500')
    filename = sprintf('%s%s-hotspots-coarse-map',...
                dir_.gdmplan_hotspots,para_.FilePrefix);
    export_fig(filename, '-pdf','-r500');

    pause(1)
    % 
    %pause
end

%-----------------------------------------------------
% hotspots on occupied cells
%-----------------------------------------------------
[ hotspots ] = fMovingOccupiedHcsToFreeSpace( hotspots,...
                                              OBS,...
                                              map_env,...
                                              cellsize_env,...
                                              para_ );

if visualize==1
    pause(1)
    for i = 1:size(hotspots,1)
        plot(hotspots(i,1),hotspots(i,2),...
             'o','color','b','MarkerSize',9);
    end
    pause(2)
    % filename = 'hotspots-coarse-map';
    % % print('-painters','-dpdf','-r500',filename); % pdf    
    % export_fig([dir_.TomographyMaps,filename], '-pdf','-r500')
    %filename = sprintf('%s%s-hotspots-coarse-map',...
    %    dir_.TomographyMaps,para_.FilePrefix);
    %export_fig(filename, '-pdf','-r500');
    
    filename = sprintf('%s%s-hotspots-coarse-map-moved-occupied-hcs',...
                dir_.gdmplan_hotspots,para_.FilePrefix);
    export_fig(filename, '-pdf','-r500');
    
    pause(1)

    %pause
end


% pause


end

function fPublishRecon(map_recon,map_env,paraPub)

% ------------- select color map for gdm --------------
gdm_colormap = flipud(hot(512)); % for main fig. 
%max_concentration = max(map_recon(:)); % main fig.
%max_concentration = 30000 %5000
max_concentration = min( [max(map_recon(:)),paraPub.UpperPPM] ) % main fig.
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
                'MarkerSize',paraPub.CellMarkerSize); %2.75); %4.85
            
            end
        end
    end
end   

end
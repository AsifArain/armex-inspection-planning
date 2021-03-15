function [ sc_true,...
           sc_estimated,...
           true_positives,...
           false_positives,...
           num_of_true_positives,...
           num_of_false_positives,...
           true_positives_dist,...
           false_positives_dist,...
           dist_estimatedSources2trueSources,...
           precision,...
           recall,...
           f_measure ] = fEvaluation4( recon_filename,...
                                       map_env,...
                                       cellsize_env,...
                                       para_,...
                                       dir_,...
                                       visualize )


% fEvaluation performs evaluation of reconstructions.
%
% Author:  Asif Arain
% Project: Simulator for the Evaluation of Sensing Geometries (JFR-2016)


load([dir_.recon,recon_filename],'mean_map','cell_coords_x','cell_coords_y')
map_recon = mean_map';


%**********************************
%   OBSTACLE LIST
%**********************************
OBSenv.ind          = find(~map_env);     % linear indices of obstacles.
[OBSenv_r,OBSenv_c] = find(~map_env);     % rows and columns of obstacles.
OBSenv.sub          = [OBSenv_r OBSenv_c];  % subscripts of obstacles.


%*********************************************************************
%
%                   TRUE GAS SOURCES
%
%*********************************************************************
[ sc_true ] = fTrueGasSources( map_env,...
                               cellsize_env,...
                               OBSenv,...
                               para_,...
                               dir_,...
                               visualize );

%*********************************************************************
%
%                   ESTIMATED GAS SOURCES
%
%*********************************************************************
[ sc_estimated ] = fEstimatedGasSources( map_recon,...
                                         map_env,...
                                         cell_coords_x,...
                                         cell_coords_y,...
                                         cellsize_env,...
                                         sc_true,...
                                         OBSenv,...
                                         para_,...
                                         dir_,...
                                         visualize );
                                     
%*********************************************************************
%
%                   POSITIVE AND FALSE POSITIVES
%
%*********************************************************************
[ true_positives,...
  false_positives,...
  num_of_true_positives,...
  num_of_false_positives,...
  true_positives_dist,...
  false_positives_dist,...
  dist_estimatedSources2trueSources] = fHits( sc_true,...
                                              sc_estimated,...
                                              map_env,...
                                              map_recon,...
                                              cell_coords_x,...
                                              cell_coords_y,...
                                              para_,...
                                              dir_,...
                                              visualize );
                                               
%*********************************************************************
%
%                   F-MEASURE
%
%*********************************************************************
[ precision,...
  recall,...
  f_measure ] = fFMeasure( num_of_true_positives,...
                           num_of_false_positives,...
                           sc_true );


end


function [ precision,...
           recall,...
           f_measure ] = fFMeasure( num_of_true_positives,...
                                    num_of_false_positives,...
                                    sc_true )
%       F-MEASURE
%--------------------------


num_of_true_sources = size(sc_true,1);

precision = num_of_true_positives/(num_of_true_positives+num_of_false_positives);
recall = num_of_true_positives/num_of_true_sources;

f_measure = 2*((precision*recall)/(precision + recall));

end


function [ true_positives,...
           false_positives,...
           num_of_true_positives,...
           num_of_false_positives,...
           true_positives_dist,...
           false_positives_dist,...
           dist_estimatedSources2trueSources] = fHits( sc_true,...
                                                       sc_estimated,...
                                                       map_env,...
                                                       map_recon,...
                                                       cell_coords_x,...
                                                       cell_coords_y,...
                                                       para_,...
                                                       dir_,...
                                                       visualize )
%   POSITIVES
%------------------------------------------

%-- positives association index (estimated source to true source)
%-------------------------------------------------------
% N-D nearest point search
positives_association_ind = dsearchn(sc_true,sc_estimated);

%-----------------
% plot
%-----------------
if visualize==1
    figure('name','Estimated gas sources - hits 1'); 
    imshow(map_env','InitialMagnification',350); hold on;
    set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
    pause(1)
    %-- estimated gas sources centers
    for i = 1:size(sc_estimated,1)
        plot(sc_estimated(i,1),sc_estimated(i,2),...
             'o','color','r','MarkerFaceColor','r','MarkerSize',4.5); %3
    end
    %-- true gas sources centers
    for i = 1:size(sc_true,1)
        plot(sc_true(i,1),sc_true(i,2),...
             'o','color','g','MarkerFaceColor','g','MarkerSize',4); %3
    end
    %-- association
    colors = linspecer(numel(positives_association_ind));
    for i = 1:size(sc_true,1)        
        pts = [sc_estimated(positives_association_ind==i,:);sc_true(i,:)];
        plot(pts(:,1),pts(:,2),...
             'o','color',colors(i,:),'MarkerSize',7.5); %3
    end
    pause(2)
    filename = sprintf('%s%s-hits-association1',...
        dir_.evaluation,para_.FilePrefix);
    %print('-painters','-dpdf','-r500',filename); % pdf
    export_fig(filename, '-pdf','-r500')
    pause(1)
end

% pause
%-----------------------------------------------------------------------
%
% true/false positives association index (estimated sources to true 
% soruces) and distance between them
%
%***********************************************************************
tf_positives_association_ind = positives_association_ind;
dist_estimatedSources2trueSources = zeros(numel(positives_association_ind),1);
for i = 1:numel(tf_positives_association_ind)
            
    pt_from = sc_true(tf_positives_association_ind(i),:);
    pt_to = sc_estimated(i,:);
    %dist_this = pdist2(pt_from,pt_to,'euclidean');
    dist_this = fPathDist2(pt_from,pt_to,map_env);
    
    if dist_this > para_.ScPositiveHits_cells
        tf_positives_association_ind(i) = 0;
    end
    dist_estimatedSources2trueSources(i) = dist_this;
end




%-----------------
% plot
%-----------------
if visualize==1
    figure('name','Estimated gas sources - hits 2'); 
    imshow(map_env','InitialMagnification',350); hold on;
    set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
    pause(1)
    %-- estimated gas sources centers
    for i = 1:size(sc_estimated,1)
        plot(sc_estimated(i,1),sc_estimated(i,2),...
             'o','color','r','MarkerFaceColor','r','MarkerSize',4.5); %3
    end
    %-- true gas sources centers
    for i = 1:size(sc_true,1)
        plot(sc_true(i,1),sc_true(i,2),...
             'o','color','g','MarkerFaceColor','g','MarkerSize',4); %3
    end
    %-- positive association
    colors = linspecer(numel(tf_positives_association_ind));
    for i = 1:size(sc_true,1)        
        pts = [sc_estimated(tf_positives_association_ind==i,:);sc_true(i,:)];
        plot(pts(:,1),pts(:,2),...
             'o','color',colors(i,:),'MarkerSize',7.5); %3
    end
    %-- negative association
    pts = sc_estimated(tf_positives_association_ind==0,:);
    plot(pts(:,1),pts(:,2),...
         'x','color','k','MarkerSize',7.5); %3    
    pause(2)
    filename = sprintf('%s%s-hits-association2',...
        dir_.evaluation,para_.FilePrefix);
    %print('-painters','-dpdf','-r500',filename); % pdf
    export_fig(filename, '-pdf','-r500')
    pause(1)
end



%--------------------------------------------------------------------------
% WHAT IF THERE IS MORE THAN ONE TRUE POSITIVE ASSOCIATED TO A SOURCE?
%                               SOLUTION-1
%
% MULTIPLE TRUE POSITIVES ASSOCIATED TO A TRUE GAS SOURCE ARE MERGED TO
% FIND ATMOST SINGLE TRUE POSITIVE FOR EACH TRUE GAS SOURCE. THIS IS TO
% AVOID OVERALL HIGH NUMBER OF TRUE POSITIVE COUNTS ASSOCITATED TO FEW
% SOURCES.
%
%----------------------------------------------------------------------------
%{


true_positive_association_ind =...
    tf_positives_association_ind(tf_positives_association_ind~=0);
% false_positive_association_ind =...
%     tf_positives_association_ind(tf_positives_association_ind==0);

true_positives = sc_estimated(tf_positives_association_ind~=0,:);
false_positives = sc_estimated(tf_positives_association_ind==0,:);



%----------------------------------------------------------
%
% merg all true positives associated to a true source
%
%**********************************************************

true_weights = scs_weights;

%-- count more than one positives associated to a source
A = true_positive_association_ind;
B = unique(A); 
num_counts = histc(A,B);
high_counts = B(num_counts>1);
%high_counts = find(num_counts>1);

%-- combined new true positives and weights
new_true_positives = zeros(numel(high_counts),2);
new_true_weights = zeros(numel(high_counts),1);
for i = 1:numel(high_counts)
    this_true_positives = true_positives(A==high_counts(i),:);
    this_true_weights = true_weights(A==high_counts(i),:);
    this_wmean_row = 0;
    this_wmean_col = 0;
    this_wmean_wgt = 0;
    for j = 1:size(this_true_positives,1)
        this_wmean_row = this_wmean_row + this_true_weights(j)*this_true_positives(j,1);
        this_wmean_col = this_wmean_col + this_true_weights(j)*this_true_positives(j,2);
        this_wmean_wgt = this_wmean_wgt + this_true_weights(j);
    end
    new_true_positive = [ this_wmean_row/this_wmean_wgt,this_wmean_col/this_wmean_wgt ];
    new_true_positives(i,:) = new_true_positive;
    new_true_weights(i) = this_wmean_wgt;
end

%-- insert new positives and weights into the existing list
for i = 1:numel(high_counts)    
    this_ind = find(A==high_counts(i));
    for j = 1:this_ind        
        true_positives(this_ind(j),:) = new_true_positive(i,:);
        true_weights(this_ind(j)) = new_true_weights(i);
    end
end

%-- now shrink the list (forget weighs, coz they are not useful anymore)
%[C,ia,ic] = unique(true_positives,'rows')
true_positives = unique(true_positives,'rows');
true_positives_association_ind = dsearchn(sc_true,true_positives);


%----------------------------------------------------------
%
% distance from sources to new true/false positives
%
%**********************************************************
tf_positives = [true_positives;false_positives];
tf_positives_association_ind = dsearchn(sc_true,tf_positives);
dist_allPositives2trueSources = zeros(size(tf_positives,1),1);
for i = 1:numel(tf_positives_association_ind)                
    pt_from = sc_true(tf_positives_association_ind(i),:);
    pt_to = tf_positives(i,:);
    %dist_this = pdist2(pt_from,pt_to,'euclidean');
    dist_this = fPathDist2(pt_from,pt_to,map_env);
    if dist_this > para_.ScPositiveHits_cells
        tf_positives_association_ind(i) = 0;
    end
    dist_allPositives2trueSources(i) = dist_this;
end


true_positives_dist = dist_allPositives2trueSources(tf_positives_association_ind~=0,:);
false_positives_dist = dist_allPositives2trueSources(tf_positives_association_ind==0,:);

num_of_true_positives = size(true_positives,1);
num_of_false_positives = size(false_positives,1);

%}




%---------------------------------------------------------------------------
% WHAT IF THERE IS MORE THAN ONE TRUE POSITIVE ASSOCIATED TO A SOURCE?
%                               SOLUTION-2
%
% MULTIPLE TRUE POSITIVES ASSOCIATED TO A TRUE SOURCE ARE TRIED TO ADJUST
% WITH THE OTHER SOURCES IF THEY ARE WITHIN THE DISTANCE AND THE TRUE
% POSTIVES FOR THE OTHER SOURCES IS NOT ALREAD WITH A SMALLER DISTANCE THAN
% THE SUBJECT POSITIVES. ELSE, THEY ARE LISTED AS FALSE POSITIVES.
%
%---------------------------------------------------------------------------

%-- initial list of true positives
true_positives = sc_estimated(tf_positives_association_ind~=0,:);

%-- try to associate redundant true positives to another sources, if not
%possible, then declear them false positives
cond_association = true;
while cond_association
    
    %-- true positive association indices
    this_tp_association_ind = tf_positives_association_ind(tf_positives_association_ind~=0);
    
    %-- find multiple true positives for a source
    A = this_tp_association_ind;    
    B = unique(A); 
    num_counts = histc(A,B);
    high_counts = B(num_counts>1);
    %high_counts = find(num_counts>1);
    
    %-- lets declare them undecided
    for i = 1:numel(high_counts)
        
        ind = find(A==high_counts(i));
        this_true_positives = true_positives(ind,:);
        this_true_source = sc_true(high_counts(i),:);
        this_dist = fPathDist2(this_true_positives,this_true_source,map_env);
        
        [~,min_ind] = min(this_dist);
        undecided_t_ind = ind;
        undecided_t_ind(min_ind) = [];        
        this_tp_association_ind(undecided_t_ind) = -1;
        tf_positives_association_ind(tf_positives_association_ind~=0) = this_tp_association_ind;
    end
    
    %-- lets try to associate them with the other sources or make them
    %false positive
    undecided_ind = find(tf_positives_association_ind==-1);    
    undecided_positives = sc_estimated(tf_positives_association_ind==-1,:);
    undecided_dists = fPathDist2(undecided_positives,sc_true,map_env);    
    for i = 1:size(undecided_positives,1)
        
        this_undecide_dists = undecided_dists(i,:);
        [dist_val,dist_ind] = sort(this_undecide_dists);        
        for j = 1:size(sc_true,1)
            current_dist =...
                dist_estimatedSources2trueSources(tf_positives_association_ind==dist_ind(j));
            if dist_val(j)<para_.ScPositiveHits_cells && dist_val(j)<current_dist
                tf_positives_association_ind(undecided_ind(i)) = dist_ind(j);
                break;
            end
        end
        
        %-- if it is still undecided then it is a false positive
        if tf_positives_association_ind(undecided_ind(i)) == -1
            tf_positives_association_ind(undecided_ind(i)) = 0;
        end
    end
    
    % -- do we need to continue?
    cond_association = any(tf_positives_association_ind(:)==-1);
end


%-- final stuff
true_positives = sc_estimated(tf_positives_association_ind~=0,:);
false_positives = sc_estimated(tf_positives_association_ind==0,:);

true_positives_dist = dist_estimatedSources2trueSources(tf_positives_association_ind~=0,:);
false_positives_dist = dist_estimatedSources2trueSources(tf_positives_association_ind==0,:);

num_of_true_positives = size(true_positives,1);
num_of_false_positives = size(false_positives,1);
%---------------------------------------------------------------------------


%*******************************************************************
%
% true/false positives over reconstruction map
%
%*******************************************************************
if visualize==1 || visualize==0
    
    fPublish_Positives( map_env,...
                        map_recon,...
                        cell_coords_x,...
                        cell_coords_y,...
                        sc_estimated,...
                        sc_true,...
                        tf_positives_association_ind,...
                        para_,...
                        dir_ )
end




end



function [ sc_true ] = fTrueGasSources( map_env,...
                                        cellsize_env,...
                                        OBSenv,...
                                        para_,...
                                        dir_,...
                                        visualize )

true_sources = dlmread([dir_.ROSLogsThis,'gas_sources.dat']);

%-------------------------------
%  agglomerative clustering
%-------------------------------
cut_off = 2;

% - cluster linkage tree
Z = linkage(true_sources,'single');
% figure; dendrogram(Z);

% - high concentration clusters
Tclusters = cluster(Z,'cutoff',cut_off,'criterion','distance');
% - number of clusters
num_of_clusters = numel(unique(Tclusters));

% -- gas sources center is the mean position
sc_true = zeros(num_of_clusters,2);
for i = 1:num_of_clusters
    Hcluster = [true_sources(Tclusters==i,1),true_sources(Tclusters==i,2)];
    sc_true(i,1) = sum(Hcluster(:,1))/size(Hcluster,1);
    sc_true(i,2) = sum(Hcluster(:,2))/size(Hcluster,1);
end



%************************************************
% moving sources to unoccupied space
%************************************************
%-- moving hcs from occupied cells to the nearest unoccupied cells
[ sc_true ] = fMovingOccupiedHcsToFreeSpace( sc_true,...
                                             OBSenv,...
                                             map_env,...
                                             cellsize_env,...
                                             para_ );



%-----------------
% plot
%-----------------
if visualize==1
    figure('name','True gas sources centers'); 
    imshow(map_env','InitialMagnification',350); hold on;
    set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
    pause(1)    
    %-- clusters
    colors = linspecer(num_of_clusters);
    for i = 1:num_of_clusters
        plot(true_sources(Tclusters==i,1)+0.5,...
             true_sources(Tclusters==i,2)+0.5,...
             's',...
             'color',colors(i,:),...
             'MarkerFaceColor',colors(i,:),...
             'MarkerSize',5); %2.75
    end
    % -- gas sources centers
    for i = 1:size(sc_true,1)
        plot(sc_true(i,1),sc_true(i,2),...
             'o','color','k','MarkerFaceColor','r','MarkerSize',5.5); %3
    end
    pause(2)
    filename = sprintf('%s%s-true-sources',...
        dir_.evaluation,para_.FilePrefix);
    %print('-painters','-dpdf','-r500',filename); % pdf
    export_fig(filename, '-pdf','-r500')
    pause(1)
end


end


function [ sc_estimated ] = fEstimatedGasSources( map_recon,...
                                                  map_env,...
                                                  cell_coords_x,...
                                                  cell_coords_y,...
                                                  cellsize_env,...
                                                  sc_true,...
                                                  OBSenv,...
                                                  para_,...
                                                  dir_,...
                                                  visualize )

%------------------------------------------
% aligned reconstruction map
%------------------------------------------
map_recon2 = zeros(size(map_env));
for i = 1:numel(cell_coords_x)
    for j = 1:numel(cell_coords_y)
        map_recon2(cell_coords_x(i),cell_coords_y(j)) = map_recon(i,j);
    end
end


%-----------------------------
% plot (reconstruction map)
%------------------------------
if visualize==1
    
    figure('name','Estimated gas sources on the recon map - final'); 
    imshow(map_env','InitialMagnification',350); hold on;
    set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])

    %------------------------------------------------------------
    cell_marker_size = 4.85; %5; %4.85;
    %color_scale = 'linear'; %
    %--------------- select color map for gdm --------------
    gdm_colormap = flipud(autumn(512)); % for main fig. 

    %------ color step size --------------------------------
    max_concentration = min( [max(map_recon(:)),para_.PublishUpperPPM] ); % main fig.
    delta = max_concentration/(size(gdm_colormap,1)-1);
    
    xyzPoints = zeros(numel(cell_coords_x)*numel(cell_coords_y),3);
    colorPoints = zeros(numel(cell_coords_x)*numel(cell_coords_y),3);
    
    reconColorR = zeros(numel(cell_coords_x)*numel(cell_coords_y),1);
    reconColorG = zeros(numel(cell_coords_x)*numel(cell_coords_y),1);
    reconColorB = zeros(numel(cell_coords_x)*numel(cell_coords_y),1);

    k = 1;
    for i = 1:numel(cell_coords_x)
        for j = 1:numel(cell_coords_y)
            % if positive concentration
            if map_recon(i,j)>=1 %para_.concentrationThresholdEvaluation%50%0
                %if map_env(i,j)>0
                %------------------- useful trash -------------------
                %col = gdm_colormap(round(map_gd(i,j)/delta)+1,:);

                % --- if concentration is greater than 1000 ppm, its still 1000 ppm
                %col = gdm_colormap( round(min(map_gd(i,j),3000)/delta)+1, :);
                %col = gdm_colormap( round(min(map_gd(i,j),1.78e+04)/delta)+1, :);
                %col = gdm_colormap( round(min(map_gd(i,j),inf)/delta)+1, :);
                %-----------------------------------------------------

                %----------------------------------------------------
                % linear color code
                %----------------------------------------------------                
                linear_color_number = round( min(map_recon(i,j),max_concentration)/delta)+1;
                col = gdm_colormap( linear_color_number,:);                
                
                %----- to avoid black borders of earch cell -----
                % {
                if sum(col) == 3
                    col = col-(1e-10);
                end
                %}
                %----- plot ------------
                xyzPoints(k,:) = [cell_coords_x(i),cell_coords_y(j),0];                
                colorPoints(k,:) = col;                
                reconColorR(k) = round(col(1)*255);
                reconColorG(k) = round(col(2)*255);
                reconColorB(k) = round(col(3)*255);

                plot(cell_coords_x(i)+0.5,cell_coords_y(j)+0.5,...
                    's',...
                    'MarkerEdgeColor',col,...
                    'MarkerFaceColor',col,...
                    'MarkerSize',cell_marker_size);
                k = k+1;
                %end
            end
        end
    end
    pause(1)
    %-- save
    filename = sprintf('%s%s-recon-map',dir_.evaluation,para_.FilePrefix);
    %print('-painters','-dpdf','-r500',filename); % pdf
    export_fig(filename, '-pdf','-r500');
    pause(1)
end



%------------------------------------------
% positive concentration map
%------------------------------------------
map_positive = zeros(size(map_env));
map_positive((map_recon2(:)>para_.concentrationThresholdEvaluation)) = 1;

%------------------------------------------
% - high concentration cells
%------------------------------------------
[high_r,high_c] = find(map_positive);
Hcells = [high_r,high_c];


%------------------------------------------
% plot
%------------------------------------------
if visualize==1
    figure('name','Estimated gas sources - positive concentrations'); 
    imshow(map_env','InitialMagnification',350); hold on;
    set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
    color_positive_con = [255,140,0]/255;
    pause(1)
    for i = 1:size(Hcells,1)
        plot(Hcells(i,1)+0.5,Hcells(i,2)+0.5,...
             's',...
             'color',color_positive_con,...
             'MarkerFaceColor',color_positive_con,...
             'MarkerSize',4.85); %2.75
    end
    pause(2)
    %filename = sprintf('%s%s-estimated-sources-positive-concentrations',...
    filename = sprintf('%s%s-estimated-sources-positive-concentrations',...
        dir_.evaluation,para_.FilePrefix);
    export_fig(filename, '-pdf','-r500')
    pause(1)
end


%******************************************************
% 
%               CLUSTERING
% 
%******************************************************

%-- linkage cutoff for clusters
% cut_off = 2 %2.5 %2;% 1; %3.5; %median([Z(end-2,3) Z(end-1,3)])


if size(Hcells,1)>1

    
    
    
    %*******************************************
    % Step-1: CLUSTERING
    % clusters based on distance
    %*******************************************

    % - cluster linkage tree
    Z = linkage(Hcells,'single');
    % figure; dendrogram(Z);
    
    % - high concentration clusters
    Hclusters1 = cluster(Z,'cutoff',para_.ScCutOffClusters,'criterion','distance');
    
    % - number of clusters
    num_of_aglo_clusters = numel(unique(Hclusters1));
    
    
    %--- plot here: clustering step 1---------------
    
    
    

    %*******************************************
    % Step-2: CLUSTERING
    % sub-divide lengthy clusters
    %*******************************************
    
    hot_clusters = zeros(size(Hcells,1),1);
    clust_num = 1;
    
    for i = 1:num_of_aglo_clusters

        cluster1_this = [Hcells(Hclusters1==i,1),Hcells(Hclusters1==i,2)];

        dist_this = zeros(size(cluster1_this,1));
        for j = 1:size(dist_this,1)
            dist_this(j,:) = pdist2(cluster1_this(j,:),cluster1_this,'euclidean');
        end
        
        % {
        %num_of_desired_clusters2_this = round(max(dist_this(:))/(SensorRange.cell/2));
        num_of_desired_clusters2_this =...
            round(max(dist_this(:))/(para_.SCsClusterDist_cells));
        if num_of_desired_clusters2_this<1
            num_of_desired_clusters2_this = 1;
        end
        %}
        %-- forced number
        %num_of_desired_clusters2_this = 1

        if num_of_desired_clusters2_this>1
            
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
    num_of_kmean_clusters = numel(unique(hot_clusters));
    
    
    %--------------------
    % plot (Step-2)
    %--------------------
    if visualize==1
        figure('name','Estimated gas sources - kmean clustering');
        imshow(map_env','InitialMagnification',350); hold on;
        set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
        pause(1)
        %-- clusters
        %colors = linspecer(num_of_clusters+num_of_clusters1);
        colors = lines(num_of_kmean_clusters);
        for i = 1:num_of_kmean_clusters
            plot(Hcells(hot_clusters==i,1)+0.5,...
                 Hcells(hot_clusters==i,2)+0.5,...
                 's',...
                 'color',colors(i,:),... colors(num_of_clusters1+i,:),...
                 'MarkerFaceColor',colors(i,:),... colors(num_of_clusters1+i,:),...
                 'MarkerSize',4.85); %2.75
             %-- cluster num
             x_text = Hcells(hot_clusters==i,1);
             y_text = Hcells(hot_clusters==i,2) ;            
             text(x_text(end),y_text(end),num2str(i),'FontSize',6);
        end
        pause(2)
        %-- save
        filename = sprintf('%s%s-estimated-sources-clustering-step2',...
            dir_.evaluation,para_.FilePrefix);
        %print('-painters','-dpdf','-r500',filename); % pdf
        export_fig(filename, '-pdf','-r500');
        pause(1)
    end

    
    
    
    %-----------------
    % plot (Step-1)
    %-----------------
    if visualize==1
        figure('name','Estimated gas sources - agglomerative clustering');
        imshow(map_env','InitialMagnification',350); hold on;
        set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
        pause(1)
        %-- clusters
        %colors = linspecer(num_of_clusters1);
        colors = lines(num_of_kmean_clusters);
        for i = 1:num_of_aglo_clusters    
            plot(Hcells(Hclusters1==i,1)+0.5,...
                 Hcells(Hclusters1==i,2)+0.5,...
                 's',...
                 'color',colors(i,:),...
                 'MarkerFaceColor',colors(i,:),...
                 'MarkerSize',4.85); %2.75
             %- cluster num
             x_text = Hcells(Hclusters1==i,1);
             y_text = Hcells(Hclusters1==i,2);
             text(x_text(end),y_text(end),num2str(i),'FontSize',6)
        end
        pause(2)
        %-- save
        filename = sprintf('%s%s-estimated-sources-clustering-step1',...
            dir_.evaluation,para_.FilePrefix);
        %print('-painters','-dpdf','-r500',filename); % pdf
        export_fig(filename, '-pdf','-r500')
        pause(1)
    end
    
    
    
    
    
    
    
    %************************************************
    % Step-3: ESTIMATE SOURCES
    % gas sources are weighted mean in each cluster
    %************************************************
    
    % -- estimated sources are weighted mean    
    %Wmap = map_recon2./max(map_recon2(:)); % max weight is 1
    Wmap = map_recon2./sum(map_recon2(:)); % total weight is 1
    sc_estimated = zeros(num_of_kmean_clusters,2);
    scs_weights = zeros(num_of_kmean_clusters,1);
    scs_ppmm = zeros(num_of_kmean_clusters,1);
    for i = 1:num_of_kmean_clusters

        Hcluster = [Hcells(hot_clusters==i,1),Hcells(hot_clusters==i,2)];

        Hcluster_wmean_row = 0;
        Hcluster_wmean_col = 0;
        Hcluster_wmean_wgt = 0;
        Hcluster_ppmm      = 0;
        for j = 1:size(Hcluster,1)
            Hcluster_wmean_row = Hcluster_wmean_row + ...
                Wmap(Hcluster(j,1),Hcluster(j,2))*Hcluster(j,1);
            Hcluster_wmean_col = Hcluster_wmean_col + ...
                Wmap(Hcluster(j,1),Hcluster(j,2))*Hcluster(j,2);
            Hcluster_wmean_wgt = Hcluster_wmean_wgt + ...
                Wmap(Hcluster(j,1),Hcluster(j,2));
            Hcluster_ppmm = Hcluster_ppmm +...
                map_recon2(Hcluster(j,1),Hcluster(j,2));
        end
        
        Hcluster_wmean = [ Hcluster_wmean_row/Hcluster_wmean_wgt,...
                           Hcluster_wmean_col/Hcluster_wmean_wgt ];
        sc_estimated(i,:) = round(Hcluster_wmean);
        scs_weights(i) = Hcluster_wmean_wgt;
        scs_ppmm(i) = Hcluster_ppmm;
        
    end
    
    
    %----------
    % plot
    %----------
    if visualize==1
        figure('name','Estimated gas sources - clustering step3');
        imshow(map_env','InitialMagnification',350); hold on;
        set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
        pause(1)
        %-- clusters
        %colors = linspecer(num_of_kmean_clusters+num_of_aglo_clusters);
        colors = lines(num_of_kmean_clusters);
        for i = 1:num_of_kmean_clusters
            plot(Hcells(hot_clusters==i,1)+0.5,...
                 Hcells(hot_clusters==i,2)+0.5,...
                 's',...
                 'color',colors(i,:),...colors(num_of_aglo_clusters+i,:),...
                 'MarkerFaceColor',colors(i,:),...colors(num_of_aglo_clusters+i,:),...
                 'MarkerSize',4.85); %2.75
        end
        %-- estimated sources
        sc_estimated_plot = sortrows(sc_estimated,1);
        for i = 1:size(sc_estimated,1)
            plot(sc_estimated_plot(i,1),sc_estimated_plot(i,2),...
                 'o','color','k','MarkerFaceColor','r','MarkerSize',4); %7 %5
             %- estimated source num             
             text(sc_estimated_plot(i,1)+1.0,sc_estimated_plot(i,2),num2str(i),'FontSize',6)
        end
        pause(2)
        %-- save
        filename = sprintf('%s%s-estimated-sources-clustering-step3',...
            dir_.evaluation,para_.FilePrefix);
        %print('-painters','-dpdf','-r500',filename); % pdf
        export_fig(filename, '-pdf','-r500');
        pause(1)
    end
    
    %pause
    
    
    %************************************************
    % Step-4: ESTIMATE SOURCES
    % move estimated sources to unoccupied space
    %************************************************
    %-- moving hcs from occupied cells to the nearest unoccupied cells
    [ sc_estimated ] = fMovingOccupiedHcsToFreeSpace( sc_estimated,...
                                                      OBSenv,...
                                                      map_env,...
                                                      cellsize_env,...
                                                      para_ );
    
    %----------
    % plot
    %----------
    if visualize==1
        figure('name','Estimated gas sources - clustering step4');
        imshow(map_env','InitialMagnification',350); hold on;
        set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
        pause(1)
        %-- clusters
        %colors = linspecer(num_of_kmean_clusters+num_of_aglo_clusters);
        colors = lines(num_of_kmean_clusters);
        for i = 1:num_of_kmean_clusters
            plot(Hcells(hot_clusters==i,1)+0.5,...
                 Hcells(hot_clusters==i,2)+0.5,...
                 's',...
                 'color',colors(i,:),... colors(num_of_aglo_clusters+i,:),...
                 'MarkerFaceColor',colors(i,:),... colors(num_of_aglo_clusters+i,:),...
                 'MarkerSize',4.85); %2.75
        end
        %-- estimated sources
        sc_estimated_plot = sortrows(sc_estimated,1);
        for i = 1:size(sc_estimated,1)            
            plot(sc_estimated_plot(i,1),sc_estimated_plot(i,2),...
                 'o',...
                 'color','k',...
                 'MarkerFaceColor','r',...
                 'MarkerSize',4); %7 %5
             %- estimated source num             
             text(sc_estimated_plot(i,1)+1.0,sc_estimated_plot(i,2),num2str(i),'FontSize',6)
        end
        pause(2)
        %-- save
        filename = sprintf('%s%s-estimated-sources-clustering-step4',...
            dir_.evaluation,para_.FilePrefix);
        %print('-painters','-dpdf','-r500',filename); % pdf
        export_fig(filename, '-pdf','-r500');
        pause(1)
    end
    
    
    %************************************************
    % Step-5: ESTIMATE SOURCES
    % nearby sources are combined
    %************************************************    
    
    
    %-- merg all hcs close enough
    dist_ind_kaamkay = 1;
    testind = 0;
    while numel(dist_ind_kaamkay)>=1

        testind = testind+1;
        %hcs_dist = fPathDist2(hotspots,hotspots,map_env);
        
        %-- combinations and pairwise distance
        %-----------------------------------------------
        scs_comb_ind = nchoosek(1:size(sc_estimated,1),2);
        scs_comb_dist = inf(size(scs_comb_ind,1),1);
        for i = 1:size(scs_comb_dist,1)
            %-- from -- to
            sc_from = sc_estimated(scs_comb_ind(i,1),:);
            sc_to   = sc_estimated(scs_comb_ind(i,2),:);
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
            
            sc_first = sc_estimated(scs_comb_ind(dist_ind_kaamkay(1),1),:);
            w_first  = scs_weights(scs_comb_ind(dist_ind_kaamkay(1),1));
            ppm_first = scs_ppmm(scs_comb_ind(dist_ind_kaamkay(1),1));
            
            sc_second = sc_estimated(scs_comb_ind(dist_ind_kaamkay(1),2),:);
            w_second  = scs_weights(scs_comb_ind(dist_ind_kaamkay(1),2));
            ppm_second = scs_ppmm(scs_comb_ind(dist_ind_kaamkay(1),2));
            
            w_mirged = w_first + w_second;
            
            sc_mirged_r = ((w_first*sc_first(1,1)) + (w_second*sc_second(1,1))) / w_mirged;
            sc_mirged_c = ((w_first*sc_first(1,2)) + (w_second*sc_second(1,2))) / w_mirged;
            
            sc_mirged = [sc_mirged_r,sc_mirged_c];
            
            ppm_mirged = ppm_first + ppm_second;
            
            sc_estimated(scs_comb_ind(dist_ind_kaamkay(1),1),:) = sc_mirged;
            sc_estimated(scs_comb_ind(dist_ind_kaamkay(1),2),:) = [];
            
            scs_weights(scs_comb_ind(dist_ind_kaamkay(1),1)) = w_mirged;
            scs_weights(scs_comb_ind(dist_ind_kaamkay(1),2)) = [];
            
            scs_ppmm(scs_comb_ind(dist_ind_kaamkay(1),1)) = ppm_mirged;
            scs_ppmm(scs_comb_ind(dist_ind_kaamkay(1),2)) = [];
            
        end
    end
    
    
    %----------
    % plot
    %----------
    if visualize==1
        figure('name','Estimated gas sources - step5');
        imshow(map_env','InitialMagnification',350); hold on;
        set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
        pause(1)
        %-- clusters
        %colors = linspecer(num_of_kmean_clusters+num_of_aglo_clusters);
        colors = lines(num_of_kmean_clusters);
        for i = 1:num_of_kmean_clusters
            plot(Hcells(hot_clusters==i,1)+0.5,...
                 Hcells(hot_clusters==i,2)+0.5,...
                 's',...
                 'color',colors(i,:),... colors(num_of_aglo_clusters+i,:),...
                 'MarkerFaceColor',colors(i,:),...colors(num_of_aglo_clusters+i,:),...
                 'MarkerSize',4.85); %2.75
        end
        %-- estimated sources
        sc_estimated_plot = sortrows(sc_estimated,1);
        for i = 1:size(sc_estimated,1)
            plot(sc_estimated_plot(i,1),sc_estimated_plot(i,2),...
                 'o',...
                 'color','k',...
                 'MarkerFaceColor','r',...
                 'MarkerSize',4); %7 %5
             %- estimated source num             
             text(sc_estimated_plot(i,1)+1.0,sc_estimated_plot(i,2),num2str(i),'FontSize',6)
        end
        pause(2)
        %-- save
        filename = sprintf('%s%s-estimated-sources-clustering-step5',...
            dir_.evaluation,para_.FilePrefix);
        %print('-painters','-dpdf','-r500',filename); % pdf
        export_fig(filename, '-pdf','-r500');
        pause(1)
    end
    
    %pause
    %----------------
    %   plot
    %-----------------
    %{
    if visualize==1
        figure('name','clusters-2 with estimated sources'); 
        imshow(map_env','InitialMagnification',350); hold on;
        set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
        pause(1)
        colors = lines(num_of_clusters);
        for i = 1:num_of_clusters
            plot(Hcells(hot_clusters==i,1)+0.5,Hcells(hot_clusters==i,2)+0.5,...
                 's',...
                 'color',colors(i,:),...
                 'MarkerFaceColor',colors(i,:),...
                 'MarkerSize',5); %5 %4
        end        
        % -- estimated sources
        for i = 1:size(sc_estimated,1)
            plot(sc_estimated(i,1),sc_estimated(i,2),...
                 'o','color','k','MarkerFaceColor','r','MarkerSize',7); %7 %5
        end
    end
    %}
    
    
    
    
    
    %************************************************
    % Step-6: ESTIMATE SOURCES
    % shortlist sources based on concentration weights
    %************************************************
    hcs_prevCluster  = sc_estimated;
    ws_prevCluster   = scs_weights;
    ppmm_prevCluster = scs_ppmm;
    %sc_estimated     = round(hcs_prevCluster);
    sc_estimated     = round(hcs_prevCluster*100)/100;
    
    
    %*************************************************
    %-- high pass filtering of hcs based on weights
    %*************************************************
    total_weight = sum(ws_prevCluster);
    scs_weights = ws_prevCluster;
    scs_ppmm = ppmm_prevCluster;
    num_of_scs = size(sc_estimated,1);
    
    %qualification_threshould = (total_weight/num_of_hcs*0.5)
    qualification_w = (total_weight/num_of_scs*para_.SCsCutoffWeightFactor);
        
    sc_estimated = sc_estimated( scs_weights>qualification_w |...
                                 scs_ppmm>para_.SCsCutoffPPMM,: );
    
    scs_weights = scs_weights( scs_weights>qualification_w |...
                               scs_ppmm>para_.SCsCutoffPPMM,: );
                             
                             
    %----------
    % plot
    %----------
    if visualize==1
        figure('name','Estimated gas sources - clustering step6');
        imshow(map_env','InitialMagnification',350); hold on;
        set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
        pause(1)
        %-- clusters
        %colors = linspecer(num_of_kmean_clusters+num_of_aglo_clusters);
        colors = lines(num_of_kmean_clusters);
        for i = 1:num_of_kmean_clusters
            plot(Hcells(hot_clusters==i,1)+0.5,...
                 Hcells(hot_clusters==i,2)+0.5,...
                 's',...
                 'color',colors(i,:),...colors(num_of_aglo_clusters+i,:),...
                 'MarkerFaceColor',colors(i,:),...colors(num_of_aglo_clusters+i,:),...
                 'MarkerSize',4.85); %2.75
        end
        %-- estimated sources
        for i = 1:size(sc_estimated,1)
            plot(sc_estimated(i,1),sc_estimated(i,2),...
                 'o',...
                 'color','k',...
                 'MarkerFaceColor','r',...
                 'MarkerSize',4); %7 %5
        end
        pause(2)
        %-- save
        filename = sprintf('%s%s-estimated-sources-clustering-step6',...
            dir_.evaluation,para_.FilePrefix);
        %print('-painters','-dpdf','-r500',filename); % pdf
        export_fig(filename, '-pdf','-r500');
        pause(1)
    end
    
    
else
    
    sc_estimated = Hcells;    
end



%************************************************
% ESTIMATED SOURCES OVER CLUSTERS
%************************************************
%{
if visualize==1
    figure('name','clusters with hotspots'); 
    imshow(map_env','InitialMagnification',350); hold on;
    set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
    pause(1)
    % color_clusters = lines(num_of_clusters);
    for i = 1:num_of_clusters
        plot(Hcells(hot_clusters==i,1)+0.5,Hcells(hot_clusters==i,2)+0.5,...
             's',...
             'color',colors(i,:),... colors(num_of_clusters1+i,:)
             'MarkerFaceColor',colors(i,:),... colors(num_of_clusters1+i,:)
             'MarkerSize',2.75); %5
    end
    % -- estimated sources
    for i = 1:size(sc_estimated,1)
        plot(sc_estimated(i,1),sc_estimated(i,2),...
             'o','color','k','MarkerFaceColor','r','MarkerSize',3); %7
    end

    pause(2)
    filename = sprintf('%s%s-hotspots-clusters',...
        dir_.Solutions,para_.ExperimentTitle);
    %print('-painters','-dpdf','-r500',filename); % pdf
    export_fig(filename, '-pdf','-r500');
    pause(1)
end
%}




%*******************************************************************
% ESTIMATED SOURCES OVER RECONSTRUCTION MAP
%*******************************************************************
if visualize==1
    
    figure('name','Estimated gas sources on the recon map - final'); 
    imshow(map_env','InitialMagnification',350); hold on;
    set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])

    %------------------------------------------------------------
    cell_marker_size = 4.85; %5; %4.85;
    color_scale = 'linear'; %
    %--------------- select color map for gdm --------------
    gdm_colormap = flipud(autumn(512)); % for main fig. 

    %------ color step size --------------------------------
    max_concentration = min( [max(map_recon(:)),para_.PublishUpperPPM] ); % main fig.
    delta = max_concentration/(size(gdm_colormap,1)-1);

    %------------------- nonlinear color code ---------------------  
    if strcmp(color_scale,'nonlinear')
        %b = (exp(256/1))/255;
        %a = 256 /(log(b*256));

        x = 1:1:256;
        epsilon = 1/(exp(0.10)-1); % slower
        normlogsum_penalty = log(1+(abs(x)/epsilon));
        normlogsum_penalty = normlogsum_penalty/(normlogsum_penalty(x==max(x))); % normalized
        normlogsum_penalty = round(256*normlogsum_penalty);
        normlogsum_penalty(1)
        normlogsum_penalty(end)
        %figure; plot(x,normlogsum_penalty); % to test

        % publish nonlinear steps
        linear_percentage = 10:10:100;
        linear_percentage_ind = round((linear_percentage/100)*numel(x));
        nonlinear_percentage_ind = normlogsum_penalty(linear_percentage_ind);
        nonlinear_percentage = (nonlinear_percentage_ind/numel(x))*100;
        nonlinear_percentage = round(nonlinear_percentage);
        %disp('nonlinear_percentage = ',num2str(nonlinear_percentage))
        disp_stng = ['nonlinear color code scale is: ',num2str(nonlinear_percentage)];
        %disp(disp_stng)
        disp(disp_stng)
    end
    %--------------------------------------------------------------

    xyzPoints = zeros(numel(cell_coords_x)*numel(cell_coords_y),3);
    colorPoints = zeros(numel(cell_coords_x)*numel(cell_coords_y),3);


    reconColorR = zeros(numel(cell_coords_x)*numel(cell_coords_y),1);
    reconColorG = zeros(numel(cell_coords_x)*numel(cell_coords_y),1);
    reconColorB = zeros(numel(cell_coords_x)*numel(cell_coords_y),1);

    k = 1;
    for i = 1:numel(cell_coords_x)
        for j = 1:numel(cell_coords_y)
            % if positive concentration
            if map_recon(i,j)>=50 %para_.concentrationThresholdEvaluation%50%0
                %if map_env(i,j)>0
                %------------------- useful trash -------------------
                %col = gdm_colormap(round(map_gd(i,j)/delta)+1,:);

                % --- if concentration is greater than 1000 ppm, its still 1000 ppm
                %col = gdm_colormap( round(min(map_gd(i,j),3000)/delta)+1, :);
                %col = gdm_colormap( round(min(map_gd(i,j),1.78e+04)/delta)+1, :);
                %col = gdm_colormap( round(min(map_gd(i,j),inf)/delta)+1, :);
                %-----------------------------------------------------

                %----------------------------------------------------
                % linear color code
                %----------------------------------------------------
                if strcmp(color_scale,'linear')
                % {
                linear_color_number =...
                    round( min(map_recon(i,j),max_concentration)/delta)+1;
                col = gdm_colormap( linear_color_number,:);                
                %}
                end
                %----------------------------------------------------

                %----------------------------------------------------
                % nonlinear color code
                %----------------------------------------------------
                if strcmp(color_scale,'nonlinear')
                % {
                linear_color_number =...
                    round( min(map_recon(i,j),max_concentration)/delta)+1;
                %nonlinear_color_number = round( a*log(b*linear_color_number) );
                %nonlinear_color_number = linear_color_number;
                nonlinear_color_number = normlogsum_penalty(linear_color_number);
                %if nonlinear_color_number <1
                %    nonlinear_color_number = 1;
                %end
                col = gdm_colormap( nonlinear_color_number,:);                
                %}
                end
                %----------------------------------------------------

                %----- to avoid black borders of earch cell -----
                % {
                if sum(col) == 3
                    col = col-(1e-10);
                end
                %}
                %----- plot ------------
                xyzPoints(k,:) = [cell_coords_x(i),cell_coords_y(j),0];                
                colorPoints(k,:) = col;                
                reconColorR(k) = round(col(1)*255);
                reconColorG(k) = round(col(2)*255);
                reconColorB(k) = round(col(3)*255);

                plot(cell_coords_x(i)+0.5,cell_coords_y(j)+0.5,...
                    's',...
                    'MarkerEdgeColor',col,...
                    'MarkerFaceColor',col,...
                    'MarkerSize',cell_marker_size);
                k = k+1;
                %end
            end
        end
    end
    pause(1)
    %--------------------------
    %--- true sources
    %--------------------------
    for i = 1:size(sc_true,1)    
        plot(sc_true(i,1),sc_true(i,2),...
             'o',...
             'color','k',...
             'MarkerFaceColor','b',...
             'MarkerSize',4.5); %7
    end
    %--------------------------
    %--- estimated sources
    %--------------------------
    for i = 1:size(sc_estimated,1)    
        plot(sc_estimated(i,1),sc_estimated(i,2),...
             'o',...
             'color','k',...
             'MarkerFaceColor','r',...
             'MarkerSize',4.5); %7
    end
    %--------------------------
    %-- LEGENDS
    %--------------------------     
    % estimated sources
    plot(160,48,'o','color',[0.1,0.1,0.1],'MarkerFaceColor','r','MarkerSize',4.5);
    % true sources
    plot(160,43,'o','color',[0.1,0.1,0.1],'MarkerFaceColor','b','MarkerSize',4.5);
    pause(2)
    %-- save
    filename = sprintf('%s%s-estimated-sources-recon-map-final',...
        dir_.evaluation,para_.FilePrefix);
    %print('-painters','-dpdf','-r500',filename); % pdf
    export_fig(filename, '-pdf','-r500');
    pause(1)
end

end




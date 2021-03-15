function [ hotspots ] = fHotspotsICRA16( map_env,...
                                         map_recon,...
                                         SensorRange,...
                                         OBS,...
                                         para_,...
                                         cellsize_env,...
                                         visualize )
%
%



%------------------------------------
%    PUBLISH PARAMETERS (PAPER)
%------------------------------------
paraPub.ZoomIn = 750;
paraPub.CellMarkerSize = 10;%7.0;
paraPub.BoundryLineWidth = 1;
paraPub.BoundryLineColor = [0.3,0.3,0.3];
paraPub.HcMarkerSize = 10;
paraPub.HcLineWidth = 2;
paraPub.UpperPPM = 20000; %30000
%--------------------------------

% %------------------------------------
% %    PUBLISH PARAMETERS (THESIS)
% %------------------------------------
% paraPub.ZoomIn = 350;
% paraPub.CellMarkerSize = 2.95 %10;%7.0;
% paraPub.BoundryLineWidth = 1;
% paraPub.BoundryLineColor = [0.3,0.3,0.3];
% paraPub.HcMarkerSize = 4 %10;
% paraPub.HcLineWidth = 1;
% paraPub.UpperPPM = 20000; %30000

%--------------------------------

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
% GAS DISTRIBUTION MAP - HOTSPOTS (HIGH CONCENTRATION):
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++

% threshold to declear a cell of high concentration
% highConcentrationThreshold_PPM = 100;

map_hotspot = zeros(size(map_recon));
% ind = find(gd_map_hotspot(:))>500;
% gd_map_hotspot(ind) = 1;
map_hotspot(find(map_recon(:)>para_.highConcentrationThreshold_PPM)) = 1;
% test_map_resized2_obs = find(test_map_resized2(:)==0);
% gd_map_coverage(test_map_resized2_obs) = 0;

% plot
figure('name',['high concentration map of ',...
    num2str(para_.highConcentrationThreshold_PPM),'ppm',100]); 
%imshow(map_env); hold on;
imshow(map_env','InitialMagnification',paraPub.ZoomIn); hold on; %350
set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
%color_positive_con = [255,140,0]/255;
color_positive_con = [218,165,32]/255;
pause(1)

for i = 1:size(map_hotspot,1)
    for j = 1:size(map_hotspot,2)        
        if map_hotspot(i,j)
            %plot(j,i,'sb','MarkerFaceColor','b');
            plot(i,j,...
             's',...
             'color',color_positive_con,...
             'MarkerFaceColor',color_positive_con,...
             'MarkerSize',paraPub.CellMarkerSize);
        end
    end
end
pause(2)
%filename = 'hotspots-positive-concentrations';
% print('-painters','-dpdf','-r500',filename); % pdf    
%export_fig([dir_.TomographyMaps,filename], '-pdf','-r500')
filename = sprintf('icra16-hotspots-high-concentrations');
export_fig(filename, '-pdf','-r500');
pause(1)
% pause

%pause
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
% hotspots on occupied cells
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++

occupiedHotspots = OBS.ind(find(ismember(OBS.ind,find(map_hotspot))));
[occupiedHotspots_r,occupiedHotspots_c] = ind2sub(size(map_env),occupiedHotspots);

for m = 1:numel(occupiedHotspots)
    
    % map to look for free cells
    map_FreeCell = zeros(size(map_env));
    map_FreeCell(occupiedHotspots(m)) = 1;
        
    % max radius to look for unoccupied cell
    %radiusFreeCell_m = 2 % thickness of candidate conf circle (meters)
    radiusFreeCell_cell = para_.radiusFreeCell_m/para_.MapRes;

    % grow area 
    se = strel('disk',radiusFreeCell_cell,0); % outer limit
    map_FreeCell = imdilate(map_FreeCell,se);
    
    % remove occupied cells from the region
    map_FreeCell(OBS.ind(find(ismember(OBS.ind,find(map_FreeCell))))) = 0;
    
    % remove other high concentration cells
    highConcentrationCells = find(map_hotspot);
    map_FreeCell(highConcentrationCells(find(ismember(highConcentrationCells,find(map_FreeCell))))) = 0;
    
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
    
    % plot
    figure('name','nearest hotspot is added'); imshow(map_env); hold on;
    for i = 1:size(map_hotspot,1)
        for j = 1:size(map_hotspot,2)        
            if map_hotspot(i,j)
                plot(j,i,'or','MarkerFaceColor','r');
            end
        end
    end
    
    % prevois centroid hotspot is removed
    map_hotspot(occupiedHotspots(m)) = 0;
    
    % plot
    figure('name','hotspot on occupied cell is removed'); imshow(map_env); hold on;
    for i = 1:size(map_hotspot,1)
        for j = 1:size(map_hotspot,2)        
            if map_hotspot(i,j)
                plot(j,i,'or','MarkerFaceColor','r');
            end
        end
    end
    
end


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
% GAS DISTRIBUTION MAP - CENTROID HOTSPOTS:
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++

% distance between hotspots to look for a centroid
% centroid_dist_m = 1.5; 
% centroid_dist_cell = para_.centroid_dist_m/cell_length_m;
centroid_dist_cell = para_.centroid_dist_m/cellsize_env;


% idices of high concentration cells
hotspot_ind = find(map_hotspot);
[hotspot_r,hotspot_c] = find(map_hotspot);

% distance matrix
hotspot_distance = zeros(numel(hotspot_ind));
for i = 1:numel(hotspot_ind)
    for j = 1:numel(hotspot_ind)
        hotspot_distance(i,j) = sqrt((hotspot_r(j)-hotspot_r(i))^2+...
                                     (hotspot_c(j)-hotspot_c(i))^2);
    end
end
% distance from a hotspot to itself is inf.
for i = 1:numel(hotspot_ind)
    hotspot_distance(i,i) = inf;
end

% finding nearest hotspots to combine them
nearest_spots = zeros(numel(hotspot_ind));
for i = 1:numel(hotspot_ind)
    nearest_spots(i,i) = 1;    
    for j = 1:numel(hotspot_ind)
        if hotspot_distance(i,j) <= centroid_dist_cell %5 % within 5 cells
            nearest_spots(i,j) = 1;
        end
    end
end

% remove redundant spots
nearest_spots = unique(nearest_spots,'rows');
% combining nearnest hotspots
A = nearest_spots;
while sum(sum(A,1),2)~=size(A,2)
    B = [];
    for i = 1:size(A,2)
        B = [B;sum(A(find(A(:,i)),:),1)];        
    end
    C = unique(B,'rows');
    D = zeros(size(C));
    D(find(C)) = 1;
    A = D;
end
combined_spots = A;
hotspot_rc = [hotspot_r,hotspot_c];
centroid_spots_rc = zeros(size(combined_spots,1),2);
for i = 1:size(combined_spots,1)
    if sum(combined_spots(i,:),2) > 1
        pts = hotspot_rc(find(combined_spots(i,:)),:);         
        centroid_spots_rc(i,:) = round(centroid(pts));
    else
        centroid_spots_rc(i,:) = hotspot_rc(find(combined_spots(i,:)),:);
    end
end
centroid_spots_ind = sub2ind(size(map_hotspot),...
    centroid_spots_rc(:,1),centroid_spots_rc(:,2));

% now create centroid hotspot map
map_centroid_hotspot = zeros(size(map_hotspot));
map_centroid_hotspot(centroid_spots_ind) = 1;
% figure; imshow(map_centroid_hotspot);

% plot
figure('name','centroid hotspots'); imshow(map_env); hold on;
for i = 1:size(map_centroid_hotspot,1)
    for j = 1:size(map_centroid_hotspot,2)        
        if map_centroid_hotspot(i,j)
            plot(j,i,'or','MarkerFaceColor','r');
        end
    end
end




% plot
figure('name','centroid hotspots'); 
%imshow(map_env); hold on;
imshow(map_env','InitialMagnification',paraPub.ZoomIn); hold on; %350
set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
%color_positive_con = [255,140,0]/255;
color_positive_con = [218,165,32]/255;
pause(1)
for i = 1:size(map_hotspot,1)
    for j = 1:size(map_hotspot,2)        
        if map_hotspot(i,j)
            %plot(j,i,'sb','MarkerFaceColor','b');
            plot(i,j,...
             's',...
             'color',color_positive_con,...
             'MarkerFaceColor',color_positive_con,...
             'MarkerSize',paraPub.CellMarkerSize);
        end
    end
end
pause(2)

for i = 1:size(map_centroid_hotspot,1)
    for j = 1:size(map_centroid_hotspot,2)        
        if map_centroid_hotspot(i,j)
            plot(i,j,'or','MarkerFaceColor','r','MarkerSize',paraPub.CellMarkerSize);
        end
    end
end

%filename = 'hotspots-positive-concentrations';
% print('-painters','-dpdf','-r500',filename); % pdf    
%export_fig([dir_.TomographyMaps,filename], '-pdf','-r500')
filename = sprintf('icra16-hotspots-centroids');
export_fig(filename, '-pdf','-r500');
pause(1)
pause



% -----------------------------------------------------------------------
% centroid hotspots on occupied cells
% -----------------------------------------------------------------------

occupiedCHotspots = OBS.ind(find(ismember(OBS.ind,find(map_centroid_hotspot))));
[occupiedCHotspots_r,occupiedCHotspots_c] = ind2sub(size(map_env),occupiedCHotspots);

for m = 1:numel(occupiedCHotspots)
    
    map_FreeCell = zeros(size(map_env));
    map_FreeCell(occupiedCHotspots(m)) = 1;
    
    % max radius to look for unoccupied cell
    %radiusFreeCell_m = 2 % thickness of candidate conf circle (meters)
    %radiusFreeCell_cell = radiusFreeCell_m/paramPreprocess.MapRes;

    % grow area 
    se = strel('disk',radiusFreeCell_cell,0); % outer limit
    map_FreeCell = imdilate(map_FreeCell,se);
    
    % remove occupied cells from the region
    map_FreeCell(OBS.ind(find(ismember(OBS.ind,find(map_FreeCell))))) = 0;
    
    % indices/subscripts of unoccupied (candidate) cells around
    map_FreeCell_ind = find(map_FreeCell);
    [map_FreeCell_r,map_FreeCell_c] = ind2sub(size(map_env),map_FreeCell_ind);
        
    % distance vector between occupied hotspot and unoccupied cells around
    dist_CHotspot2FreeCells = zeros(size(map_FreeCell_ind));
    
    for n = 1:numel(map_FreeCell_ind)
        dist_CHotspot2FreeCells(n) = sqrt(((map_FreeCell_c(n)-occupiedCHotspots_c(m))^2)+...
                                          ((map_FreeCell_r(n)-occupiedCHotspots_r(m))^2));
    end
    
    % nearest cell
    [~,cell_ind] = min(dist_CHotspot2FreeCells);
    
    % nearnest cell is new centroid hotspot.
    map_centroid_hotspot(map_FreeCell_ind(cell_ind)) = 1;
    
    % plot
    figure('name','nearest centroid hotspot is added'); imshow(map_env); hold on;
    for i = 1:size(map_centroid_hotspot,1)
        for j = 1:size(map_centroid_hotspot,2)        
            if map_centroid_hotspot(i,j)
                plot(j,i,'or','MarkerFaceColor','r');
            end
        end
    end
    
    % prevois centroid hotspot is removed
    map_centroid_hotspot(occupiedCHotspots(m)) = 0;
    
    % plot
    figure('name','centroid hotspot on occupied cell is removed'); imshow(map_env); hold on;
    for i = 1:size(map_centroid_hotspot,1)
        for j = 1:size(map_centroid_hotspot,2)        
            if map_centroid_hotspot(i,j)
                plot(j,i,'or','MarkerFaceColor','r');
            end
        end
    end
    
end


% hotspots_indices(ismember(hotspots_indices,occupiedHotspots)) =...
%     occupiedHotspots_relocated_ind;
[hotspots_r,hotspots_c] = ind2sub(size(map_env),find(map_centroid_hotspot));

hotspots = [hotspots_r,hotspots_c];





%///////////////////////////////////////////// 
% HOTSPOTS OVER COARSE MAP
%/////////////////////////////////////////////
% visualize = 1
if visualize==1
    
figure('name','coarse map with hotspots'); 
imshow(map_env','InitialMagnification',750); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])

% ------------- select color map for gdm --------------
gdm_colormap = flipud(hot(512)); % for main fig. 
max_concentration = max(map_recon(:)); % main fig.
delta = max_concentration/(size(gdm_colormap,1)-1);

[gd_row,gd_col] = size(map_recon);

for i = 1:gd_row
    for j = 1:gd_col
        % if positive concentration

        if map_recon(i,j)>=0

            linear_color_number = round( min(map_recon(i,j),max_concentration)/delta)+1;
            col = gdm_colormap( linear_color_number,:);                
            %----- to avoid black borders of earch cell -----
            % {
            if sum(col) == 3
                col = col-(1e-10);
            end
            plot(i,j,...
                's',...
                'MarkerEdgeColor',col,...
                'MarkerFaceColor',col,...
                'MarkerSize',4.85);
        end
    end
end   

pause(1)

for i = 1:size(hotspots,1)    
    plot(hotspots(i,1),hotspots(i,2),...
         'o','color','k','MarkerFaceColor','r','MarkerSize',7);
end

pause(2)
print('-painters','-dpdf','-r500','hotspots-over-coarse-map-icra16-old'); % pdf
pause(1)
end
 


pause


end

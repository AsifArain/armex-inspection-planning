function [ hc_s ] = fMovingOccupiedHcsToFreeSpace( hc_s,...
                                                   OBS,...
                                                   map_env,...
                                                   cellsize_m,...
                                                   para_ )
% moving Hcs on the occupied cells to the nearest unoccupied cells
%

%-- round off
%-----------------------------
hc_s = round(hc_s);

%-- indices of hc cells
%-----------------------------
hc_i = sub2ind(size(map_env),hc_s(:,1),hc_s(:,2));
%-- occupied cells
%------------------------------
%occupiedHotspots = OBS.ind(ismember(OBS.ind,hotspots_indices)) % buggy
occupiedHc_i = hc_i(ismember(hc_i,OBS.ind)); % correct
[occupiedHc_r,occupiedHc_c] = ind2sub(size(map_env),occupiedHc_i);
occupiedHc_s = [occupiedHc_r,occupiedHc_c];

%-- occupied cells - relocated
%------------------------------
occupiedHc_relocated_i = zeros(numel(occupiedHc_i),1);

for i = 1:numel(occupiedHc_i)
    
    %--- map to look for free cells
    %-------------------------------------------
    map_FreeCell = zeros(size(map_env));
    map_FreeCell(occupiedHc_i(i)) = 1;
        
    %-- max radius to look for unoccupied cell
    %--------------------------------------------    
    %radiusFreeCell_cell = para_.radiusFreeCell_m/dataYAML.resolution
    radiusFreeCell_cell = para_.radiusFreeCell_m/cellsize_m;
    
    %-- grow area
    %-------------------------------------------
    se = strel('disk',radiusFreeCell_cell,0); % outer limit
    map_FreeCell = imdilate(map_FreeCell,se);
    
    % remove occupied cells from the region
    %-------------------------------------------
    map_FreeCell(OBS.ind(ismember(OBS.ind,find(map_FreeCell)))) = 0;
    
    %{
    figure('name','This free cells map'); 
    imshow(map_FreeCell','InitialMagnification',350); hold on;
    set(gca,'YDir','normal'); colormap('gray'); caxis([-1 1.25])
    pause(1)    
    plot(occupiedHc_r(i)+0.0,occupiedHc_c(i)+0.0,...
        'xb','MarkerSize',7);
    %}
    
    % indices/subscripts of unoccupied cells around
    %-------------------------------------------------------
    freeCells_i = find(map_FreeCell);
    [freeCells_r,freeCells_c] = find(map_FreeCell); %ind2sub(size(map_env),map_FreeCell_ind);
    freeCells_s = [freeCells_r,freeCells_c];
        
    % distance vector between occupied hotspot and unoccupied cells around
    %-----------------------------------------------------------
    %dist_this = fPathDist2(sc_from,sc_to,map_env);
    dist_hc2FreeCells = pdist2(occupiedHc_s(i,:),freeCells_s,'euclidean');
    
    %-- nearest cell
    %-----------------------------
    [~,cell_ind] = min(dist_hc2FreeCells);
    
    % nearnest cell is new centroid hotspot.
    %-------------------------------------------------
    %map_hotspot(freeCells_i(cell_ind)) = 1;    
    occupiedHc_relocated_i(i) = freeCells_i(cell_ind);
end


hc_i(ismember(hc_i,occupiedHc_i)) = occupiedHc_relocated_i;
[hc_r,hc_c] = ind2sub(size(map_env),hc_i);

hc_s = [hc_r,hc_c];

end
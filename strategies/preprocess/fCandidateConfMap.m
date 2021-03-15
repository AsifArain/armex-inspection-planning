function [ map_candConf ] = fCandidateConfMap( hc_sub,...
                                               map_env,...
                                               SensorRange,...
                                               OBSconf,...
                                               visualize )
% 
% 


[ confMap_RangeBased ] = fConfMapRangeBased( hc_sub,...
                                             map_env,...
                                             SensorRange,...
                                             OBSconf,...
                                             visualize );
                                         
confMapR_ind = find(confMap_RangeBased(:));
[confMapR_r,confMapR_c] = ind2sub(size(map_env),confMapR_ind);

validInd = zeros(numel(confMapR_ind),1);

for i = 1:numel(confMapR_ind)
    
    startCell_sub = hc_sub;
    endCell_sub = [confMapR_r(i),confMapR_c(i)];

    % validatation for the visibility
    [ validation ] = fValidation( startCell_sub,endCell_sub,OBSconf );

    % update 
    validInd(i) = validation;
    
end

confMap_RangeBased(confMapR_ind(validInd==0)) = 0;
map_candConf = confMap_RangeBased;

                                              
end


function [ confMap_RangeBased ] = fConfMapRangeBased( hc_sub,...
                                                      map_env,...
                                                      SensorRange,...
                                                      OBSconf,...
                                                      visualize )
% 
% 
% [hotspot_r,hotspot_c] = ind2sub(size(map_env),hc_this);
hc_r = hc_sub(1);
hc_c = hc_sub(2);


% initialize map for candidate conf. -- ring only
confMap_RangeBased = zeros(size(map_env));    
confMap_RangeBased(hc_r,hc_c) = 1;

% -- sample candidate conf within circle radius of sensing range
% candCircleThickness_m = 0 % thickness of candidate conf circle (meters)
% candCircleThickness_cell = candCircleThickness_m/paramPreprocess.MapRes
radius_candCircleOut = ceil(SensorRange.cell);


% grow hotspot to max radius
se = strel('disk',radius_candCircleOut,0); % outer limit
confMap_RangeBased = imdilate(confMap_RangeBased,se);

% -- plot
if visualize == 1
figure('name','cand conf - max limit and obs are removed');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
map_candConf_ind = find(confMap_RangeBased(:)==1);
[map_candConf_r,map_candConf_c] = ind2sub(size(map_env),map_candConf_ind);
plot(map_candConf_r,map_candConf_c,'xb','MarkerFaceColor','b');        
plot(hc_r,hc_c,'xr','MarkerFaceColor','r');
end

% -- remove obs cells
confMap_RangeBased(OBSconf.ind(ismember(OBSconf.ind,find(confMap_RangeBased)))) = 0;

% plot
if visualize == 1
figure('name','cand conf - max limit and obs are removed');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
map_candConf_ind = find(confMap_RangeBased(:)==1);
[map_candConf_r,map_candConf_c] = ind2sub(size(map_env),map_candConf_ind);
plot(map_candConf_r,map_candConf_c,'xb','MarkerFaceColor','b');        
plot(hc_r,hc_c,'xr','MarkerFaceColor','r');
end




end



function [ validation ] = fValidation( startCell_sub,endCell_sub,OBSconf )
% 
% 

% angle to the center of the cell
ang = atan2(((endCell_sub(2))-(startCell_sub(2))),...
    ((endCell_sub(1))-(startCell_sub(1))))*(180/pi);
ang = ang+(360*(ang<0));

% conditions to project the line along (+/-x,+/-y).
cond_x1 = ang>=000&&ang<=045;
cond_x2 = ang>=135&&ang<=225;
cond_x3 = ang>=315&&ang<=360;

cond_y1 = ang>=045&&ang<=135;
cond_y2 = ang>=225&&ang<=315;

% sampling length.
delta = 0.01;
% Note: fine/coarse sampling changes the results for intersection points
% and hence changes the list of cells interseted. TODO: Find the optimal
% value for delta.

if cond_x1 || cond_x3
    
    % intermediate x and y points.
    x = startCell_sub(1):delta:endCell_sub(1);
    y = x*tand(ang);
    
    % 4 decimal precision
    x = round(x*1e+4)/1e+4;
    y = round(y*1e+4)/1e+4;
    
    % corresponding cell numbers.
    cell_xy = round([x',y']);
    
    % cells that belong to the intersection point (x.5,y.5)
    ind_x  = find(abs(fix(x)-x)==0.5);
    ind_yy = find(abs(fix(y(ind_x))-y(ind_x))==0.5);
    ind_y  = ind_x(ind_yy);
    
    for i = 1:numel(ind_y)
        cell_xy = [cell_xy;...
                   ceil(x(ind_y(i))),  ceil(y(ind_y(i)))  ;...
                   ceil(x(ind_y(i))),  floor(y(ind_y(i))) ;...
                   floor(x(ind_y(i))), ceil(y(ind_y(i)))  ;...
                   floor(x(ind_y(i))), floor(y(ind_y(i)))];
    end    
    
elseif cond_x2
    
    % intermediate x and y points.
    x = startCell_sub(1):-delta:endCell_sub(1);
    y = x*tand(ang);
    
    % 4 decimal precision
    x = round(x*1e+4)/1e+4;
    y = round(y*1e+4)/1e+4;
    
    % corresponding cell numbers.
    cell_xy = round([x',y']);
    
    % cells that belong to the intersection point (x.5,y.5)
    ind_x  = find(abs(fix(x)-x)==0.5);
    ind_yy = find(abs(fix(y(ind_x))-y(ind_x))==0.5);
    ind_y  = ind_x(ind_yy);
    
    for i = 1:numel(ind_y)
        cell_xy = [cell_xy;...
                   ceil(x(ind_y(i))),  ceil(y(ind_y(i)))  ;...
                   ceil(x(ind_y(i))),  floor(y(ind_y(i))) ;...
                   floor(x(ind_y(i))), ceil(y(ind_y(i)))  ;...
                   floor(x(ind_y(i))), floor(y(ind_y(i)))];
    end
    
    
elseif cond_y1
    
    % intermediate x and y points.
    y = startCell_sub(2):delta:endCell_sub(2);
    x = y/tand(ang);
    
    % 4 decimal precision
    x = round(x*1e+4)/1e+4;
    y = round(y*1e+4)/1e+4;
    
    % corresponding cell numbers.
    cell_xy = round([x',y']);
    
    % cells that belong to the intersection point (x.5,y.5)
    ind_x  = find(abs(fix(x)-x)==0.5);
    ind_yy = find(abs(fix(y(ind_x))-y(ind_x))==0.5);
    ind_y  = ind_x(ind_yy);
    
    for i = 1:numel(ind_y)
        cell_xy = [cell_xy;...
                   ceil(x(ind_y(i))),  ceil(y(ind_y(i)))  ;...
                   ceil(x(ind_y(i))),  floor(y(ind_y(i))) ;...
                   floor(x(ind_y(i))), ceil(y(ind_y(i)))  ;...
                   floor(x(ind_y(i))), floor(y(ind_y(i)))];
    end
    
    
elseif cond_y2
    
    % intermediate x and y points.
    y = startCell_sub(2):-delta:endCell_sub(2);
    x = y/tand(ang);
    
    % 4 decimal precision
    x = round(x*1e+4)/1e+4;
    y = round(y*1e+4)/1e+4;
    
    % corresponding cell numbers.
    cell_xy = round([x',y']);
    
    % cells that belong to the intersection point (x.5,y.5)
    ind_x  = find(abs(fix(x)-x)==0.5);
    ind_yy = find(abs(fix(y(ind_x))-y(ind_x))==0.5);
    ind_y  = ind_x(ind_yy);
    
    for i = 1:numel(ind_y)
        cell_xy = [cell_xy;...
                   ceil(x(ind_y(i))),  ceil(y(ind_y(i)))  ;...
                   ceil(x(ind_y(i))),  floor(y(ind_y(i))) ;...
                   floor(x(ind_y(i))), ceil(y(ind_y(i)))  ;...
                   floor(x(ind_y(i))), floor(y(ind_y(i)))];
    end
    
end

% remove repeated cells
% cell_xy = unique(cell_xy,'rows');

validation = 1;
for i = 1:numel(OBSconf.ind)
    obs = [OBSconf.sub(i,1),OBSconf.sub(i,2)];
    
    % validation for the visibility
    ind = find(cell_xy(:,1)==obs(1)&cell_xy(:,2)==obs(2));
    if size(ind,1) > 0
        validation = 0;
        break;
    end
end


end


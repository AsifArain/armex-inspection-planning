function [ tDomain ] = fDomain( FoV,FoVfaces,SensorRange )
% 

tDomaini = tic; % initialize pre computational time.

f = FoVfaces.num; % for simplification.

%**********************************************************************
% COMPUTE DOMAIN
%**********************************************************************
fComputeDomain( f,...
                FoV,...
                FoVfaces,...
                SensorRange );

%**********************************************************************
% Total computation time.
%**********************************************************************
tDomain = 1e-4*(round(toc(tDomaini)*1e4));

fprintf(1,'\nComputation time: %0.4f sec \n',tDomain);




end



function [ ] = fComputeDomain( f,...
                               FoV,...
                               FoVfaces,...
                               SensorRange )
%fV generates Visibility matrix V.

% --- ANGLES:
FoVfaces.lgt = [1/FoVfaces.num:1/FoVfaces.num:1-(1/FoVfaces.num) 0]'*360;
FoVfaces.alf = FoVfaces.lgt-(FoV/2);
FoVfaces.bay = FoVfaces.lgt+(FoV/2);
% compensate negative angles
%FoVfaces.alf(FoVfaces.alf<0) = 360+FoVfaces.alf(FoVfaces.alf<0);
FoVfaces.alf = wrapTo360(FoVfaces.alf);
FoVfaces.bay = wrapTo360(FoVfaces.bay);
FoVfaces.ang = [FoVfaces.alf FoVfaces.bay];
FoVfaces.ang = [FoVfaces.alf FoVfaces.bay];


%***********************************************************
%                     Sensor Domain
%***********************************************************
disp('.. Sensor Domain')

        
        % For each cell of the map.
        sD_all = cell(1,f); % sensor domain (visible).
        oD_all = cell(1,f); % sensor's osbtacle domain (used to check the 
                        % visibility of "sD" for each cell in "oD"]

        % a single reference cell in the map.
        for i = 1
            
            % sensor position.
            sensor = [0,0];
            
            % For each sensing configuration.
            for j = 1:f
                                
                % zawia := initial, middle and final (opening) angle.
                zawia = [FoVfaces.ang(j,1) FoVfaces.lgt(j) FoVfaces.ang(j,2)];

                % Finding sensor domain, i.e. visible cells for the
                % configuration. 
                [ visibleRange ] = fSenDomain( sensor,SensorRange,zawia );
                % visibleRange := visible range (row,col).

                % post-processing for sorting and removing duplicates
                visibleRange = unique(visibleRange,'rows'); % remove repeated numbers
                visibleRange = sortrows(visibleRange);      % sort rows
                
                % here we have jth sensor's domain.
                sD_all{j} = visibleRange;
                
                
                % since, to varify if a cell in sensor's domain is visible,
                % we need to check all the cells that effect the
                % visibility. Some of them are even not in the sensing
                % range, therefore, we need to compute obstacle domain that
                % contains the cells in the sensor domain and connecting.
                % "obsDomain" will be used to check visibility of the
                % sensor's domain.

                % obstacle domain
                [ obsRange ] = fSenDomainOBS( sensor,SensorRange,zawia );
                % obsRange := obstacle range (row,col).
                
                % post-processing for sorting and removing duplicates
                obsRange = unique(obsRange,'rows'); % remove repeated numbers
                obsRange = sortrows(obsRange);      % sort rows
                
                % here we have jth sensor's domain.
                oD_all{j} = obsRange;
                
            end
        end

        %--------------------------------------
        % EFFECTED VISIBILITY BY EACH OBSTACLE
        %--------------------------------------
        
        disp('.. Visibility vs Obstacles')
        %VObs = cell(1,f); % Visibility vs Obstacle
        
                
        
        for i = 1
            sensor = [0,0];
            % For each configuration of a cell.
            for j = 1:f
                
                j
                
                % let take sensor's mask
                sensorDom  = sD_all{j};
                % and obstacle domain
                obstacDom  = oD_all{j};
                
                %VObs{j} = zeros(size(obstacDom,1),size(sensorDom,1)); 
                VObs = zeros(size(obstacDom,1),size(sensorDom,1)); 
                % rows    := obstacle cell,
                % columns := visibility effect.
                
                %for k = 1:size(sensorDom,1)
                for k = 1:size(obstacDom,1)
                    
                    % kth cell is an obstacle
                    %obstacDomOBS(k,3)  = 0; 
                    obs = obstacDom(k,:);
                    
                    % visibility effect on all the cells
                    [ vDom ] = fVisSDvObs( sensor,sensorDom,obs,SensorRange );
                    
                    % update visibility against kth obstacle
                    %VObs{j}(k,:) = vDom;
                    VObs(k,:) = vDom;
                   
                end
                
                
                %filename = sprintf('sD_oD_VObs_Tomography_range%02d_fov%03d_fovfaces%03d.mat',...
                %                    SensorRange.cell,FoV,f);
                %save(filename,'sD','oD','VObs','-V7.3');
                
                
                this_file = sprintf('sD_oD_VObs_rng%02d_fov%03d_orn%03d.mat',...
                                    SensorRange.cell,FoV,j);
                                
                this_dir = sprintf('../sD_oD_VObs/sensing_range_%02dcell/fov_%03ddeg/orientation_%03ddeg/',SensorRange.cell,FoV,j);
    
                if ~exist(this_dir,'dir')
                    mkdir(this_dir);
                end

                sD = sD_all{j};
                oD = oD_all{j};
                VObs;


                save([this_dir,this_file],'sD','oD','VObs','-V7.3');

    
               
            end
        end
        
        
        %filename = sprintf('sD_oD_VObs_Tomography_range%02d_fov%03d_fovfaces%03d.mat',...
        %                    SensorRange.cell,FoV,f);
        %save(filename,'sD','oD','VObs','-V7.3');

end




function [ visibleRange ] = fSenDomain( sensor,SensorRange,zawia )
%fSenDomain returns visible cells (domain/mask) for a given sensing
%configuration. 
%

% INTERPOLATING INTERMEDIATE ANGLES:
% -------------------------------------
%{
DELTA_ANGLE = 1;
% --- Intermediate angles of visibile window. ---
% First half:
angl_1h = zawia(1):DELTA_ANGLE:zawia(2);
if zawia(1) > zawia(2) % if final angle is greater than 360 degrees.
    angl_1h = zawia(1):DELTA_ANGLE:360;
end
% Second half:
angl_2h = zawia(2):DELTA_ANGLE:zawia(3);
if zawia(2) > zawia(3) % if final angle is greater than 360 degrees.
    angl_2h = [zawia(2):DELTA_ANGLE:360 0:DELTA_ANGLE:zawia(3)];
end
angl = unique([angl_1h angl_2h]);
%}

DELTA_ANGLE = 1;
if     zawia(2) > zawia(1)
    angl_1h = zawia(1):DELTA_ANGLE:zawia(2);
elseif zawia(2) < zawia(1) % if final angle is greater than 360 degrees.
    angl_1h = [zawia(1):DELTA_ANGLE:360 0:DELTA_ANGLE:zawia(2)];
end

% Second half:
if     zawia(3) > zawia(2)
    angl_2h = zawia(2):DELTA_ANGLE:zawia(3);    
elseif zawia(3) < zawia(2) % if final angle is greater than 360 degrees.
    angl_2h = [zawia(2):DELTA_ANGLE:360 0:DELTA_ANGLE:zawia(3)];
end

angl = unique([angl_1h angl_2h]);
angl = wrapTo360(angl);

% angl = zawia
visibleRange = [];        %
% subCellANGspecial = cell(1,4); % subject cells of special angles.

for i = 1:numel(angl)
    ALF = angl(i);
    % Find visible cells along one edge point.
    [ visCells ] = fRayTracSD(sensor,SensorRange,ALF);
    
    % Updating visible range.
    visibleRange = unique([visibleRange; visCells],'rows');   

end

end

function [ effective_visCells ] = fRayTracSD(sensor,SensorRange,ALF)
% fRayTracSD is ray tracing for sensor domain.

% CONSTANTS:
stepSize = 0.1; % step size of the ray length while sampling from origion to max. 

% initialize some high number of indices
effective_visCells = zeros(round(2*length(0.5:stepSize:SensorRange.cell)),2);
num = 0;
for d = [0.5:stepSize:SensorRange.cell SensorRange.cell]
    dr = 0.001*round((d*cosd(ALF))/0.001);  % row
    dc = 0.001*round((d*sind(ALF))/0.001);  % column
    
    vis_r = round(sensor(1)+round(dr)); % visible cell if validated
    vis_c = round(sensor(2)+round(dc)); % visible cell if validated
    
    if abs((dc-0.5)-fix(dc-0.5))>0.5 && abs((dr-0.5)-fix(dr-0.5))>0.5
        num = num+1;
        effective_visCells(num,:) = [vis_r vis_c];
    end
end
effective_visCells(num+1:end,:) = [];
effective_visCells = unique(effective_visCells,'rows');


end

function [ obsRange ] = fSenDomainOBS( sensor,SensorRange,zawia )
%fSenDomainOBS returns sensor's domain to check visibility against given
%obstacles. It contains more cells than the sensor's domain (mask).
%

% INTERPOLATING INTERMEDIATE ANGLES:
% ----------------------------------
%{
DELTA_ANGLE = 1;
% --- Intermediate angles of visibile window. ---
% First half:
angl_1h = zawia(1):DELTA_ANGLE:zawia(2);
if zawia(1) > zawia(2) % if final angle is greater than 360 degrees.
    angl_1h = zawia(1):DELTA_ANGLE:360;
end
% Second half:
angl_2h = zawia(2):DELTA_ANGLE:zawia(3);
if zawia(2) > zawia(3) % if final angle is greater than 360 degrees.
    angl_2h = [zawia(2):DELTA_ANGLE:360 0:DELTA_ANGLE:zawia(3)];
end
angl = unique([angl_1h angl_2h]);
%}

DELTA_ANGLE = 1;
if     zawia(2) > zawia(1)
    angl_1h = zawia(1):DELTA_ANGLE:zawia(2);
elseif zawia(2) < zawia(1) % if final angle is greater than 360 degrees.
    angl_1h = [zawia(1):DELTA_ANGLE:360 0:DELTA_ANGLE:zawia(2)];
end

% Second half:
if     zawia(3) > zawia(2)
    angl_2h = zawia(2):DELTA_ANGLE:zawia(3);    
elseif zawia(3) < zawia(2) % if final angle is greater than 360 degrees.
    angl_2h = [zawia(2):DELTA_ANGLE:360 0:DELTA_ANGLE:zawia(3)];
end

angl = unique([angl_1h angl_2h]);
angl = wrapTo360(angl);

% angl = zawia;

% initialize the obstacle range.
obsRange = []; 

% find the list of cells (for obstacles) in the sensing range.
for i = 1:numel(angl)
    
    % angle.
    ALF = angl(i);
    
    % Finding cells that are intersected along the ray.
    [ obsCells ] = fRayTracOD(sensor,SensorRange,ALF);
    
    % Updating the list of cells to be checked for obstacles.
    %obsRange = unique([obsRange; obsCells],'rows');
    obsRange = [obsRange; obsCells];

end


end

function [ obsCells ] = fRayTracOD(sensor,SensorRange,ALF)
% fRayTracOD is ray tracing to find sensor domain for visibility check
% against listed obstacles

% point corresponding to the sensor range and ray angle.
subPoint_x = sensor(1) + (SensorRange.cell*cosd(ALF));
subPoint_y = sensor(2) + (SensorRange.cell*sind(ALF));

subPoint = [subPoint_x,subPoint_y];


% conditions.
cond_x1 = ALF>=0&&ALF<=45;
cond_x2 = ALF>=135&&ALF<=225;
cond_x3 = ALF>=315&&ALF<=360;

cond_y1 = ALF>=45&&ALF<=135;
cond_y2 = ALF>=225&&ALF<=315;

% sampling length.
delta = 0.01;
% Note: fine/coarse sampling changes the results for intersection points
% and hence changes the list of cells interseted. TODO: Find the optimal
% value for delta.


if cond_x1 || cond_x3
    
    % intermediate x and y points.
    x = [sensor(1):delta:subPoint(1) subPoint(1)];
    y = x*tand(ALF);
    
    % 4 decimal precision
    x = round(x*1e+4)/1e+4;
    y = round(y*1e+4)/1e+4;
    
    % corresponding cell numbers.
    cell_xy = round([x',y']);
    
    % cells that belong to the intersection point (x.5,y.5)
    ind_x  = find(abs(fix(x)-x)==0.5);
    ind_yy = abs(fix(y(ind_x))-y(ind_x))==0.5;
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
    x = [sensor(1):-delta:subPoint(1) subPoint(1)];
    y = x*tand(ALF);
    
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
    y = [sensor(2):delta:subPoint(2) subPoint(2)];
    x = y/tand(ALF);
    
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
    y = [sensor(2):-delta:subPoint(2) subPoint(2)];
    x = y/tand(ALF);
    
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


 obsCells = unique(cell_xy,'rows');

end

function [ vDom_new ] = fVisSDvObs( sensor,sensorDom,obs,SensorRange )
%fVisibileArea returns visible cells for a given sensing configuration.
%


% initialize the list of visible cells.
% vDom = zeros(size(sensorDom,1),1);
vDom_new = ones(size(sensorDom,1),1);
% vDom_old = ones(size(sensorDom,1),1);


%-- efficient method
%=============================================
%{
for i = 1:size(sensorDom,1)
    %i
    % cell to be validated for visibility
    subCell = [sensorDom(i,1) sensorDom(i,2)];
    %subCell = [eSensorDom(i,1) eSensorDom(i,2)];

    % validatation for the visibility
    %[ cellValid ] = fCellValidSD( sensor,subCell,obstacDomOBS );
    [ validation ] = fCellValidSD( sensor,subCell,obs );

    % update visDom
    vDom_old(i) = validation;
        
end
%}


%-- efficient method
%=============================================

dist_sensor2obs = distancePoints(sensor,obs);

%-- apply the efficient technique if obs is reasonable far from the sensor
if dist_sensor2obs >= 2
    
    %-- efficient sensor domain
    eSensorDom = fEffectiveSensorDomain( sensor,obs,SensorRange );
    
    %-- consider efficient domain present in sensor domain only
    eSensorDom_validation = zeros(size(eSensorDom,1),1);
    for i = 1:size(eSensorDom,1)
        
        ind = find(eSensorDom(i,1)==sensorDom(:,1) & eSensorDom(i,2)==sensorDom(:,2));
        eSensorDom_validation(i) = size(ind,1);
    end
    eSensorDom = eSensorDom(find(eSensorDom_validation),:);
    
    
else
    eSensorDom = sensorDom;
    
end

%-- troubleshooting
%==============================================
%{
figure; hold on;
xlim([-SensorRange.cell-4,SensorRange.cell+4])
ylim([-SensorRange.cell-4,SensorRange.cell+4])
axis equal

plot(sensorDom(:,1),sensorDom(:,2),'sb','MarkerSize',8);

if size(eSensorDom,1)>0
plot(eSensorDom(:,1),eSensorDom(:,2),'sr','MarkerSize',5.5);
end

plot(0,0,'ok')
plot(obs(1),obs(2),'xk')

axis equal
%}


% for all cells in the sensor domain, check visibility against obs.
for i = 1:size(eSensorDom,1)
    %i
    % cell to be validated for visibility
    subCell = [eSensorDom(i,1) eSensorDom(i,2)];

    % validatation for the visibility
    [ validation ] = fCellValidSD( sensor,subCell,obs );

    % update visDom
    ind = find(eSensorDom(i,1)==sensorDom(:,1) & eSensorDom(i,2)==sensorDom(:,2));
    vDom_new(ind) = validation;
    
end


%{
test = find(vDom_new(:)~=vDom_old(:));
if size(test,1) >0
    
    [vDom_new'; vDom_old']
    
    disp('issue occured')
    
    plot(sensorDom(test,1),sensorDom(test,2),'dg','MarkerSize',4);    
    pause
end

% close all;
% clc;
%}


end


function [ validation ] = fCellValidSD( sensor,subCell,obs )
% 
% 

% angle to the center of the cell
ang = atan2(((subCell(2))-(sensor(2))),((subCell(1))-(sensor(1))))*(180/pi);
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
    x = sensor(1):delta:subCell(1);
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
    x = sensor(1):-delta:subCell(1);
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
    y = sensor(2):delta:subCell(2);
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
    y = sensor(2):-delta:subCell(2);
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

% validation for the visibility
ind = find(cell_xy(:,1)==obs(1)&cell_xy(:,2)==obs(2));
if size(ind,1) > 0
    validation = 0;
else
    validation = 1;
end
% validation = (size(ind,1) < 1);

end

%--------------------------- End of the document ------------------------------------

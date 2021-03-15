function [ cellsScanned ] = fSingleConfCoverage( sensing_conf,...
                                                 map_env,...
                                                 FoV,...
                                                 SensorRange )
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  V
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -- orientation
sensing_conf(3) = sensing_conf(3) + (sensing_conf(3)==0)*360;
sD = []; oD = []; VObs = [];
this_orn = sensing_conf(3);

% this_dir = sprintf(['sD_oD_VObs/sensing_range_%02dcell/fov_%03ddeg/',...
%                     'orientation_%03ddeg/'],SensorRange.cell,FoV,this_orn);
this_dir = sprintf('visibility/templates/sensing_range_%02dcell/fov_%03ddeg/',...
        SensorRange.cell,FoV);
this_file = sprintf('sD_oD_VObs_rng%02d_fov%03d_orn%03d.mat',...
            SensorRange.cell,FoV,this_orn);
load([this_dir,this_file],'sD','oD','VObs');

% sensorDom := sensor domain translation w.r.t sensor position.
sensorDom = [sD(:,1)+round(sensing_conf(1)-0.5),...
             sD(:,2)+round(sensing_conf(2)-0.5)];

% obstacDom := obstalce domain translated w.r.t sensor's
% position. 
obstacDom = [oD(:,1)+round(sensing_conf(1)-0.5),...
             oD(:,2)+round(sensing_conf(2)-0.5)];

% VOBS := Visibility vs Obstacles matrix for jth FoV.
VOBS = VObs;

% visibile cells 
[ visDom ] = fVisibileDomain( map_env,sensorDom,obstacDom,VOBS );

cellsScanned = visDom;


end

%------------------------ End of the document --------------------------

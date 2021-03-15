function [V] = fVisbilityMatrixDetection( n,...
                                          f,...
                                          c,...
                                          FoVfaces,...
                                          SensorRange,...
                                          FoV,...
                                          map_conf,...
                                          map_env )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               V
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('.. Visibility Matrix')
V = sparse(n,c); % initializing
% V = zeros(n,c); % initializing
     
% for each orientation
for j = 1:f
    
    this_orn = FoVfaces.lgt(j);
    this_orn = this_orn + (this_orn==0)*360;
    sD   = [];
    oD   = [];
    VObs = [];
    %this_dir = sprintf(['sD_oD_VObs/sensing_range_%02dcell/fov_%03ddeg/',...
    %                    'orientation_%03ddeg/'],SensorRange.cell,FoV,this_orn);
    this_dir = sprintf('visibility/templates/sensing_range_%02dcell/fov_%03ddeg/',...
        SensorRange.cell,FoV);
    this_file = sprintf('sD_oD_VObs_rng%02d_fov%03d_orn%03d.mat',...
                SensorRange.cell,FoV,this_orn);
    load([this_dir,this_file],'sD','oD','VObs');
    % For each configuration of a cell.
    for i = 1:n
        if map_conf(i) % if its not occupied
            [rr,cc] = ind2sub(size(map_conf),i);
            sensor = [rr,cc];
            % sensorDom := sensor domain translation w.r.t sensor position.
            %sensorDom = [sD{j}(:,1)+sensor(1) sD{j}(:,2)+sensor(2)];
            sensorDom = [sD(:,1)+sensor(1) sD(:,2)+sensor(2)];
            % obstacDom := obstalce domain translated w.r.t sensor's
            % position.
            %obstacDom = [oD{j}(:,1)+sensor(1) oD{j}(:,2)+sensor(2)];
            obstacDom = [oD(:,1)+sensor(1) oD(:,2)+sensor(2)];
            % VOBS := Visibility vs Obstacles matrix for jth FoV.
            %VOBS = VObs{j};
            VOBS = VObs;
            % visibile cells 
            [ visDom ] = fVisibileDomain( map_env,sensorDom,obstacDom,VOBS );
            % Update Visibility matrix.
            V(visDom,((i-1)*f)+j) = 1;
        end
    end
end
    

% sensing config can not visualize own cell
for i = 1:numel(map_conf)    
    V(i,(i-1)*f+1:(i-1)*f+f) = 0;
end



end

function fFixedMeasurements( FoV,...
                             dir_,...
                             para_,...
                             visualize )


                         
%/////////////////////////////////////
% -- directories
%/////////////////////////////////////
dir_.logs = dir_.TomographyLogsFix;

%/////////////////////////////////////
% -- load previous results
%/////////////////////////////////////

% -- preprocessing
load([dir_.TomographyMaps,'preprocess_maps.mat'],...
    'hotspotSEQ',...    
    'OBSenv',...
    'OBSconf',...
    'map_coverage',...
    'map_env',...
    'SensorRange');
% -- fusion
load([dir_.TomographyFusion,'spp_fusion.mat'],'refPointDist','selectedConfFusion_ioh')
% -- tsp
load([dir_.TomographyResults,'tsp_tomography_spp.mat'],'confSequenceNum');



%////////////////////////////////////////////////////////
% -- global list of planned conf by sequence
%////////////////////////////////////////////////////////

% -- add sequence number to the fused conf
selectedConfPlanned_sioh = [confSequenceNum,selectedConfFusion_ioh];
% -- sorted by sequence number
selectedConfPlanned_sioh = sortrows(selectedConfPlanned_sioh,1);



if visualize == 1

    % - ground truth        
    %[map_conc_______gt] = fGroundTruth( map_env,...
    %                                    para_,...
    %                                    visualize,...
    %                                    dir_);
    %gt_filename = sprintf('map_conc_%d.dat',para_.GroundTruthTimeStamp);
    gt_filename = sprintf('map_con.dat');
    map_conc_______gt = dlmread([dir_.GroundTruthLogs,gt_filename]);


    h = figure; hold on;
    %imshow(map_conc_______gt,'InitialMagnification',800); hold on;
    %set(gca,'YDir','normal')
    %colormap('hot'); % not go back to gray scale.

    %{
    image(map_conc_______gt,'CDataMapping','scaled')
    colormap hot;
    colormap(flipud(colormap))
    %colorbar;
    xlim([0,size(map_conc_______gt,1)])
    ylim([0,size(map_conc_______gt,2)])
    %}

    % {
    imshow(map_env','InitialMagnification',500); hold on;
    set(gca,'YDir','normal')
    colormap('gray'); % not go back to gray scale.
    caxis([-2.75 1]) % for sensor placement solution
    %}

    % -- ground truth
    % {
    gt_colormap = flipud(hot(256));
    max_concentration = max(map_conc_______gt(:));
    delta = max_concentration/(size(gt_colormap,1)-1);    

    for i = 1:size(map_conc_______gt,1)
        for j = 1:size(map_conc_______gt,2)
            % if positive concentration

            if map_conc_______gt(i,j)>=0

                linear_color_number = round( min(map_conc_______gt(i,j),...
                    max_concentration)/delta)+1;
                col = gt_colormap( linear_color_number,:);                

                %----- to avoid black borders of earch cell -----
                % {
                if sum(col) == 3
                    col = col-(1e-10);
                end
                plot(i,j,...
                    's',...
                    'MarkerEdgeColor',col,...
                    'MarkerFaceColor',col,...
                    'MarkerSize',2.5);
            end
        end
    end
    %}

    set(gca, ...
          'Box'         , 'on'      , ...
          'TickLength'  , [.0 .0]   , ...
          'XMinorTick'  , 'off'     , ...
          'YMinorTick'  , 'off'     , ...
          'XGrid'       , 'off'     , ...
          'YGrid'       , 'off'     , ...
          'XtickLabel'  , []        , ...
          'YtickLabel'  , []        );

    color_hotspots = [1,1,0;...
                      0,1,0;...
                      1,0,1];

end




% -- sensing conf to execute
[confs_x,confs_y] = ind2sub(size(map_env),selectedConfPlanned_sioh(:,2));
confs_a = selectedConfPlanned_sioh(:,3);

% -- other sensing parameters
conf_fov = FoV;
conf_r   = SensorRange.cell;
conf_frq = FoV*2;
conf_tim = para_.GroundTruthTimeStamp;

sensing_confs = [confs_x,...
                 confs_y,...
                 confs_a,...
                 conf_fov*ones(size(confs_x,1),1),...
                 conf_r*ones(size(confs_x,1),1),...
                 conf_frq*ones(size(confs_x,1),1),...
                 conf_tim*ones(size(confs_x,1),1)];
             
             
             
             
switch para_.SensingSystem
        
        %=======================================
    case 'robot'
        %=======================================

        sensing_confs = sensing_confs(:,1:3);
        conf_sequence_num = size(sensing_confs,1);
        %para_.RobotPlanFileName = 'plan_two-step-exploration-tomography-fixed.txt';
        %para_.RobotPlanFileName = 'robot_plan.dat';
        dir_.ROSLogsThis = sprintf('%s%s/',dir_.ROSLogs,...
            'two-step-exploration-tomography-fixed');
        dir_.logs;
        ros_service_conf = rossvcclient(para_.ROSServiceExecutedConfNum);

        fRobotMeasurement( sensing_confs,...
                           conf_sequence_num,...
                           ros_service_conf,...
                           dir_,...
                           para_ );


        %=======================================
    case 'simulation'
        %=======================================


        % for conf_num = 1:size(selectedConfPlanned_sioh,1)
        for conf_num = 1:size(sensing_confs,1)

            disp('************************************************************')
            disp('************************************************************')
            disp(' CONFIGURATION ')
            disp('************************************************************')
            disp('************************************************************')



            % % -- sensing conf to execute
            % [conf_x,conf_y] = ind2sub(size(map_env),selectedConfPlanned_sioh(conf_num,2));
            % conf_a = selectedConfPlanned_sioh(conf_num,3);
            % 
            % % -- other sensing parameters
            % conf_fov = FoV;
            % conf_r   = 15;
            % conf_frq = 90;
            % conf_tim = para_.GroundTruthTimeStamp;
            % 
            % sensing_conf = [conf_x,conf_y,conf_a,conf_fov,conf_r,conf_frq,conf_tim];

            sensing_conf = sensing_confs(conf_num,:);

            % -- output measurement file name
            %measu_file = sprintf('measurement_conf%02d_fixed.dat',conf_num);
            measu_file = sprintf('measurement_conf%02d.dat',conf_num);

            % -- if output file does not exist
            %fname_out = [dir_.meas{time_stamp_number},measu_file];
            %if exist(fname_out, 'file') == 0                
                % -- collect integral measurements
                %fprintf('Remote Sensor:\n');
                %fprintf('Set# %02d, %s:\n',i, conf_file);
                %fMeasurements( sensing_conf,...
                %               measu_file,...
                %               para_,...
                %               dir_,...
                %               visualize );
                fSimulatedMeasurement( sensing_conf,...
                                       measu_file,...
                                       map_env,...
                                       para_,...
                                       dir_,...
                                       visualize );
            %end
            fprintf('- Measurements collected for conf #%02d.\n',conf_num);




            if visualize == 1
            conf_num_hot = mod(conf_num,para_.NumberOfConf_SPPH)+...
                (~mod(conf_num,para_.NumberOfConf_SPPH)*para_.NumberOfConf_SPPH);

            hotspot_num = find(selectedConfPlanned_sioh(conf_num,4:end));
            [hot_r,hot_c] = ind2sub(size(map_env),hotspotSEQ(hotspot_num));
            hotspot_current = [hot_r,hot_c];

            fPublish_MeasurementProcess( conf_num,...
                                         conf_num_hot,...
                                         measu_file,...
                                         map_env,...
                                         hotspot_current,...
                                         selectedConfPlanned_sioh,...
                                         selectedConfPlanned_sioh,...                                     
                                         dir_ )

            end
            fprintf('\n');
        end

        if visualize == 1
        print_file_name = sprintf([dir_.logs,'sampling_fixed_coverage']);
        print(print_file_name, '-painters', '-dpdf', '-r1200');
        end
            
end



end

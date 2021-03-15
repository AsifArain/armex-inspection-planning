function [ hc_num,...
           hc_prev,...
           hc_block_num,...
           confs_executed,...
           conf_sequence_num,...
           confs_planned,...
           cellsToBeScanned,...
           retrieve1SE_pts,...
           hcs_list_sub ] = fAdaptiveReplanning1t( FoV,...
                                                   map_env,...
                                                   f,...
                                                   o,...
                                                   FoVfaces,...
                                                   E,...
                                                   T,...
                                                   oV,...
                                                   confKept,...
                                                   hc_num,...
                                                   hc_prev,...
                                                   hc_block_num,...
                                                   SensorRange,...
                                                   confs_executed,...
                                                   confs_executed_hc,...
                                                   conf_sequence_list_hc,...
                                                   conf_sequence_num,...
                                                   confs_planned,...
                                                   cellsToBeScanned,...
                                                   hcs_list_sub,...
                                                   OBSenv,...
                                                   OBSconf,...
                                                   cellsize_env,...
                                                   retrieve1SE_pts,...
                                                   dir_,...
                                                   para_,...
                                                   visualize )

%
%
%


%----- TOMOGRAPHY PART IS MOVED TO ANOTHER FILE
% {

flag_tomo = 1; % flag to perform tomography
reselectedConf2Exe_ind = [];
confNumToExecute_thisSolution = 1;


%/////////////////////////////////////////////////////////////
%
%   ESTIMATE HC
%
%/////////////////////////////////////////////////////////////
% fprintf('\n\n------ INITIAL ESTIMATED HCs IN THIS PART OF TOMOGRAPHY ------\n')
fprintf('\n\n------ INITIALLY ESTIMATED HCs ------\n')
%fprintf('--> Estimate hotspot\n');

para_.HotspotReEstimationMethod = 'ClusteringBased';
% para_.HotspotReEstimationMethod = 'SimpleWeightedMean';

[ hcs_estimated,...
  recon_local_x,...
  recon_local_y,...
  map_recon_local ] = fReEstimatedHcs1t( conf_sequence_list_hc,...
                                         SensorRange,...
                                         OBSenv,...
                                         para_,...
                                         dir_,...
                                         visualize );

fprintf('Num of initial estimated hcs (before check): %02d\n',...
    size(hcs_estimated,1));

%----------------------------------
% plot initially estimated hcs
%----------------------------------
H1 = plot(hcs_estimated(:,1),hcs_estimated(:,2),...
    'ob','MarkerSize',6,'MarkerFaceColor','y');
%H2 = plot(hcs_estimated(:,1),hcs_estimated(:,2),'xb','MarkerSize',6);
%pause;
% delete(H1)
%delete(H2)
if strcmp(para_.SensingSystem,'simulated-sampling')
    %pause
end
%------------------

%*************************************************************************
% Remove the estimated hcs close enough to the existing hcs
%*************************************************************************
if size(hcs_list_sub,1) >= 1
        
    pathDists = fPathDist2(hcs_estimated,hcs_list_sub,map_env);
    euclDists = pdist2(hcs_estimated,hcs_list_sub,'euclidean');

    [minPath,indPath_c] = min(pathDists,[],2);
    indPath_r = (1:size(pathDists,1))';
    indPath = sub2ind(size(pathDists),indPath_r,indPath_c);
    minEucl = euclDists(indPath);
    path2euclRatio = minPath./minEucl;
    
    cond_vec = minPath>=para_.ResEstimatedHcOffset |...
        path2euclRatio>=para_.ReEstimatedHcpath2euclRatio;
    hcs_estimated = hcs_estimated(cond_vec,:);
       
    
    %dist_estimated2all = fPathDist2(hcs_estimated,hcs_list_sub,map_env);    
    %min_dist_estimated2all = min(dist_estimated2all,[],2);
    %hcs_estimated = hcs_estimated(min_dist_estimated2all>=para_.ReplanningThreshould,:);
    
else
    %min_dist_estimated2all = inf;
    %hcs_estimated;
end
%fprintf('Num of initial estimated hcs (after check): %02d\n',size(hcs_estimated,1));
fprintf('Num of initial estimated hcs: %02d\n',size(hcs_estimated,1));
% pause

%/////////////////////////////////////////////////////////////
%
%   IF THIS HOTSPOT IS NOT ALREADY RECONSTRUCTED
%
%/////////////////////////////////////////////////////////////
%if min_dist_estimated2all >= SensorRange.cell
% if any(min_dist_estimated2all >= SensorRange.cell)


%-- Data backup/retrieval point
%--------------------------------------
pt_num = 2;
if any(retrieve1SE_pts == pt_num)
    fprintf("\nData retrieval point %02d....\n",pt_num);
    load([dir_.str,sprintf('LastSave1SE_Pnt%02d',pt_num)]);
    retrieve1SE_pts(retrieve1SE_pts == pt_num) = [];
else
    fprintf("\nData backup point %02d....\n",pt_num);
    save([dir_.str,sprintf('LastSave1SE_Pnt%02d',pt_num)],...
            '-regexp','^(?!(retrieve1SE_pts)$).');
    
    dlmwrite([dir_.str,'BackupPtsSequence1SE.dat'],...
         pt_num,'-append','delimiter',' ');
    fileID = fopen([dir_.str,'BackupPtsSequence1SE.dat'],'a'); 
    fclose(fileID);
end
% pause

if size(hcs_estimated,1) >= 1
    
    
    %-- Data backup/retrieval point
    %--------------------------------------
    pt_num = 3;
    if any(retrieve1SE_pts == pt_num)
        fprintf("\nData retrieval point %02d....\n",pt_num);
        load([dir_.str,sprintf('LastSave1SE_Pnt%02d',pt_num)]);
        retrieve1SE_pts(retrieve1SE_pts == pt_num) = [];
    else
        fprintf("\nData backup point %02d....\n",pt_num);
        save([dir_.str,sprintf('LastSave1SE_Pnt%02d',pt_num)],...
            '-regexp','^(?!(retrieve1SE_pts)$).');
        
        dlmwrite([dir_.str,'BackupPtsSequence1SE.dat'],...
             pt_num,'-append','delimiter',' ');
        fileID = fopen([dir_.str,'BackupPtsSequence1SE.dat'],'a'); 
        fclose(fileID);
    end
    %pause
    
    %--------------------------    
    confNumToExecute_thisSolution = 1;
    hcs_list_sub = [hcs_list_sub; hcs_estimated];
    
    %*******************************
    %   Hcs to publish
    %*******************************
    switch para_.SensingSystem
        case 'robot'
            %-- create hcs file, if not already exists
            if exist([dir_.ROSLogsThis,para_.RobotHcsFileName],'file') == 0    
                fileID = fopen([dir_.ROSLogsThis,para_.RobotHcsFileName],'wt'); 
                fclose(fileID);
            end
            %-- write hcs file
            dlmwrite([dir_.ROSLogsThis,para_.RobotHcsFileName],...
                     hcs_list_sub,'delimiter',' ');
            fileID = fopen([dir_.ROSLogsThis,para_.RobotHcsFileName],'a'); 
            fclose(fileID);

            %-- copy
            thisFile = ([dir_.ROSLogsThis,para_.RobotHcsFileName]);
            thatFile = ([dir_.logs,para_.RobotHcsFileName]);
            syncCommand = sprintf('rsync -raz %s %s',thisFile,thatFile);
            system(syncCommand);
            
        case 'simulated-sampling'
            % do nothing
    end
    %----------------------------------------------
        
    fprintf('\n\n ****************************************************\n');
    fprintf('                  ADAPTIVE TOMOGRAPHY \n');
    fprintf(' ****************************************************\n');
    fprintf('Num of estimated hcs for tomography planning: %02d\n',size(hcs_estimated,1));
    
                                          
    %-- Data backup/retrieval point
    %--------------------------------------
    pt_num = 4;
    if any(retrieve1SE_pts == pt_num)
        fprintf("\nData retrieval point %02d....\n",pt_num);
        load([dir_.str,sprintf('LastSave1SE_Pnt%02d',pt_num)]);
        retrieve1SE_pts(retrieve1SE_pts == pt_num) = [];
    else
        
        [ hc_num,...
          reselectedConf2Exe_ind,...
          reselectedConf2Exe_orn] = fPlan_1tARMEx_GDM( hc_num,...
                                                       hcs_estimated,...
                                                       hc_block_num,...
                                                       conf_sequence_num,...
                                                       conf_sequence_list_hc,...
                                                       confs_executed,...
                                                       map_env,...
                                                       recon_local_x,...
                                                       recon_local_y,...
                                                       map_recon_local,...
                                                       FoV,...
                                                       FoVfaces,...
                                                       confs_executed_hc,...
                                                       OBSenv,...
                                                       OBSconf,...
                                                       SensorRange,...
                                                       cellsize_env,...
                                                       dir_,...
                                                       para_ );
        
        
        
        fprintf("\nData backup point %02d....\n",pt_num);
        save([dir_.str,sprintf('LastSave1SE_Pnt%02d',pt_num)],...
            '-regexp', '^(?!(retrieve1SE_pts)$).')
        
        dlmwrite([dir_.str,'BackupPtsSequence1SE.dat'],...
             pt_num,'-append','delimiter',' ');
        fileID = fopen([dir_.str,'BackupPtsSequence1SE.dat'],'a'); 
        fclose(fileID);
    end
    %pause
        
    hc_block_num = hc_block_num+1;
    %--------------------------
    fprintf('\n-- Num of confs in the list: %02d \n',numel(reselectedConf2Exe_ind));
    
    cond_PossibleHcBasedOnLastConfInList = 1;

    %while (numel(reselectedConf2Exe_ind)>=1 || flag_tomo)
    while ( numel(reselectedConf2Exe_ind)>=1 ||...
            cond_PossibleHcBasedOnLastConfInList ||...
            flag_tomo )

        
        
        %-- Data backup/retrieval point
        %--------------------------------------
        pt_num = 5;
        if any(retrieve1SE_pts == pt_num)
            fprintf("\nData retrieval point %02d....\n",pt_num);
            load([dir_.str,sprintf('LastSave1SE_Pnt%02d',pt_num)]);
            retrieve1SE_pts(retrieve1SE_pts == pt_num) = [];
        else
            fprintf("\nData backup point %02d....\n",pt_num);
            save([dir_.str,sprintf('LastSave1SE_Pnt%02d',pt_num)],...
                '-regexp', '^(?!(retrieve1SE_pts)$).')
            
            dlmwrite([dir_.str,'BackupPtsSequence1SE.dat'],...
                 pt_num,'-append','delimiter',' ');
            fileID = fopen([dir_.str,'BackupPtsSequence1SE.dat'],'a'); 
            fclose(fileID);
        end
        %pause
        
        %loop_num = loop_num + 1;
        %sprintf('************ LOOP TEST #%d *************',loop_num)
        flag_tomo = 0;

        %/////////////////////////////////////////////////////////////
        %
        %   ESTIMATE HC
        %
        %/////////////////////////////////////////////////////////////
        %fprintf('**************** ESTIMATED HOTSPOT *******************\n')
        %fprintf('\n\n------ ESTIMATED HCs DURING TOMOGRAPHIC PLAN EXECUTION ------\n');
        fprintf('\n\n------ ESTIMATED HCs DURING TOMOGRAPHY ------\n');
        
        % -- re-estimated hotspot
        % ----------------------------
        [ hcs_reestimated,...
          recon_local_x,...
          recon_local_y,...
          map_recon_local ] = fReEstimatedHcs1t( conf_sequence_list_hc,...
                                                  SensorRange,...
                                                  OBSenv,...
                                                  para_,...
                                                  dir_,...
                                                  visualize );
                                              
        %hcs_reestimated                                      
        fprintf('Num of estimated hcs (before check): %02d\n',...
            size(hcs_reestimated,1));
        
        %-------------------------------
        %-- plot
        %-------------------------------
        H1 = plot(hcs_reestimated(:,1),hcs_reestimated(:,2),...
            'ob','MarkerSize',6,'MarkerFaceColor','y');
        %H2 = plot(hcs_reestimated(:,1),hcs_reestimated(:,2),'xb','MarkerSize',6);
        %pause;
        delete(H1)
        %delete(H2)
        if strcmp(para_.SensingSystem,'simulated-sampling')
            %pause
        end
        %--------------------------------
        
        %*************************************************************************
        % Remove the estimated hcs close enough to the existing hcs
        %*************************************************************************
        if size(hcs_list_sub,1) >= 1
            
            
            pathDists = fPathDist2(hcs_reestimated,hcs_list_sub,map_env);
            euclDists = pdist2(hcs_reestimated,hcs_list_sub,'euclidean');

            [minPath,indPath_c] = min(pathDists,[],2);
            indPath_r = (1:size(pathDists,1))';
            indPath = sub2ind(size(pathDists),indPath_r,indPath_c);
            minEucl = euclDists(indPath);
            path2euclRatio = minPath./minEucl;

            cond_vec = minPath>=para_.ResEstimatedHcOffset |...
                path2euclRatio>=para_.ReEstimatedHcpath2euclRatio;
            hcs_reestimated = hcs_reestimated(cond_vec,:);


            %dist_estimated2all = fPathDist2(hcs_estimated,hcs_list_sub,map_env);    
            %min_dist_estimated2all = min(dist_estimated2all,[],2);
            %hcs_estimated = hcs_estimated(min_dist_estimated2all>=para_.ReplanningThreshould,:);
            
            %{           
            %dist_estimated2prv = fPathDist(hc_estimated,hc_prev,map_env)
            dist_estimated2all = fPathDist2(hcs_reestimated,hcs_list_sub,map_env);
            %min_dist_estimated2all = min(dist_estimated2all);
            min_dist_estimated2all = min(dist_estimated2all,[],2);            
            hcs_reestimated = ...
                hcs_reestimated(min_dist_estimated2all>para_.ReplanningThreshould,:);
            %}
            
        end        
        %fprintf('Num of estimated hcs (after check): %02d\n',size(hcs_reestimated,1));
        fprintf('Num of newly estimated hcs: %02d\n',size(hcs_reestimated,1));
        %pause
        
        
        
        %-- Data backup/retrieval point
        %--------------------------------------
        pt_num = 6;
        if any(retrieve1SE_pts == pt_num)
            fprintf("\nData retrieval point %02d....\n",pt_num);
            load([dir_.str,sprintf('LastSave1SE_Pnt%02d',pt_num)]);
            retrieve1SE_pts(retrieve1SE_pts == pt_num) = [];
        else
            fprintf("\nData backup point %02d....\n",pt_num);
            save([dir_.str,sprintf('LastSave1SE_Pnt%02d',pt_num)],...
                '-regexp', '^(?!(retrieve1SE_pts)$).')
            
            dlmwrite([dir_.str,'BackupPtsSequence1SE.dat'],...
                 pt_num,'-append','delimiter',' ');
            fileID = fopen([dir_.str,'BackupPtsSequence1SE.dat'],'a'); 
            fclose(fileID);
        end
        %pause
        
        %/////////////////////////////////////////////////////////////
        %
        %   IF DISTANCE BETWEEN THE ESTIMATED HC AND THE PREVIOUS
        %   IS GREATER THAN A THERSHOULD
        %
        %/////////////////////////////////////////////////////////////
        %if dist_estimated2prv > para_.ReplanningThreshould
        
        %/////////////////////////////////////////////////////////////
        %
        %   IF THERE IS A NEWLY ESITMATED HCs
        %
        %/////////////////////////////////////////////////////////////
        
        
        
        
        %-- Data backup/retrieval point
        %--------------------------------------
        pt_num = 8;
        if any(retrieve1SE_pts == pt_num)
            fprintf("\nData retrieval point %02d....\n",pt_num);
            load([dir_.str,sprintf('LastSave1SE_Pnt%02d',pt_num)]);
            retrieve1SE_pts(retrieve1SE_pts == pt_num) = [];
        else
            
            
            %****************************************************************                        
            if size(hcs_reestimated,1) >= 1
                
                %-- Data backup/retrieval point
                %--------------------------------------
                pt_num = 7;
                if any(retrieve1SE_pts == pt_num)
                    fprintf("\nData retrieval point %02d....\n",pt_num);
                    load([dir_.str,sprintf('LastSave1SE_Pnt%02d',pt_num)]);
                    retrieve1SE_pts(retrieve1SE_pts == pt_num) = [];
                else
                    fprintf("\nData backup point %02d....\n",pt_num);
                    save([dir_.str,sprintf('LastSave1SE_Pnt%02d',pt_num)],...
                        '-regexp', '^(?!(retrieve1SE_pts)$).')

                    dlmwrite([dir_.str,'BackupPtsSequence1SE.dat'],...
                         pt_num,'-append','delimiter',' ');
                    fileID = fopen([dir_.str,'BackupPtsSequence1SE.dat'],'a'); 
                    fclose(fileID);
                end
                %pause

                confNumToExecute_thisSolution = 1;
                hcs_list_sub = [hcs_list_sub; hcs_reestimated];
                hcs_estimated = [hcs_estimated;hcs_reestimated];


                %*******************************
                %   Hcs to publish
                %*******************************
                switch para_.SensingSystem
                    case 'robot-sampling'
                        %-- create hcs file, if not already exists
                        if exist([dir_.ROSLogsThis,para_.RobotHcsFileName],'file') == 0    
                            fileID = fopen([dir_.ROSLogsThis,para_.RobotHcsFileName],'wt'); 
                            fclose(fileID);
                        end
                        %-- write hcs file
                        dlmwrite([dir_.ROSLogsThis,para_.RobotHcsFileName],...
                                 hcs_list_sub,'delimiter',' ');
                        fileID = fopen([dir_.ROSLogsThis,para_.RobotHcsFileName],'a'); 
                        fclose(fileID);

                        %-- copy
                        thisFile = ([dir_.ROSLogsThis,para_.RobotHcsFileName]);
                        thatFile = ([dir_.logs,para_.RobotHcsFileName]);
                        syncCommand = sprintf('rsync -raz %s %s',thisFile,thatFile);
                        system(syncCommand);

                    case 'simulated-sampling'
                        % do nothing
                end
                %----------------------------------------------
                fprintf('\n\n ****************************************************\n');
                fprintf('                  ADAPTIVE TOMOGRAPHY \n');
                fprintf(' ****************************************************\n');


                fprintf('Num of estimated hcs for tomography planning: %02d\n',size(hcs_estimated,1));
                [ hc_num,...
                  reselectedConf2Exe_ind,...
                  reselectedConf2Exe_orn] = fPlan_1tARMEx_GDM( hc_num,...
                                                               hcs_estimated,...
                                                               hc_block_num,...
                                                               conf_sequence_num,...
                                                               conf_sequence_list_hc,...
                                                               confs_executed,...
                                                               map_env,...
                                                               recon_local_x,...
                                                               recon_local_y,...
                                                               map_recon_local,...
                                                               FoV,...
                                                               FoVfaces,...
                                                               confs_executed_hc,...
                                                               OBSenv,...
                                                               OBSconf,...
                                                               SensorRange,...
                                                               cellsize_env,...
                                                               dir_,...
                                                               para_ );
                hc_block_num = hc_block_num+1;
            end
            %****************************************************************
            
            
            
            
            
            
            fprintf("\nData backup point %02d....\n",pt_num);
            save([dir_.str,sprintf('LastSave1SE_Pnt%02d',pt_num)],...
                '-regexp', '^(?!(retrieve1SE_pts)$).')
            
            dlmwrite([dir_.str,'BackupPtsSequence1SE.dat'],...
                 pt_num,'-append','delimiter',' ');
            fileID = fopen([dir_.str,'BackupPtsSequence1SE.dat'],'a'); 
            fclose(fileID);
        end
        %pause
        
        
        % --- plot reselected conf 
        % ----------------------------
        [ h_reselct ] = fPlotReselectedConf( reselectedConf2Exe_ind,...
                                             reselectedConf2Exe_orn,...
                                             map_env,...
                                             FoV,...
                                             para_ );
        %/////////////////////////////////////////////////////////////
        %
        %  SENSING CONF FOR MEASUREMENTS
        %
        %/////////////////////////////////////////////////////////////

        %fprintf('\n******** measurement conf #%02d \n',conf_sequence_num);
        
        fprintf('\n-- Num of confs in the list: %02d \n',numel(reselectedConf2Exe_ind));

        [ reselectedConf2Exe_ind,...
          reselectedConf2Exe_orn,...
          sensing_conf] = fConf2Execute( reselectedConf2Exe_ind,...
                                         reselectedConf2Exe_orn,...
                                         SensorRange,...
                                         map_env,...
                                         FoV,...
                                         dir_,...
                                         para_ );


        conf_sequence_num = conf_sequence_num+1;
        conf_sequence_list_hc = [conf_sequence_list_hc,conf_sequence_num];
        confNumToExecute_thisSolution = confNumToExecute_thisSolution+1;
        
        % conf for
        sensing_conf = [sensing_conf 1]
        
        [ cellsToBeScanned,...
          cellsScanned ] = fExecuteThisConf( sensing_conf,...
                                             conf_sequence_num,...
                                             map_env,...
                                             FoV,...
                                             SensorRange,...
                                             cellsToBeScanned,...
                                             para_,...
                                             dir_ );
        %pause
        confs_executed
        confs_executed    = [confs_executed;sensing_conf]
        %confs_executed    = [confs_executed;[sensing_conf,1]];
        %confs_executed_hc = [confs_executed_hc;sensing_conf];
        confs_executed_hc = [confs_executed_hc;sensing_conf(:,1:7)];
        
        
        % -- plot 
        % ----------------------------
        confFor = 'tomography';
        %[ ppmm ] = fPlotExecutedConf(conf_sequence_num,scannedInd,map_env,confFor,dir_);
        [ ppmm ] = fPublish_ExecutedConf( conf_sequence_num,...
                                          sensing_conf,...
                                          cellsScanned,...
                                          SensorRange,...
                                          map_env,...
                                          confFor,...
                                          FoV,...
                                          para_,...
                                          dir_ );
        %pause
        
        
        
        
        
        %/////////////////////////////////////////////////////////////
        %
        %   CHECK: GOING BACK TO CONTINUE ESTIMATING NEW HC(S) (AND 
        %   REPLANNING) IF THIS IS THE LAST CONF IN THE LIST AND A NEW
        %   ESTIMATION OF HC(S) IS POSSIBLE.
        %
        %/////////////////////////////////////////////////////////////
        
        
        if numel(reselectedConf2Exe_ind)<1
                        
            fprintf('\n-- ESTIMATED HCs FOR THE LAST CONF IN THE TOMOGRAPHIC LIST --\n');

            % -- re-estimated hotspot
            % ----------------------------
            [ hcs_estimated_cond,~,~,~ ] = fReEstimatedHcs1t( conf_sequence_list_hc,...
                                                              SensorRange,...
                                                              OBSenv,...
                                                              para_,...
                                                              dir_,...
                                                              visualize );
            %{
            fprintf('Num of conditional hcs (before check): %02d\n',size(hcs_estimated_cond,1));
            H1 = plot(hcs_estimated_cond(:,1),hcs_estimated_cond(:,2),'ob','MarkerSize',6,'MarkerFaceColor','y');
            H2 = plot(hcs_estimated_cond(:,1),hcs_estimated_cond(:,2),'xb','MarkerSize',6);
            pause;
            %delete(H1)
            delete(H2)
            %}
                                                               
            % -- distances and validity
            % ----------------------------
            if size(hcs_list_sub,1) >= 1
                
                
                pathDists = fPathDist2(hcs_estimated_cond,hcs_list_sub,map_env);
                euclDists = pdist2(hcs_estimated_cond,hcs_list_sub,'euclidean');

                [minPath,indPath_c] = min(pathDists,[],2);
                indPath_r = (1:size(pathDists,1))';
                indPath = sub2ind(size(pathDists),indPath_r,indPath_c);
                minEucl = euclDists(indPath);
                path2euclRatio = minPath./minEucl;

                cond_vec = minPath>=para_.ResEstimatedHcOffset |...
                    path2euclRatio>=para_.ReEstimatedHcpath2euclRatio;
                hcs_estimated_cond = hcs_estimated_cond(cond_vec,:);
                
                %dist_estimated2all = fPathDist2(hcs_estimated,hcs_list_sub,map_env);    
                %min_dist_estimated2all = min(dist_estimated2all,[],2);
                %hcs_estimated = hcs_estimated(min_dist_estimated2all>=para_.ReplanningThreshould,:);
                
                %{
                dist_estimated2all = fPathDist2(hcs_estimated_cond,hcs_list_sub,map_env);
                %min_dist_estimated2all = min(dist_estimated2all);
                min_dist_estimated2all = min(dist_estimated2all,[],2);            
                hcs_estimated_cond = ...
                    hcs_estimated_cond(min_dist_estimated2all>para_.ReplanningThreshould,:);
                %}
                
            end
            %fprintf('Num of conditional hcs (after check): %02d\n',size(hcs_estimated_cond,1));
            fprintf('Num of newly estimated Hcs for possibly last conf: %02d\n',...
                size(hcs_estimated_cond,1));            
            cond_PossibleHcBasedOnLastConfInList = size(hcs_estimated_cond,1)>=1;
            
        else
            cond_PossibleHcBasedOnLastConfInList = 0;
        end
        
        

        % ------------------------------------------------------------ 

        % -- distance from previous hotspot to relocated hotspot
        %dist_prev2relocated = fPathDist(hotspot_prev,hotspot_estimated,map_env);
        hc_prev = hcs_estimated;
        
        
        %-- Data backup/retrieval point
        %--------------------------------------
        pt_num = 9;
        if any(retrieve1SE_pts == pt_num)
            fprintf("\nData retrieval point %02d....\n",pt_num);
            load([dir_.str,sprintf('LastSave1SE_Pnt%02d',pt_num)]);
            retrieve1SE_pts(retrieve1SE_pts == pt_num) = [];
        else
            fprintf("\nData backup point %02d....\n",pt_num);
            save([dir_.str,sprintf('LastSave1SE_Pnt%02d',pt_num)],...
                '-regexp', '^(?!(retrieve1SE_pts)$).')
            
            dlmwrite([dir_.str,'BackupPtsSequence1SE.dat'],...
                 pt_num,'-append','delimiter',' ');
            fileID = fopen([dir_.str,'BackupPtsSequence1SE.dat'],'a'); 
            fclose(fileID);
        end
        %pause

    end


    %hc_num = hc_num+1;


    %//////////////////////////////////////////////////////////
    %
    %   REPLANNING FOR COVERAGE
    %
    %//////////////////////////////////////////////////////////
    %fprintf('\n*************** REPLANNING FOR COVERAGE ***************\n');
    fprintf('\n\n ****************************************************\n');
    fprintf('             REPLANNING FOR GAS DETECTION \n');
    fprintf(' ****************************************************\n');
    
    [ C,...
      confSequenceNum,...
      coverageCells ] = fPlan_1tARMEx_GD( hc_num,...
                                          f,...
                                          o,...
                                          FoVfaces,...
                                          cellsToBeScanned,...
                                          confs_executed,...
                                          E,...
                                          T,...
                                          map_env,...
                                          OBSenv,...
                                          dir_,...
                                          para_,...
                                          visualize );
    
                                          
    %-- Data backup/retrieval point
    %--------------------------------------
    pt_num = 10;
    if any(retrieve1SE_pts == pt_num)
        fprintf("\nData retrieval point %02d....\n",pt_num);
        load([dir_.str,sprintf('LastSave1SE_Pnt%02d',pt_num)]);
        retrieve1SE_pts(retrieve1SE_pts == pt_num) = [];
    else
        fprintf("\nData backup point %02d....\n",pt_num);
        save([dir_.str,sprintf('LastSave1SE_Pnt%02d',pt_num)],...
            '-regexp', '^(?!(retrieve1SE_pts)$).')
        
        dlmwrite([dir_.str,'BackupPtsSequence1SE.dat'],...
             pt_num,'-append','delimiter',' ');
        fileID = fopen([dir_.str,'BackupPtsSequence1SE.dat'],'a'); 
        fclose(fileID);
    end
    %pause
    
    %//////////////////////////////////////////////////////////
    %
    %   UPDATED COVERAGE PLAN
    %
    %//////////////////////////////////////////////////////////
    [confs_detection_iterative] = fExecutableConfsGD( FoVfaces,...
                                                      SensorRange,...
                                                      FoV,...
                                                      oV,...
                                                      C,...
                                                      confSequenceNum,...
                                                      confKept,...
                                                      map_env,...
                                                      dir_,...
                                                      para_ );
    % -- plot
    % ----------------------------------------
    % FIXME
    % -- commented for smokebot only    
    %{
    fPlotDetectionPlan1SE( confs_detection_iterative,...
                           FoV,...
                           SensorRange,...
                           coverageCells,...
                           map_env,...
                           para_ );
    %}
    
    %coverageCells = [];
    fPublish_GDPlan1t( confs_detection_iterative,...
                       FoV,...
                       SensorRange,...
                       coverageCells,...
                       map_env,...
                       para_ );
    
                       
    % -- update list of planned confs
    % ----------------------------------------
    confs_planned  = confs_detection_iterative;


end
conf_sequence_num = conf_sequence_num+1; %FIXME

%}


end






function [ C,...
           confSequenceNum,...
           coverageCells ] = fPlan_1tARMEx_GD( hc_num,...
                                               f,...
                                               o,...
                                               FoVfaces,...
                                               cellsToBeScanned,...
                                               confs_executed,...
                                               E,...
                                               T,...
                                               map_env,...
                                               OBSenv,...
                                               dir_,...
                                               para_,...
                                               visualize )
%
%*******************************************************
%    RE-PLANNING FOR GAS DETECTION
%*******************************************************
%
% @AsifArain 08-Apr-2018 

            
%**********************************************************
%
%               SENSOR PLACEMENT PROBLEM
%
%**********************************************************
fprintf('\n****** SPP FULL COVERAGE ******\n');
% -- update visibility matrix
% ----------------------------
load([dir_.sppdetection,'preprocess_detection.mat'],...
    'V','oV','confKept','cellKept','FoVfaces','SensorRange')
V = sparse(oV(cellsToBeScanned,confKept));

% -- spp parameters
% ----------------------------
para_.Algorithm = 'convSPP';  % Options: 'Comb','convSPP'
para_.AlgoRev   = 5;       % Revision.
para_.Prog      = 'GRB';   % Options: 'CVX','GRB','MOT'
para_.numHistIR = 5;       % 
para_.maxIterat = 150;     % Maximum number of iterations.
para_.numCu     = 80;      % Number of uncertain confs. to terminate RWL1.

% optimization
% ----------------------------
%[ C,tSPP ] = fOptCoverageS2( V,FoVfaces,para_ );
[ C,...
  C0,...
  C1,...
  Cu,...
  combC,...
  lBound,...
  uBound,...
  logs_rwl1,...
  tSPP_rwl1,...
  tSPP_Comb,...
  tSPP ] = fPlan_GD_FullCoverage( V,para_ );


%**********************************************************
%
%       CONFS FOR REDUCED COVERAGE 
%
%**********************************************************
fprintf('\n****** SPP REDUCED COVERAGE ******\n');

%-- coverage index based confs
%--------------------------------
% cellsToBeScanned
cellsNOtToBeScanned = cellKept(ismember(cellKept,cellsToBeScanned)==0);
oV(cellsNOtToBeScanned',confKept) = 1;
V = sparse(oV(cellKept,confKept));

% size(cellsToBeScanned)
% size(cellsNOtToBeScanned)
% size(cellKept)

FoV = para_.FieldOfViewDEG;
% visualize = 0
[ C,...
  reducedCoverage ] = fPlan_GD_ReducedCoverage( C,...
                                                 V,...
                                                 oV,...
                                                 FoV,...
                                                 FoVfaces,...
                                                 confKept,...
                                                 SensorRange,...
                                                 map_env,...
                                                 para_,...
                                                 dir_,...
                                                 visualize );
                             
V = sparse(oV(cellKept,confKept));
% V*C
% numel(V*C)
% numel(find((V*C)~=0))
% numel(find((V*C)==0))
percentage_coverage = numel(find((V*C)~=0))/numel(V*C)

oC = zeros(size(oV,2),1);
oC(confKept(C)) = 1;
coverageCells = find(oV*oC);

% pause

% pause
%file_name = sprintf('coverage_s2_hotspot%02d.mat',hc_num);
file_name = sprintf('detection_hc%02d_spp.mat',hc_num);

save([dir_.sppdetection,file_name],'C','reducedCoverage','coverageCells','tSPP');


%##############################################################
%
%               TSP
%
%##############################################################

if strcmp(para_.TravelingType1SE,'TSP2')
    
fprintf('\n******** TSP ********\n');
tsp_time_init = tic;

Conf = zeros(size(oV,2),1);
Conf(confKept) = C;

% -- conf vector for tsp.
ConfTSP = Conf';


% robotStartPosition_ind = sub2ind( size(map_env),...
%                                   round(confs_executed(end,1)-0.5),...
%                                   round(confs_executed(end,2)-0.5) );
% robotStartConf_ind = ((robotStartPosition_ind-1)*f)+1;
% ConfTSP(robotStartConf_ind) = 1;


% -- removing conf on occupied cells
% --------------------------------
% f
% obs-conf nums
confOBS = zeros(numel(OBSenv.ind)*f,1);
for i = 1:numel(OBSenv.ind)
    confOBS(((i-1)*f)+1:((i-1)*f)+f) = ...
        (((OBSenv.ind(i)-1)*f)+1):(((OBSenv.ind(i)-1)*f)+f);
end
ConfTSP(:,confOBS) = [];
ConfTSP = ConfTSP';

% -- TSP PARAMETERS
%-------------------------
para_.TravCost   = 10; %ParamPreprocess.TravCost;
para_.ConfOBS    = 1;
para_.EoTinAstar = 1; %0; 1 for simplified path A*, 0 for complex path A*

[ tRoute,...
  tCost,...
  tTSP ] = fTSP2( o,...
                  FoVfaces,...
                  ConfTSP,...
                  E,...
                  T,...
                  map_env,...
                  OBSenv,...
                  para_ );
              
              
%   disp('after fTSP2,,,,')            
%--------------------------------------------------------------
% the start conf is the nearnest conf to the executed conf
%--------------------------------------------------------------
% Index numbers of Conf.
thisConf = find(ConfTSP);
for i = 1:numel(OBSenv.ind)
    thisConf(thisConf>((OBSenv.ind(i)-1)*f)) =...
        thisConf(thisConf>((OBSenv.ind(i)-1)*f))+f;
end

% thisConf

% Cell number of each Conf.
iCELL = fix(thisConf/f)+(~~mod(thisConf,f)); 
[rCELL,cCELL] = ind2sub(size(map_env),iCELL);
thisConfXY = [rCELL cCELL];

dist = fPathDist2(confs_executed(end,1:2),thisConfXY,map_env);
[~,ind] = min(dist);

% thisConfXY(ind,1)
% thisConfXY(ind,2)

startConf_ind = thisConf(ind); 
%sub2ind(size(map_env),thisConfXY(ind,1),thisConfXY(ind,2))

%---------------------------------------------------------------


% -- traveling route with sequence starting from initial position            
% -------------------------------------
% tRouteOrderedWRTStart = tRoute(1:end-1,:);
% tRouteOrderedWRTStart = circshift( tRouteOrderedWRTStart,...
%     numel(tRouteOrderedWRTStart)-...
%     find(tRouteOrderedWRTStart==robotStartConf_ind)+1);
tRouteOrderedWRTStart = tRoute(1:end-1,:);
tRouteOrderedWRTStart = circshift( tRouteOrderedWRTStart,...
    numel(tRouteOrderedWRTStart)-...
    find(tRouteOrderedWRTStart==startConf_ind)+1);

% -- conf sequence number in the order of conf vector index
% (start position is not a conf) 
% -------------------------------------
[~,confSequenceNum,~] = intersect(tRouteOrderedWRTStart,find(Conf));

% -- since start position is not a conf,
% -------------------------------------
% confSequenceNum = confSequenceNum-1;

% -- save results
% -------------------------------------
file_name = sprintf('detection_hc%02d_tsp.mat',hc_num);
save([dir_.sppdetection,file_name],...
    'tRoute',...
    'tCost',...
    'tTSP',...
    'confSequenceNum',...
    'tRouteOrderedWRTStart');

% -- computation time
% -------------------------------------
tsp_time =  1e-4*(round(toc(tsp_time_init)*1e4));
disp(['Compuation time (TSP): ',num2str(tsp_time),' sec']);
fprintf('\n--------------------------------------------\n');

end


%**************************************************************
%
%               EXPLOITATION ENABLED SORT
%
%**************************************************************

if strcmp(para_.TravelingType1SE,'ExploitationEnabledSort')
    
    fprintf('\n******** EXPLOITATION ENABLED SORT ********\n');
    ees_time_init = tic;

    
    confs_vec = zeros(size(oV,2),1);
    confs_vec(confKept) = C;
    Confs_num = find(confs_vec);
    
    f = FoVfaces.num;
    confCell_ind = fix(Confs_num/f)+(~~mod(Confs_num,f));
    [confCell_r,confCell_c] = ind2sub(size(map_env),confCell_ind);
    
    sensing_confs = [confCell_r,confCell_c];
    initial_position = confs_executed(end,1:2);

    [ confSequenceNum,...
      sensing_confs_sorted ] = fExploitationSortedConfs( sensing_confs,...
                                                         initial_position,...
                                                         map_env );

    % -- save results
    % -------------------------------------
    file_name = sprintf('detection_hc%02d_ees.mat',hc_num);
    save([dir_.sppdetection,file_name],...
        'confSequenceNum',...
        'sensing_confs_sorted');

    % -- computation time
    % -------------------------------------
    ees_time =  1e-4*(round(toc(ees_time_init)*1e4));
    disp(['Compuation time (EES): ',num2str(ees_time),' sec']);
    

end


fprintf('\n--------------------------------------------\n');

end


% {
function [ reselectedConf2Exe_ind,...
           reselectedConf2Exe_orn,...
           sensing_conf ] = fConf2Execute( reselectedConf2Exe_ind,...
                                           reselectedConf2Exe_orn,...
                                           SensorRange,...
                                           map_env,...
                                           FoV,...
                                           dir_,...
                                           para_ )
% 
% 


%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


%-- all confs in the plan
%=================================================
confs_ind = reselectedConf2Exe_ind;
[confs_x,confs_y] = ind2sub(size(map_env),confs_ind);
confs_a = reselectedConf2Exe_orn;


%-- publish all confs in the plan
%=================================================

switch para_.SensingSystem
    
    case 'robot'        
        %
        %-- create plan file, if not already exists
        if exist([dir_.ROSLogsThis,para_.RobotLocalPlanFileName],'file') == 0    
            fileID = fopen([dir_.ROSLogsThis,para_.RobotLocalPlanFileName],'wt'); 
            fclose(fileID);
        end
        %-- write local/current plan
        dlmwrite([dir_.ROSLogsThis,para_.RobotLocalPlanFileName],...
                 [confs_x,... % x-position
                  confs_y,... % y-position
                  confs_a,... % orientation
                  ],...               
                  'delimiter',' ');
        fileID = fopen([dir_.ROSLogsThis,para_.RobotLocalPlanFileName],'a'); 
        fclose(fileID);

        %-- 
        %-- copy command
        thisFile = ([dir_.ROSLogsThis,para_.RobotLocalPlanFileName]);
        thatFile = ([dir_.logs,para_.RobotLocalPlanFileName]);
        syncCommand = sprintf('rsync -raz %s %s',thisFile,thatFile);
        system(syncCommand);
        
    case 'simulated-sampling'
        %
end


%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx






%-- sensing conf to execute
%================================================
%conf_ind = reselectedConf2Exe_ind(conf_num_this_hc) %executedConfReplanned_ioh(end,1);
conf_ind = reselectedConf2Exe_ind(1);
[conf_x,conf_y] = ind2sub(size(map_env),conf_ind);
%conf_a   = reselectedConf2Exe_orn(conf_num_this_hc);
conf_a   = reselectedConf2Exe_orn(1);

% -- other sensing parameters
conf_fov = FoV;
conf_r   = SensorRange.cell;
conf_frq = FoV*2;
conf_tim = para_.GroundTruthTimeStamp;

%sensing_conf = [conf_x+0.5,conf_y+0.5,conf_a,conf_fov,conf_r,conf_frq,conf_tim];
switch para_.SensingSystem
    case 'robot'
        sensing_conf = [conf_x+0.0,conf_y+0.0,conf_a,conf_fov,conf_r,conf_frq,conf_tim];

    case 'simulated-sampling'
        sensing_conf = [conf_x+0.5,conf_y+0.5,conf_a,conf_fov,conf_r,conf_frq,conf_tim];
end

reselectedConf2Exe_ind(1) = [];
reselectedConf2Exe_orn(1) = [];

                
                

end

%}





% {
function [ h_reselct ] = fPlotReselectedConf( reselectedConf2Exe_ind,...
                                              reselectedConf2Exe_orn,...
                                              map_env,...
                                              FoV,...
                                              para_ )
% 
% 


% ----------------------------
%{
for i = 1:numel(selectdConf_r)
    start_angle  = reselectedConf2Exe_orn(i)-FoV/2; % first/start angle.
    sector_angle = FoV;              % increment in the first/start angle.
    MaxPts       = 1000;             % maximum points to plot the FoV.
    [xx,yy]      = SecDraw(start_angle,sector_angle,SensorRange.cell,MaxPts); % x,y points.
    plot(xx+selectdConf_r(i)-0.5,yy+selectdConf_c(i)-0.5,...
        'color',[0.5,0.8,0.8],...
        'LineWidth',2);
end
%}


% col = [235,190,190]/256;
col = [1,0,0];


[selectdConf_r,selectdConf_c] = ind2sub(size(map_env),reselectedConf2Exe_ind);

n = 0;
for i = 1:numel(selectdConf_r)
    
    %-- symoblic fov
    %------------------------------
    if para_.Publish_1tARMEx_SymbolicFOV == 1
        
        start_angle  = reselectedConf2Exe_orn(i)-FoV/2;
        end___angle  = reselectedConf2Exe_orn(i)+FoV/2;
        r            = 2; %5;
        x1           = r*cosd(start_angle);
        y1           = r*sind(start_angle);
        x2           = r*cosd(end___angle);
        y2           = r*sind(end___angle);

        h_reselct(n+1) = plot([selectdConf_r(i)-0.0,x1+selectdConf_r(i)-0.0],...
             [selectdConf_c(i)-0.0,y1+selectdConf_c(i)-0.0],...
              'color',col,'LineWidth',1.5);      
        h_reselct(n+2) = plot([selectdConf_r(i)-0.0,x2+selectdConf_r(i)-0.0],...
             [selectdConf_c(i)-0.0,y2+selectdConf_c(i)-0.0],...
          'color',col,'LineWidth',1.5);
    end
    
    %-- position
    %-----------------------
    if para_.Publish_1tARMEx_ConfPosition == 1
        
        h_reselct(n+3) = plot(selectdConf_r(i)-0.0,selectdConf_c(i)-0.0,...
            'o','color',[235,190,190]/256,'MarkerFacecolor',...
            col,'MarkerSize',2);
    end
    
    %-- orientation arrow
    %-----------------------
    if para_.Publish_1tARMEx_OrientationArrow == 1
        
        arrow_length = 3;
        line_width_arrow = 1;%2; %0.25
        arrow_head_size = 30; %20;
        marker_size_arrow = 2; %1; %2
        h_reselct(n+4) = quiver(selectdConf_r(i)-0.0,selectdConf_c(i)-0.0,...
               arrow_length*cosd(reselectedConf2Exe_orn(i)),...
               arrow_length*sind(reselectedConf2Exe_orn(i)),...
               'LineWidth',line_width_arrow,...
               'color',col,...
               'MaxHeadSize',arrow_head_size,... %'Marker','o',...
               'MarkerSize',marker_size_arrow);
    end
    
    n = n+4;
end


file_name = sprintf('finalseminar/%s-1t-armex-gdm-plan',para_.ExperimentTitle);
export_fig([file_name], '-pdf','-r500')



if strcmp(para_.SensingSystem,'simulated-sampling')
    %pause
end

if numel(h_reselct) > 0
    delete(h_reselct);
end
pause(1)

end

%}





function fPlotMaps(map_env,dir_)
%
%

% -- environment
% -----------------------
figure; hold on;
% [map_row,map_col] = size(map_env);
% for i = 1:map_row
%     for j = 1:map_col
%         % if obstacle/free
%         if ~map_env(i,j) % occupied
%             col = [0.4 0.4 0.4];
%             faceCol = [0.4 0.4 0.4];
%         elseif map_env(i,j) % free
%             col = [0.9 0.9 0.9];
%             faceCol = [1 1 1];
%             %col = [1 1 1];
%         end
%         plot(i+0.5,j+0.5,'S','Color',col,...
%                          'MarkerFaceColor',faceCol,...
%                          'MarkerSize',20);
%     end
% end

imshow(map_env','InitialMagnification',1400);
hold on; set(gca,'YDir','normal'); 
colormap('gray'); caxis([-1 1.25])
%pause
% -- ground truth
% -------------------------
map_con = dlmread([dir_.gt,'map_con.dat']);
[map_row,map_col] = size(map_con);
gdm_colormap = (flipud(colormap('hot'))); % for main fig. 
max_concentration = max(map_con(:)); % main fig.
delta = max_concentration/(size(gdm_colormap,1)-1);    
for i = 1:map_row
    for j = 1:map_col
        % if positive concentration
        if map_con(i,j)>0
            linear_color_number = round( min(map_con(i,j),max_concentration)/delta)+1;
            col = gdm_colormap( linear_color_number, :);
            % ----------------------------------------------------
            %plot(i+0.5,j+0.5,'s','Color',col,'MarkerFaceColor',col,'MarkerSize',15);
            plot(i+0.0,j+0.0,'s','Color',col,'MarkerFaceColor',col,'MarkerSize',10);
        end
    end
end
set(gca,'xtick',1:1:map_row);
set(gca,'ytick',1:1:map_col);
grid on; axis equal;


end






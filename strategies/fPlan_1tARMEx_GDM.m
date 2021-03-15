function [ hc_num,...
           reselectedConf2Exe_ind,...
           reselectedConf2Exe_orn] = fPlan_1tARMEx_GDM( hc_num,...
                                                        hc_estimated_sub,...
                                                        hc_block_num,...
                                                        conf_sequence_num,...
                                                        conf_sequence_list_hc,...
                                                        confs_executed,...
                                                        map_env,...
                                                        recon_local_x,...
                                                        recon_local_y,...
                                                        map_gd_local,...
                                                        FoV,...
                                                        FoVfaces,...
                                                        confs_executed_hc,...
                                                        OBSenv,...
                                                        OBSconf,...
                                                        SensorRange,...
                                                        cellsize_env,...
                                                        dir_,...
                                                        para_ )
                                                       
%
%*************************************************************************
%  A tomographic plan for the estimated hcs.
%************************************************************************* 
% 
%



h_reselct_prev = [];


%/////////////////////////////////
% -- update current hotspot
%/////////////////////////////////
hc_current_sub = hc_estimated_sub;

% -- plot
% ----------------------------
% plot(hc_current_sub(:,1),hc_current_sub(:,2),...
%     'ob','MarkerSize',8,'MarkerFaceColor','g');
% pause

%/////////////////////////////////
% -- reference point for distance
%/////////////////////////////////            
refPointDistance = confs_executed(conf_sequence_num,1:2);


%/////////////////////////////////
% -- current gdm
%/////////////////////////////////
map_gd_current = zeros(size(map_env));
% map_gd_current(round(recon_local_x),round(recon_local_y)) = map_gd_local';
map_gd_current(round(recon_local_x),round(recon_local_y)) = map_gd_local;


%////////////////////////////////////////////////////
% 
%               LOCAL SOLUTIONS
% 
%////////////////////////////////////////////////////

fprintf(1,'\n   ****************************************************\n')
fprintf(1,'                     LOCAL SOLUTIONS         \n')
fprintf(1,'   ****************************************************\n\n')

[ hc_num,...
  hc_nums_local_this,...
  map_coverage,...
  refDistPoint] = fOptimalHotspotGeometries1t( hc_current_sub,...
                                               hc_num,...
                                               hc_block_num,...
                                               FoV,...
                                               map_env,...
                                               OBSenv,...
                                               OBSconf,...
                                               SensorRange,...
                                               map_gd_current,...
                                               confs_executed_hc,...
                                               para_,...
                                               dir_ );


%////////////////////////////////////////////////////
% 
%           FUSION OF THE LOCAL SOLUTIONS
% 
%////////////////////////////////////////////////////

fprintf(1,'\n   ****************************************************\n')
fprintf(1,'             FUSION OF THE LOCAL SOLUTIONS         \n')
fprintf(1,'   ****************************************************\n\n')

visualize = 0;


% hc_current_sub
% hc_nums_local_this
% hc_block_num

% preprocess parameters:
% --------------------------
% number of elementray sensing action per cell.
para_.numFoVs = 360; 
para_.SensorDomainComputation = 'UseEarlier'; % 'ComputeNew'
% circular area for the candidate conf around redundant confs
para_.candConfRadiusForRedConfs_m = 1.5*SensorRange.cell; 
% min radius where sensing conf are not allowed.
para_.innerRadius_m = 1; 
% number of allowed cand conf.
%para_.NumOfCandConf = 200; 
% min distance to check redundant confs.
para_.DistanceForRedundantConfs_m = 2*para_.SensingRangeM; %3; 

%para_.toleranceERQ_Fusion_percentage


[ reselectedConf2Exe_ind,...
  reselectedConf2Exe_orn,...
  tFusion ] = fFusionHotspotGeometries1t( FoV,...
                                          hc_current_sub,...
                                          hc_nums_local_this,...
                                          hc_block_num,...
                                          SensorRange,...
                                          FoVfaces,...
                                          map_env,...
                                          map_gd_current,...
                                          map_coverage,...
                                          confs_executed_hc,...
                                          refDistPoint,...
                                          OBSenv,...
                                          OBSconf,...
                                          cellsize_env,...
                                          dir_,...
                                          para_,...
                                          visualize );
                      

filename = sprintf('spp_solution_block%02d.mat',hc_block_num);
save([dir_.sppgdm,filename],...
    'reselectedConf2Exe_ind',...
    'reselectedConf2Exe_orn',...
    'tFusion');

                                
end

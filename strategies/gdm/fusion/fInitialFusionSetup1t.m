function [ localConfs,...
           localERGall,...
           localERQall,...
           localNumOfSelectedConfs] = fInitialFusionSetup1t( hcs_sub,...
                                                              hc_nums_local_this,...
                                                              map_env,...
                                                              dir_,...
                                                              para_,...
                                                              visualize )
%fFusionSetup identify the redundant pairs of selected conf from the local
%solutions.



%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
%                RETRIEVE INFO FROM THE LOCAL SOLUTIONS (SPP-HOTSPOTS)
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤

% Lets retrieve objective values (G+U, G, and U) of subproblems, reference
% distance points, maximum Ds and 1Vs.

% -- initialize.
localConfs.ind = cell(1,numel(hc_nums_local_this));
localConfs.orn = cell(1,numel(hc_nums_local_this));

localNumOfSelectedConfs = zeros(numel(hc_nums_local_this),1);

localERGall  = zeros(numel(hc_nums_local_this),1);
% localERGGall = zeros(numel(hc_nums_local_this),1);
% localERGUall = zeros(numel(hc_nums_local_this),1);
localERQall  = zeros(numel(hc_nums_local_this),1);

% refPointDist = zeros(numel(hc_nums_local_this),2); % row col of ref point

if visualize == 1
figure('name','redundant'); 
imshow(map_env','InitialMagnification',700); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
end

col_lines = lines(numel(hc_nums_local_this));

% for each hotspot
for i = 1:numel(hc_nums_local_this)
    
    
    %-- load preprocessing
    filename = sprintf('optimal_geometry_hc%02d_xvd.mat',hc_nums_local_this(i));
    this = load([dir_.sppgdm,filename],...
                  'wV',...
                  'Conf_DistOrnGainV',...
                  'confKept',...
                  'ConfOrns',...
                  'NumOfAllowedConf',...
                  'conf_crossAngles_G',...
                  'D');
    % load SPP solution for the center of mass          
    filename = sprintf('optimal_geometry_hc%02d_spp.dat',hc_nums_local_this(i));
    selectedConf = load([dir_.sppgdm,filename]);
    selectedConf = round(selectedConf);
    
    localConfs.ind{i} = this.confKept(selectedConf==1);
    localConfs.orn{i} = this.ConfOrns(selectedConf==1);
    
    localNumOfSelectedConfs(i) = this.NumOfAllowedConf;
    
    
    
    % debug
    if visualize == 1
        conf_debug_ind = this.confKept(selectedConf==1);
        confcell_debug_ind = conf_debug_ind; %fix(conf_debug_ind/f)+(~~mod(conf_debug_ind,f)); 

        [confcell_debug_r,confcell_debug_c] = ind2sub(size(map_env),confcell_debug_ind);
        %plot(confcell_debug_c,confcell_debug_r,'x','Color',col_lines(i,:));
        plot(confcell_debug_r,confcell_debug_c,'x','Color',col_lines(i,:));

        %plot(hot_debug_c,hot_debug_r,'o','Color',col_lines(i,:));
        plot(hcs_sub(:,1),hcs_sub(:,2),'o','Color',col_lines(i,:));
    end
    %}
    
    %-- current ERG (gain)
    %---------------------------------------------
    [ thisERG ] = fCalculateERG( this.NumOfAllowedConf,...
                                 selectedConf,...
                                 this.wV,...
                                 this.D,...
                                 this.conf_crossAngles_G,...
                                 this.Conf_DistOrnGainV,...
                                 para_ );
    
    %-- current ERQ (map quality)
    %---------------------------------------------
    nc = this.NumOfAllowedConf;
    if     nc == 2
        thisERQ = para_.MapQuality2Conf*thisERG;
    elseif nc == 3
        thisERQ = para_.MapQuality3Conf*thisERG;
    elseif nc == 4
        thisERQ = para_.MapQuality4Conf*thisERG;
    elseif nc == 5
        thisERQ = para_.MapQuality5Conf*thisERG;
    end
    
    %-- update the results
    %---------------------------------------------    
    %localERGall(i,1) = C'*(G*C+U);
    %localERGGall(i,1) = C'*(G*C);
    %localERGUall(i,1) = C'*U;
    
    localERGall(i,1) = thisERG;
    localERQall(i,1) = thisERQ;
    
    % reference distance point [r,c] in each hotspot
    %refPointDist(i,:) = confs_executed_hc(1,1:2);
    
    
    
end





end

%------------------------------- End of the document -------------------------------------


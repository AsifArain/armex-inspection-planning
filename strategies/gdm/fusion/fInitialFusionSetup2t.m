function [ localConfs,...
           localERGall,...
           localERQall,...
           localNumOfSelectedConfs,...
           refPointDist ] = fInitialFusionSetup2t( hcsSEQ_ind,...
                                                   map_env,...
                                                   dir_,...
                                                   para_,...
                                                   visualize )
%fFusionSetup identify the redundant pairs of selected conf from the local
%solutions.



%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
%                RETRIEVE INFO FROM THE LOCAL SOLUTIONS (SPP-HOTSPOTS)
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤

% Lets retrieve objective values (G+U, G, and U) of subproblems, reference
% distance points, maximum Ds and 1Vs.

% -- initialize.
% selectedConf_ind = zeros(size(hotspotSEQ,1),para_.NumberOfConf_SPPH);
% selectedConf_ont = zeros(size(hotspotSEQ,1),para_.NumberOfConf_SPPH);

localConfs.ind = cell(1,size(hcsSEQ_ind,1));
localConfs.orn = cell(1,size(hcsSEQ_ind,1));

localNumOfSelectedConfs = zeros(size(hcsSEQ_ind,1),1);

localERGall  = zeros(size(hcsSEQ_ind,1),1);
% localERGGall = zeros(size(hcsSEQ_ind,1),1);
% localERGUall = zeros(size(hcsSEQ_ind,1),1);
localERQall  = zeros(size(hcsSEQ_ind,1),1);


refPointDist = zeros(size(hcsSEQ_ind,1),2); % row col of ref point

if visualize == 1
figure('name','redundant'); 
imshow(map_env','InitialMagnification',700); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
end

col_lines = lines(numel(hcsSEQ_ind));

% for each hotspot
for hc_num = 1:numel(hcsSEQ_ind)
    
    % load preprocessing
    filename = sprintf('optimal_geometry_hc%02d_xvd.mat',hc_num);
    this = load([dir_.gdmplan_spp,filename],...
                  'wV',...
                  'Conf_DistOrnGainV',...
                  'confKept',...
                  'ConfOrns',...
                  'conf_crossAngles_G',...
                  'refDistPoint',...
                  'NumOfAllowedConf',...
                  'D');
    
    % load SPP solution for the center of mass
    filename = sprintf('optimal_geometry_hc%02d_spp.dat',hc_num);
    selectedConf = load([dir_.gdmplan_spp,filename]);
    selectedConf = round(selectedConf);

    localConfs.ind{hc_num} = this.confKept(selectedConf==1);
    localConfs.orn{hc_num} = this.ConfOrns(selectedConf==1);
    
    localNumOfSelectedConfs(hc_num) = this.NumOfAllowedConf;
    
    
    % debug
    if visualize == 1
    conf_debug_ind = this.confKept(selectedConf==1);
    confcell_debug_ind = conf_debug_ind; %fix(conf_debug_ind/f)+(~~mod(conf_debug_ind,f)); 
        
    [confcell_debug_r,confcell_debug_c] = ind2sub(size(map_env),confcell_debug_ind);
    %plot(confcell_debug_c,confcell_debug_r,'x','Color',col_lines(i,:));
    plot(confcell_debug_r,confcell_debug_c,'x','Color',col_lines(hc_num,:));
    
    hot_debug_ind = hcsSEQ_ind(hc_num);
    [hot_debug_r,hot_debug_c] = ind2sub(size(map_env),hot_debug_ind);
    %plot(hot_debug_c,hot_debug_r,'o','Color',col_lines(i,:));
    plot(hot_debug_r,hot_debug_c,'o','Color',col_lines(hc_num,:));
    end
    %}
    
    
    
    %/////////////////////////////////
    % 	SOLUTION QUALITY (ERQ)
    %/////////////////////////////////
    
    % -- current ERQ
    % ----------------------------    
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
       
    
    localERGall(hc_num,1) = thisERG;
    localERQall(hc_num,1) = thisERQ;
    
    
    % reference distance point [r,c] in each hotspot
    refPointDist(hc_num,:) = this.refDistPoint;
    
    
end





end

%------------------------------- End of the document -------------------------------------


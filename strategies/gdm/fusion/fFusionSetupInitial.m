function [ ConfFusionPairs,...
           selectedConf_ind,...
           selectedConf_ont,...
           ERQ___local,...
           ERQ_G_local,...
           ERQ_U_local,...
           refPointDist ] = fFusionSetupInitial( hotspotSEQ,...
                                                 map_env,...
                                                 alpha,...
                                                 beta,...
                                                 gamma,...
                                                 NumberOfConf_SPPH,...
                                                 dataYAML,...
                                                 dir_,...
                                                 para_ )
%fFusionSetup identify the redundant pairs of selected conf from the local
%solutions.



%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
%                RETRIEVE INFO FROM THE LOCAL SOLUTIONS (SPP-HOTSPOTS)
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤

% Lets retrieve objective values (G+U, G, and U) of subproblems, reference
% distance points, maximum Ds and 1Vs.

% -- initialize.
selectedConf_ind = zeros(size(hotspotSEQ,1),NumberOfConf_SPPH);
selectedConf_ont = zeros(size(hotspotSEQ,1),NumberOfConf_SPPH);

% -- selected conf index orientation and corresponding hotspots
selectedConf_ioh = zeros(size(hotspotSEQ,1)*NumberOfConf_SPPH,2+size(hotspotSEQ,1)); 

ERQ___local = zeros(size(hotspotSEQ,1),1);
ERQ_G_local = zeros(size(hotspotSEQ,1),1);
ERQ_U_local = zeros(size(hotspotSEQ,1),1);

refPointDist = zeros(size(hotspotSEQ,1),2); % row col of ref point


figure('name','redundant'); 
imshow(map_env','InitialMagnification',700); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])

col_lines = lines(numel(hotspotSEQ));

% for each hotspot
for i = 1:numel(hotspotSEQ)
    
    % load preprocessing
    this = load([dir_.results,'preprocess_XVD_hotspot',num2str(i),'.mat'],...
                  'oV',...
                  'V',...
                  'wV',...
                  'confKept',...
                  'ConfOrientations',...
                  'oConfOrientations',...
                  'conf_crossAngles_G',...
                  'refDistPoint',...
                  'Hot2Conf_Angles',...
                  'D');
    
    % load SPP solution for the center of mass
    selectedConf = load([dir_.results,'spp_selectedconf_hotspot',num2str(i),'.dat']);
    selectedConf = round(selectedConf);

    selectedConf_ind(i,:) = this.confKept(selectedConf==1);
    selectedConf_ont(i,:) = this.ConfOrientations(selectedConf==1);
    
    % -- update new selected conf for ind ont and hotspots
    selectedConf_ioh(((i-1)*NumberOfConf_SPPH)+(1:3),1) = selectedConf_ind(i,:);
    selectedConf_ioh(((i-1)*NumberOfConf_SPPH)+(1:3),2) = selectedConf_ont(i,:);
    selectedConf_ioh(((i-1)*NumberOfConf_SPPH)+(1:3),i+2) = 1;
    
    
    % debug
    %    
    conf_debug_ind = this.confKept(selectedConf==1);
    confcell_debug_ind = conf_debug_ind; %fix(conf_debug_ind/f)+(~~mod(conf_debug_ind,f)); 
        
    [confcell_debug_r,confcell_debug_c] = ind2sub(size(map_env),confcell_debug_ind);
    %plot(confcell_debug_c,confcell_debug_r,'x','Color',col_lines(i,:));
    plot(confcell_debug_r,confcell_debug_c,'x','Color',col_lines(i,:));
    
    hot_debug_ind = hotspotSEQ(i);
    [hot_debug_r,hot_debug_c] = ind2sub(size(map_env),hot_debug_ind);
    %plot(hot_debug_c,hot_debug_r,'o','Color',col_lines(i,:));
    plot(hot_debug_r,hot_debug_c,'o','Color',col_lines(i,:));
    
    %}
    
    
    % spp-hotspots solution/objective values
    C = selectedConf;
    G = alpha*this.conf_crossAngles_G;
    
    %Uc = gamma*((ones(size(priorH.V,1),1)'*priorH.V)/max((ones(size(priorH.V,1),1)'*priorH.V)));
    Uc = gamma*((ones(size(this.wV,1),1)'*this.wV)/max((ones(size(this.wV,1),1)'*this.wV)));
    Ud = (1-gamma)*(1-(this.D/max(this.D)));
    U  = beta*(Uc'+Ud);
    
    ERQ___local(i,1) = C'*(G*C+U);    
    ERQ_G_local(i,1) = C'*(G*C);
    ERQ_U_local(i,1) = C'*U;
    
    % max distances in Ds of each hotpsot
    %maxDs_sppH(i) = max(this.D);
    
    % max visibilites in 1Vs of each hotpsot
    %max1Vs_sppH(i) = max(sum(this.V,1));
    
    % reference distance point [r,c] in each hotspot
    refPointDist(i,:) = this.refDistPoint;
    
    
end


selectedConf_ioh

conf_comb = nchoosek(1:size(selectedConf_ioh,1),2)

% -- selected conf index, orientation, corresponding hotspots and dist

selectedConfComb_iohd = zeros(size(conf_comb,1),4+size(hotspotSEQ,1)+1)

selectedConf_ioh

for i = 1:size(conf_comb,1)
    
    selectedConfComb_iohd(i,1) = selectedConf_ioh(conf_comb(i,1),1);
    
    selectedConfComb_iohd(i,2) = selectedConf_ioh(conf_comb(i,1),2);
    
    selectedConfComb_iohd(i,3) = selectedConf_ioh(conf_comb(i,2),1);
    
    selectedConfComb_iohd(i,4) = selectedConf_ioh(conf_comb(i,2),2);
    
    selectedConfComb_iohd(i,4+(1:size(hotspotSEQ,1))) = logical( ... 
        selectedConf_ioh(conf_comb(i,1),2+(1:size(hotspotSEQ,1))) + ...
        selectedConf_ioh(conf_comb(i,2),2+(1:size(hotspotSEQ,1))) );
    
    
    [conf1_r,conf1_c] = ind2sub(size(map_env),selectedConfComb_iohd(i,1))
    [conf2_r,conf2_c] = ind2sub(size(map_env),selectedConfComb_iohd(i,3));
    
    
    selectedConfComb_iohd(i,4+(size(hotspotSEQ,1))+1) = ...
        pdist2( [conf1_r,conf1_c],[conf2_r,conf2_c],'euclidean' );
    
    
end

    
selectedConfComb_iohd

conf_hots = sum(selectedConfComb_iohd(:,4+(1:size(hotspotSEQ,1))),2)

selectedConfComb_iohd = selectedConfComb_iohd(conf_hots>1,:)

conf_dist = selectedConfComb_iohd(:,end)
DistForRedConfs_cell = para_.DistanceForRedundantConfs_m/dataYAML.resolution;

selectedConfComb_iohd = selectedConfComb_iohd(conf_dist<=DistForRedConfs_cell,:)


[~,ind] = sortrows(selectedConfComb_iohd,4+(size(hotspotSEQ,1))+1);

selectedConfComb_iohd = selectedConfComb_iohd(ind,:)

pause




%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
%                               CONF PAIRS FOR FUSION
%¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤

% ConfFusionPairs = [conf_ind1,conf_ind2,hotspot1,hotspot2,distance]
ConfFusionPairs = [];
for this_hot = 1:numel(hotspotSEQ)-1
    for this_conf = 1:NumberOfConf_SPPH
        
        for other_hot = this_hot+1:1:numel(hotspotSEQ)
            for other_conf = 1:NumberOfConf_SPPH
                
                ConfFusionPairs = [ ConfFusionPairs; ...
                                    selectedConf_ind(this_hot,this_conf),...
                                    selectedConf_ind(other_hot,other_conf),...
                                    this_hot,other_hot];
            end
        end
    end
end



% -- Distances between the conf pairs
ConfFusionPairs(:,5) = zeros(size(ConfFusionPairs,1),1);
[ConfPairs1_r,ConfPairs1_c] = ind2sub(size(map_env),ConfFusionPairs(:,1));
[ConfPairs2_r,ConfPairs2_c] = ind2sub(size(map_env),ConfFusionPairs(:,2));

for i = 1:size(ConfFusionPairs,1)
    
    ConfFusionPairs(i,5) = pdist2( [ConfPairs1_r(i),ConfPairs1_c(i)],...
                             [ConfPairs2_r(i),ConfPairs2_c(i)],...
                             'euclidean' );
end


% -- Conf paris with sorted min distance (in cells)
DistForRedConfs_cell = para_.DistanceForRedundantConfs_m/dataYAML.resolution;
ConfFusionPairs = ConfFusionPairs(ConfFusionPairs(:,5)<=DistForRedConfs_cell,:);
[~,ind] = sortrows(ConfFusionPairs,5);

ConfFusionPairs = ConfFusionPairs(ind,:);



end

%------------------------------- End of the document -------------------------------------
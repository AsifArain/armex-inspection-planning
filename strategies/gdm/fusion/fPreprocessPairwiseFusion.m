function [ V,...
           oV,...
           constV,...
           confKept,...
           confRemoved,...
           cellKept,...
           cellRemoved,...
           FoVfaces,...
           map_candConfFusion,...
           confTheta,...
           confCrossAngles,...
           confCrossAngles_G,...
           redundantConf_ind,...
           fixedConf_subH,...
           subHotspots,...
           infoGainSubH_sppH,...
           subH_RedundantConfPairsAll_sppH,...
           Z ] = fPreprocessPairwiseFusion( FoVfaces,...
                                            hotspotSEQ,...
                                            map_env,...
                                            FoV,...
                                            map_coverage,...
                                            SensorRange,...
                                            OBS,...
                                            alpha,...
                                            beta,...
                                            infoGain_sppH,...
                                            infoGainG_sppH,...
                                            infoGainU_sppH,...
                                            maxDs_sppH,...
                                            max1Vs_sppH,...
                                            selectedConf_sppH,...
                                            selectedConfGlobal_sppH,...
                                            redundantConfPairsAll_sppH,...
                                            subH_RedundantConfPairsAll_sppH,...
                                            referenceDistPoint_sppH,...
                                            PairNumToHandle,...
                                            dist_redundantConfPairsAll_sppH,...
                                            dataYAML,...
                                            dir_,...
                                            para_,...
                                            visualize )
%
% 

%% HOTSPOTS

hotspotSEQ;
Hotspots.ind = hotspotSEQ;
[rr,cc] = ind2sub(size(map_env),Hotspots.ind);
Hotspots.sub = [rr,cc];


%% SENSING PARAMETERS AND CONSTANTS:
% ---------------------------------------------------------------------------------------------------------------------------
FoVfaces.num        = para_.numFoVs;   % number of conf per cell.
n                   = numel(map_env);         % number of cells in the map.
f                   = FoVfaces.num;       % for simplification.
c                   = n*f;                % total conf.


%% UPDATE SELECTED CONFS BASED ON PREVIOUS FUSIONS.
% ---------------------------------------------------------------------------------------------------------------------------

selectedConf_sppHF.vec       = selectedConf_sppH.vec;
selectedConfGlobal_sppHF.vec = selectedConfGlobal_sppH.vec;
    

for i = 1:PairNumToHandle-1
    %i
    % load SPP solution for the center of mass
    selectedConf = load([dir_.TomographyFusion,'selectedConf_PairwiseFusion',...
        num2str(i),'.dat']);
    selectedConf = round(selectedConf);
    
    %numel(find(selectedConf))
    
    if numel(find(selectedConf))<=1
        
        % load preprocessing of first center of mass for initialization.
        priorF = load([dir_.TomographyFusion,'preprocess_pairwise_fusion',num2str(i),'.mat']);

        %subH = find(subH_RedundantConfPairsAll_sppH.vec(i,:));
        subH = find(priorF.subH_RedundantConfPairsAll_sppH.vec(i,:));
        %selectedConf_sppHF.vec(priorF.redundantConfPairs1stSPP_sppH.ind(i,:)) = 0;    
        
        %selectedConfGlobal_sppHF.vec(priorF.redundantConfPairs1stSPP_sppH.ind(i,:)) = 0;
        selectedConfGlobal_sppHF.vec(priorF.redundantConf_ind) = 0;
        for j = 1:numel(subH)
            selectedConf_sppHF.vec(priorF.redundantConf_ind,subH(j)) = 0;
        end

        % load SPP solution for the center of mass
        selectedConf = load([dir_.TomographyFusion,'selectedConf_PairwiseFusion',...
            num2str(i),'.dat']);
        selectedConf = round(selectedConf);

        % conf vector with all cells in the map.
        fusedConf_vec = zeros(size(priorF.oV,2),1);
        %priorF.confKept
        fusedConf_vec(priorF.confKept) = selectedConf;
        fusedConf_ind = find(fusedConf_vec);

        %selectedConf_sppHF.vec(fusedConf_ind,subH)  = 1;
        % chack for validity based on the coverage    
        for j = 1:numel(subH) 
            %j            
            valid_fusedConfSubH_vec = zeros(size(priorF.oV,2),1);
            %valid_fusedConfSubH_vec(priorF.confKept) = selectedConf;
            valid_fusedConfSubH_vec(priorF.confKept(find(selectedConf))) = priorF.constV(j,find(selectedConf));
            selectedConf_sppHF.vec(:,subH(j)) = selectedConf_sppHF.vec(:,subH(j)) + valid_fusedConfSubH_vec;
        end
        selectedConfGlobal_sppHF.vec(fusedConf_ind) = 1;
        
        % now, update redundant confs with fused conf all over.
        %fusedConf_ind
        for j = 1:2
            
            %find(redundantConfPairsAll_sppH.ind(:)== priorF.redundantConf_ind(j))
            
            %size(redundantConfPairsAll_sppH.ind)
            
            %priorF.redundantConf_ind
            redundantConfPairsAll_sppH.ind(find(redundantConfPairsAll_sppH.ind(:)== priorF.redundantConf_ind(j)));
            
            all_conf = redundantConfPairsAll_sppH.ind(find(redundantConfPairsAll_sppH.ind(:)== priorF.redundantConf_ind(j)));
            
            %confToReplace(j) = redundantConfPairsAll_sppH.ind(find(redundantConfPairsAll_sppH.ind(:)== priorF.redundantConf_ind(j)));
            confToReplace(j) = all_conf(1);
            %confToReplace2 = find(redundantConfPairsAll_sppH.ind(:)== priorF.redundantConf_ind(2))
        end
        
        % replace confs
        %redundantConfPairsAll_sppH.ind
        for j = 1:2
            
            % find pair numbers going to be replaced
            [pair_nums,~] = find(redundantConfPairsAll_sppH.ind==confToReplace(j));
            
            % replace current redundant confs appear elsewehere.
            redundantConfPairsAll_sppH.ind(find(redundantConfPairsAll_sppH.ind(:)==confToReplace(j))) = fusedConf_ind;
            
            % also update the subject hotspots.
            
            %subH_RedundantConfPairsAll_sppH.vec(pair_nums,:)
            
            %((subH_RedundantConfPairsAll_sppH.vec(i,:))'*ones(1,size( subH_RedundantConfPairsAll_sppH.vec(pair_nums,:),2)))'
            
            subH_RedundantConfPairsAll_sppH.vec(pair_nums,:) =...
                subH_RedundantConfPairsAll_sppH.vec(pair_nums,:) +...
                ((subH_RedundantConfPairsAll_sppH.vec(i,:))'*ones(1,size( subH_RedundantConfPairsAll_sppH.vec(pair_nums,:),1)))';
            %subH_RedundantConfPairsAll_sppH.vec
        end        
        
        
    end
end


selectedConfGlobal_sppHF.ind = find(selectedConfGlobal_sppHF.vec);
selectedConfGlobalCell_sppHF.ind = fix(selectedConfGlobal_sppHF.ind/f)+(~~mod(selectedConfGlobal_sppHF.ind,f));


%% INFORMATION GAIN FOR THE CONDENSED CONFIGURATIONS.
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% get G, Uc and Ud of condensed conf

% initialize


subHotspots.num = find((subH_RedundantConfPairsAll_sppH.vec(PairNumToHandle,:)));
subHotspots.ind = hotspotSEQ(subHotspots.num);
[rr,cc] = ind2sub(size(map_env),subHotspots.ind);
subHotspots.sub = [rr,cc];


infoGainG_condConf_sppH  = zeros(size(subHotspots.num));
infoGainUc_condConf_sppH = zeros(size(subHotspots.num));
infoGainUd_condConf_sppH = zeros(size(subHotspots.num));

fixedConf_subH.num = zeros(1,numel(subHotspots.num)); % number of fixed conf for each sub hotspot.
fixedConf_subH.vec = zeros(size(selectedConfGlobal_sppHF.vec,1),numel(subHotspots.num)); % vector of fixed conf for each sub hospot.

% for each center of mass
for i = 1:numel(subHotspots.num)
    
    % --- condensed configurations
    %condConf = find(condensedConf_sppH.vec(:,subHotspots_num(i)))
    
    % initialize condensed conf as selected conf.s
    %fixedConf = find(selectedConf_sppH.vec(:,subHotspots_num(i)))
    %fixedConf = selectedConf_sppHF.ind(:,subHotspots_num(i))    
    [fixedConf,~] = find(selectedConf_sppHF.vec(:,subHotspots.num(i)));
    
    % and then remove redundant conf
    %redundantConfPairsAll_sppH.ind(PairNumToHandle,:)
    fixedConf(find(ismember(fixedConf,redundantConfPairsAll_sppH.ind(PairNumToHandle,:)))) = [];
    
    
    % --- lets calculate cross angles fromt the scratch.
    fixedConf_cell_ind = fix(fixedConf/f)+(~~mod(fixedConf,f));
    [fixedConf_cell_r,fixedConf_cell_c] = ind2sub(size(map_env),fixedConf_cell_ind);
    
    % orientation
    %Hotspots.sub(subHotspots.num(i),:)
    deltaR = fixedConf_cell_r - Hotspots.sub(subHotspots.num(i),1);
    deltaC = fixedConf_cell_c - Hotspots.sub(subHotspots.num(i),2);
    
    fixConfTheta = atan2(deltaC,deltaR);
    
    % cross angles
    fixConfCrossAngles = zeros(size(fixedConf));
    % cross angles
    for j = 1:numel(fixedConf)
        for k = 1:numel(fixedConf)
            fixConfCrossAngles(j,k) = angleAbsDiff(fixConfTheta(j),fixConfTheta(k));
        end
    end
    

    % cross angles gain
    fixConfCrossAnglesG = zeros(size(fixConfCrossAngles));
    gain_sigma = deg2rad(para_.xGainSigma_DEG); % standard deviation    
    for gain_meu = deg2rad([para_.xGainMeu1_DEG,para_.xGainMeu2_DEG]) % expected value or center
        %gain_meu
        fixConfCrossAnglesG = fixConfCrossAnglesG +...
            gaussmf(fixConfCrossAngles,[gain_sigma,gain_meu]);
    end
    
    % visibility indes
    [ Vindex,Vvalidity ] = fCoverageIndex( fixedConf,...
                                          map_coverage,...
                                          map_env,...
                                          n,...
                                          f,...
                                          FoV,...
                                          FoVfaces,...
                                          SensorRange,...
                                          OBS,...
                                          Hotspots.ind(subHotspots.num(i)),...
                                          para_ );
    
    % Distance from the reference points
    D = zeros(1,numel(fixedConf));
    
    for j = 1:numel(fixedConf)
        manhattan_distance = abs(fixedConf_cell_r(j)-referenceDistPoint_sppH(subHotspots.num(i),1))+...
                             abs(fixedConf_cell_c(j)-referenceDistPoint_sppH(subHotspots.num(i),2));
        D(1,j) = manhattan_distance;
    end
      
    
    G  = alpha     * fixConfCrossAnglesG;
    Uc = beta      * ((Vindex)/max1Vs_sppH(subHotspots.num(i)));
    Ud = (1-beta)  * (1-(D/maxDs_sppH(subHotspots.num(i))));
    
    
    %C  = ones(numel(fixedConf),1);
    % info gains
    %infoGainG_condConf_sppH(i)  = C'*(G*C)
    %infoGainUc_condConf_sppH(i) = sum(Uc)
    %infoGainUd_condConf_sppH(i) = sum(Ud)
    
    C  = Vvalidity';
    % info gains
    infoGainG_condConf_sppH(i)  = C'*(G*C);
    infoGainUc_condConf_sppH(i) = C'*Uc';
    infoGainUd_condConf_sppH(i) = C'*Ud';
        
    % keep fixed conf for sub hotspots
    %fixedConf_subH(:,i) = fixedConf
    fixedConf_subH.num(1,i) = numel(fixedConf);
    fixedConf_subH.vec(fixedConf,i) = 1;
       
    
end




%% CANDIDATE CONF MAP (INITIAL)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

redundantConfCellsPairs1stSPP_sppH.ind =...
    fix(redundantConfPairsAll_sppH.ind/f)+(~~mod(redundantConfPairsAll_sppH.ind,f));

% redundantConfCellsPairs1stSPP_sppH.ind(PairNumToHandle,:)

% map for the candidate sensing conf
map_candConfFusion_1st = zeros(size(map_env));
map_candConfFusion_1st(redundantConfCellsPairs1stSPP_sppH.ind(PairNumToHandle,:)) = 1;

% plot
if visualize==1
pause(0.1)
figure('name','cand nearest conf');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
for i = 1:size(map_candConfFusion_1st,1)
    for j = 1:size(map_candConfFusion_1st,2)        
        if map_candConfFusion_1st(i,j)
            %plot(j,i,'*b','MarkerFaceColor','b');
            plot(i,j,'*b','MarkerFaceColor','b');
        end
    end
end
% plot(subHotspots.sub(:,2),subHotspots.sub(:,1),'or','MarkerFaceColor','r');
plot(subHotspots.sub(:,1),subHotspots.sub(:,2),'or','MarkerFaceColor','r');
end

% grow area
radiusFreeCell_cell = ceil(para_.candConfRadiusForRedConfs_m/dataYAML.resolution); %8

se = strel('disk',radiusFreeCell_cell,0); % outer limit
map_candConfFusion_1st = imdilate(map_candConfFusion_1st,se);

% remove cand conf on the occupied cells.
%map_candConfFusion_1st(OBS.ind(find(ismember(OBS.ind,find(map_candConfFusion_1st))))) = 0;

% ---- lets have new obs list by growing the obstacles one meter
thikness_cells = ceil(1.5/dataYAML.resolution);
se = strel('square',thikness_cells); % outer limit
map_grown_obs_list = ~map_env;
map_grown_obs_list = imdilate(map_grown_obs_list,se);
map_grown_obs_list = ~map_grown_obs_list;


if visualize==1
pause(0.1)
figure('name','map_grown_obs_list'); 
imshow((map_env-map_grown_obs_list)','InitialMagnification',800);
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
end

obs_list_grown = find(~map_grown_obs_list);
map_candConfFusion_1st(obs_list_grown(find(ismember(obs_list_grown,find(map_candConfFusion_1st))))) = 0;



% plot
if visualize==1
pause(0.1)
figure('name','cand nearest conf enlarged area');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
for i = 1:size(map_candConfFusion_1st,1)
    for j = 1:size(map_candConfFusion_1st,2)
        if map_candConfFusion_1st(i,j)
            %plot(j,i,'xb','MarkerFaceColor','b');
            plot(i,j,'xb','MarkerFaceColor','b');
        end
    end
end
% plot(subHotspots.sub(:,2),subHotspots.sub(:,1),'or','MarkerFaceColor','r');
plot(subHotspots.sub(:,1),subHotspots.sub(:,2),'or','MarkerFaceColor','r');
end


%% VISIBILITY MATRIX:
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

disp('V...')
[ V,...
  oV,...
  constV,...
  confKept,...
  confRemoved,...
  cellKept,...
  cellRemoved,...
  FoVfaces ] = fV( map_candConfFusion_1st,...
                   map_coverage,...
                   map_env,...
                   n,...
                   f,...
                   c,...
                   FoV,...
                   FoVfaces,...
                   SensorRange,...
                   OBS,...
                   subHotspots,...
                   para_ );



%% UPDATED CANDIDATE CONF MAP :
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% cand conf that do not cover hotspot are removed from the second map.

% initialize final candidate conf map
map_candConfFusion = zeros(size(map_candConfFusion_1st));

% corresponding cell numbers for conf kept.
confCELL = fix(confKept/f)+(~~mod(confKept,f)); % cell num for each conf num.

% update the map
map_candConfFusion(confCELL) = 1;
% figure; imshow(map_candConf)
% numel(find(map_candConf))*4

% plot on environment map
if visualize==1
pause(0.1)
figure('name','cand conf for fusion - final');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
% [candCell_r,candCell_c] = find(map_candConf);
for i = 1:size(map_candConfFusion,1)
    for j = 1:size(map_candConfFusion,2)        
        if map_candConfFusion(i,j)
            %plot(j,i,'xb','MarkerFaceColor','b');
            plot(i,j,'xb','MarkerFaceColor','b');
        end
    end
end
% plot(subHotspots.sub(:,2),subHotspots.sub(:,1),'or','MarkerFaceColor','r');
plot(subHotspots.sub(:,1),subHotspots.sub(:,2),'or','MarkerFaceColor','r');
% numel(find(map_candConfFusion))
end


%% CROSS ANGLES:
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

disp('Cross angles...')
[ confTheta,...
  confCrossAngles,...
  confCrossAngles_G] = fCrossAngles( f,...
                                     map_env,...
                                     confKept,...
                                     fixedConf_subH,...
                                     subHotspots,...
                                     para_,...
                                     visualize );

                           
%% DISTANCE MATRIX:
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
disp('Distance matrix...')

% ---- D (distance) matrix is manhattan distance between cand conf position
% and reference position
confKeptCells_ind = fix(confKept/f)+(~~mod(confKept,f));
[confKeptCells_r,confKeptCells_c] = ind2sub(size(map_env),confKeptCells_ind);

% subjectHotspots_num = find(subjectHotspots)

D = zeros(numel(subHotspots.ind),numel(confKeptCells_ind));

for i = 1:numel(subHotspots.ind)    
    for j = 1:numel(confKeptCells_ind)
        manhattan_distance = abs(confKeptCells_r(j)-referenceDistPoint_sppH(subHotspots.num(i),1))+...
                             abs(confKeptCells_c(j)-referenceDistPoint_sppH(subHotspots.num(i),2));        
        D(i,j) = manhattan_distance;
    end
end

                           
%% Z:
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% U = sum(V,1)
% U = beta*(U/max(U))
% 
% Z = alpha*confCrossAngles_G + [U;U]

% infoGain_sppH
% infoGainG_sppH
% infoGainU_sppH

infoGainSubH_sppH = infoGain_sppH(subHotspots.num);

% infoGainG_condConf_sppH
% infoGainUc_condConf_sppH
% infoGainUd_condConf_sppH
% infoGainU_sppH = infoGainUc_sppH + infoGainUd_condConf_sppH

% subHotspots.num
% infoGain_sppH(subHotspots.num)
% infoGainG_sppH(subHotspots.num)
% infoGainU_sppH(subHotspots.num)


redundantConf_ind = redundantConfPairsAll_sppH.ind(PairNumToHandle,:);

% max 1Vs to narmalize 1Vs
norm1Vs = 1./max1Vs_sppH(subHotspots.num);

% max Ds to narmalize Ds
normDs = 1./maxDs_sppH(subHotspots.num);

% Uc = gamma * (ones(numel(subjectHotspotsCELL_ind),1)*(sum(V,1)/max(sum(V,1))))


% Uc = gamma     * (norm1Vs*(sum(V,1))) + (infoGainUc_condConf_sppH'*ones(1,size(D,2)));
% 
% Ud = (1-gamma) * (1-((normDs*ones(1,size(D,2))).*D)) + infoGainUd_condConf_sppH'*ones(1,size(D,2));

% --- info gain is valid only if subject hotspot is covered by cand conf.
% so, lets do it!

Uc = beta     * ((norm1Vs*(sum(V,1))).*constV)                + (infoGainUc_condConf_sppH'*ones(1,size(D,2)));
Ud = (1-beta) * ((1-((normDs*ones(1,size(D,2))).*D)).*constV) + (infoGainUd_condConf_sppH'*ones(1,size(D,2)));

U = (1-alpha)*(Uc+Ud);

% G = alpha*(2*confCrossAngles_G) + infoGainG_condConf_sppH'*ones(1,size(D,2));
G = alpha*(2*(confCrossAngles_G.*constV)) + infoGainG_condConf_sppH'*ones(1,size(D,2));

Z = G+U;


% -----
% size(confCrossAngles_G)
% size((norm1Vs*(sum(V,1))))
% size((1-((normDs*ones(1,size(D,2))).*D)))
% size(constV)
% ----


for i = 1:numel(redundantConf_ind)
    C = ismember(confKept,redundantConf_ind(i));
    infoAll = Z*C';
    infoG   = G*C';
    infoU   = U*C';
    
end

%{
% check the solution - 1
C = ismember(confKept,redundantConf_ind(1));
Z*C'
G*C'
U*C'

% check the solution - 2
C = ismember(confKept,redundantConf_ind(2));
Z*C'
G*C'
U*C'

% check the solution - 3
C = ismember(confKept,redundantConf_ind(3));
Z*C'
% G*C'

% check the solution - 4
C = ismember(confKept,redundantConf_ind(4));
Z*C'
% G*C'

% info gain required for the spp-fusion
infoGain_sppHotspots = infoGain_sppH(subHotspots_num)

% check the solution - XX
C = zeros(size(confKept));
C(5) = 1;
Z*C'
% G*C'

%}


end

function [ confTheta,...
           confCrossAngles,...
           confCrossAngles_G ] = fCrossAngles( f,...
                                               map_env,...
                                               confKept,...
                                               fixedConf_subH,...
                                               subHotspots,...
                                               para_,...
                                               visualize )

% conf_theta            := conf orientations [1,c]
% conf_crossAngles      := cross angles between conf [c,c]
% conf_crossAngles_G    := cross angles gain [c,c]

%{
center_of_mass = paramPreprocess.CoM; % reference cell
[center_of_mass_y,center_of_mass_x] = ind2sub(size(map_candConf),center_of_mass);

% ----- conf orientations
conf_theta = 1000*ones(1,size(oV,2)); % vector contains orientations of candidate conf.

for conf = 1:size(oV,2)
    conf_cell_num = fix(conf/f)+(~~mod(conf,f));   % cell num for each conf num.
    [conf_cell_y,conf_cell_x] = ind2sub(size(map_candConf),conf_cell_num); % row and col of cell num.
    
    if map_candConf(conf_cell_y,conf_cell_x)
        deltaX = conf_cell_x - center_of_mass_x;
        deltaY = conf_cell_y - center_of_mass_y;
        conf_theta(conf) = atan2(deltaY,deltaX);
    end
    
end
%}

% conf orientations
% ---------------------------------------------------------------------------------
% Since, facing angle to a candidate conf is different for each hotspot.
disp('----- conf orientations -------')

confTheta = 1000*ones(numel(subHotspots.ind),numel(confKept)); % matrix contains orientations of candidate conf.

% plot for debugging
if visualize==1
pause(0.1)
figure('name','orientation angles debug');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
end

for i = 1:numel(subHotspots.ind)
    
    % index of current hotspot    
    if visualize==1    
    currentHotspot_H = plot(subHotspots.sub(i,1),subHotspots.sub(i,2),'or','MarkerFaceColor','r');
    end
    
    for j = 1:numel(confKept)
        
        confKeptCELL_ind = fix(confKept(j)/f)+(~~mod(confKept(j),f));   % cell num for each conf num.
        [confKeptCELL_r,confKeptCELL_c] = ind2sub(size(map_env),confKeptCELL_ind); % row and col of cell num.
        
        deltaR = confKeptCELL_r - subHotspots.sub(i,1); %TODO make sure subHotspotsCELL row and col are correct
        deltaC = confKeptCELL_c - subHotspots.sub(i,2);
        
        confTheta(i,j) = atan2(deltaC,deltaR);        
        
        % -- for debugging
        %{
        subConfAngle_degree = rad2deg(confTheta(i,j))
        subConfCell_ind = fix(confKept(j)/f)+(~~mod(confKept(j),f));
        [subConfCell_r,subConfCell_c] = ind2sub(size(map_env),subConfCell_ind);
        subConfCell_H = plot(subConfCell_c,subConfCell_r,'+b','MarkerFaceColor','b');        
        pause;
        delete(subConfCell_H);
        clc;
        %}
    end
    if visualize==1
    delete(currentHotspot_H);
    end
end

% confTheta

disp('........... end of orientation angle debug ---------------')


% cross angles
% ---------------------------------------------------------------------------------

disp('----- cross angles ------')


confCrossAngles = 1000*ones(sum(fixedConf_subH.num),numel(confKept));
% size(confCrossAngles)


% plot for debugging
if visualize==1
pause(0.1)
figure('name','cross angles debug');
imshow(map_env','InitialMagnification',800); hold on;
set(gca,'YDir','normal'); colormap('gray'); caxis([-2 2])
end

fixConfNum = 1;
for i = 1:numel(subHotspots.ind) % for subject hotspots
    
    % indices of fixed conf.s
    %currentFixedConfs_ind = fixedConf_subH(:,i); %find(condensedConfHotspots_vec(:,subHotspots_ind(i)))
    currentFixedConfs_ind = find(fixedConf_subH.vec(:,i)); %find(condensedConfHotspots_vec(:,subHotspots_ind(i)))
    
    %{ debug
    
    currentFixedConfsCell_ind = fix(currentFixedConfs_ind/f)+(~~mod(currentFixedConfs_ind,f));
    [currentFixedConfsCell_r,currentFixedConfsCell_c] = ind2sub(size(map_env),currentFixedConfsCell_ind);
    
    if visualize==1
    %plot(currentFixedConfsCell_c,currentFixedConfsCell_r,'xr')    
    plot(currentFixedConfsCell_r,currentFixedConfsCell_c,'xr')    
    %pause
    
    %}
    
    % index of current hotspot
    currentHotspot_ind = subHotspots.ind(i) %hotspotSEQ(subHotspots_ind(i))
    [currentHotspot_r,currentHotspot_c] = ind2sub(size(map_env),currentHotspot_ind);
    %currentHotspot_H = plot(currentHotspot_c,currentHotspot_r,'or','MarkerFaceColor','r');
    currentHotspot_H = plot(currentHotspot_r,currentHotspot_c,'or','MarkerFaceColor','r');
    end
    
    
    for j = 1:fixedConf_subH.num(1,i) 
        
        % orientation of fixed conf.
        %fixConfAngle = selectedConfTheta_sppH.vec(currentFixedConfs_ind(j),find(hotspotSEQ==subHotspots_ind(i)))
        %fixConfAngle = confTheta(i,find(ismember(confKept,currentFixedConfs_ind(j))))
        deltaR = currentFixedConfsCell_r(j) - subHotspots.sub(i,1);
        deltaC = currentFixedConfsCell_c(j) - subHotspots.sub(i,2);
        fixConfAngle = atan2(deltaC,deltaR);
         
        
        
        if visualize==1
        
        currentFixedConfCELL_ind = fix(currentFixedConfs_ind(j)/f)+(~~mod(currentFixedConfs_ind(j),f));
        [currentFixedConfCELL_r,currentFixedConfCELL_c] = ind2sub(size(map_env),currentFixedConfCELL_ind);
        
        %currentFixedConfCELL_H = plot(currentFixedConfCELL_c,currentFixedConfCELL_r,'ob');
        currentFixedConfCELL_H = plot(currentFixedConfCELL_r,currentFixedConfCELL_c,'ob');
        end
        
        for k = 1:numel(confKept)
            
            %fixConfNum
            %confKept(k)
            confCrossAngles(fixConfNum,k) = angleAbsDiff(fixConfAngle,confTheta(i,k));
            %confCrossAngles(((i-1)* numel(subjectHotspots_ind))+j,confKept(k)) = 
            
            % -- for debugging            
            %{
            fixConfAngle_degree = rad2deg(fixConfAngle)
            subConfAngle_degree = rad2deg(confTheta(i,k))
            cross_angles_degree = rad2deg(angleAbsDiff(fixConfAngle,confTheta(i,k)))
            
            subConfCell_ind = fix(confKept(k)/f)+(~~mod(confKept(k),f));
            [subConfCell_r,subConfCell_c] = ind2sub(size(map_env),subConfCell_ind);
            subConfCell_H = plot(subConfCell_c,subConfCell_r,'+b','MarkerFaceColor','b');
            
            pause;
            delete(subConfCell_H);
            clc;
            %}
            
        end
        
        fixConfNum = fixConfNum+1;
        if visualize==1
        delete(currentFixedConfCELL_H);
        end
    end
    if visualize==1
    delete(currentHotspot_H);
    end
end

% size(confCrossAngles)



% % conf variables corresponding to candidate conf map only: -- OLD --
% % ------------------------------------------------------------------------------
% 
% % --- instead, use info from fV:
% % update conf orientation vector
% confTheta = confTheta(confKept);
% 
% % update conf cross angles matrix
% confCrossAngles = confCrossAngles(:,confKept);
% % confCrossAngles = confCrossAngles(confKept,:);
% confCrossAngles = full(confCrossAngles);



% ------ cross angles gain
disp('----- gain')
gain_sigma = deg2rad(para_.xGainSigma_DEG); % standard deviation
gain_meu   = deg2rad(para_.xGainMeu1_DEG); % expected value or center
confCrossAngles_G1 = gaussmf(confCrossAngles,[gain_sigma,gain_meu]);


gain_sigma = deg2rad(para_.xGainSigma_DEG); % standard deviation
gain_meu   = deg2rad(para_.xGainMeu2_DEG); % expected value or center
confCrossAngles_G2 = gaussmf(confCrossAngles,[gain_sigma,gain_meu]);

confCrossAngles_G = confCrossAngles_G1+confCrossAngles_G2;
 
% confCrossAngles_G
% size(confCrossAngles_G)

I = zeros(numel(subHotspots.ind),size(confCrossAngles_G,1));
% fixedConf_subH

% for i = 1:size(fixedConf_subH,2)
%     i
%     ((i-1)*size(fixedConf_subH,1))+1:((i-1)*size(fixedConf_subH,1))+size(fixedConf_subH,1)
%     
%     I(i,((i-1)*size(fixedConf_subH,1))+1:((i-1)*size(fixedConf_subH,1))+size(fixedConf_subH,1)) = 1
% end

for i = 1:size(fixedConf_subH.num,2)
    if i == 1
        I(i,1:fixedConf_subH.num(i)) = 1;
    else
        I(i,sum(fixedConf_subH.num(1:i-1))+1:sum(fixedConf_subH.num(1:i-1))+fixedConf_subH.num(1,i) ) = 1;        
    end
end


confCrossAngles_G = I*confCrossAngles_G;

% size(confCrossAngles_G)



end


function [ Vindex,Vvalidity ] = fCoverageIndex( fixedConf,...
                                                  map_coverage,...
                                                  map_env,...
                                                  n,...
                                                  f,...
                                                  FoV,...
                                                  FoVfaces,...
                                                  SensorRange,...
                                                  OBS,...
                                                  subHot,...
                                                  para_ )
%fV generates Visibility matrix V.

% --- ANGLES:
FoVfaces.lgt = [1/FoVfaces.num:1/FoVfaces.num:1-(1/FoVfaces.num) 0]'*360;
FoVfaces.alf = FoVfaces.lgt-(FoV/2);
FoVfaces.bay = FoVfaces.lgt+(FoV/2);
% compensate negative alf
FoVfaces.alf(FoVfaces.alf<0) = 360+FoVfaces.alf(FoVfaces.alf<0);
FoVfaces.ang = [FoVfaces.alf FoVfaces.bay];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SensorDomain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('.. Sensor Domain')

switch para_.SenDom
    
    case 'ComputeNew'
        
        % For each cell of the map.
        sD = cell(1,f); % sensor domain (visible).
        oD = cell(1,f); % sensor's osbtacle domain (used to check the 
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
                sD{j} = visibleRange;
                
                
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
                oD{j} = obsRange;
                
            end
        end

        %--------------------------------------
        % EFFECTED VISIBILITY BY EACH OBSTACLE
        %--------------------------------------
        
        disp('.. Visibility vs Obstacles')
        VObs = cell(1,f); % Visibility vs Obstacle
        
        for i = 1
            sensor = [0,0];
            % For each configuration of a cell.
            for j = 1:f
                
                % let take sensor's mask
                sensorDom  = sD{j};
                % and obstacle domain
                obstacDom  = oD{j};
                
                VObs{j} = zeros(size(obstacDom,1),size(sensorDom,1)); % no more a square matrix
                % rows    := obstacle cell,
                % columns := visibility effect.
                
                %for k = 1:size(sensorDom,1)
                for k = 1:size(obstacDom,1)
                    
                    % kth cell is an obstacle
                    %obstacDomOBS(k,3)  = 0; 
                    obs = obstacDom(k,:);
                    
                    % visibility effect on all the cells
                    [ vDom ] = fVisSDvObs( sensor,sensorDom,obs );
                    
                    % update visibility against kth obstacle
                    VObs{j}(k,:) = vDom;
                   
                end
               
            end
        end
        
        % save the mask and visibility info to a file.
        save(['SenDom_Rng',num2str(SensorRange.cell),'FoV',num2str(FoV)],...
            'sD','oD','VObs');
        
    case 'UseEarlier'
        
        load(['SenDom_Rng',num2str(SensorRange.cell),'FoV',num2str(FoV)]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('.. Visibility Matrix')

% fixedConf

fixedConf_cell_ind = fix(fixedConf/f)+(~~mod(fixedConf,f));
[fixedConf_cell_r,fixedConf_cell_c] = ind2sub(size(map_env),fixedConf_cell_ind);

% V = sparse(n,c); % initializing
V = zeros(n,numel(fixedConf)); % initializing

% For each cell of the map.
for i = 1:numel(fixedConf)
    
    sensor = [fixedConf_cell_r(i),fixedConf_cell_c(i)];
    j = mod(fixedConf(i),f)+(~mod(fixedConf(i),f)*f);   

    % sensorDom := sensor domain translation w.r.t sensor position.
    sensorDom = [sD{j}(:,1)+sensor(1) sD{j}(:,2)+sensor(2)];

    % obstacDom := obstalce domain translated w.r.t sensor's
    % position.
    obstacDom = [oD{j}(:,1)+sensor(1) oD{j}(:,2)+sensor(2)];

    % VOBS := Visibility vs Obstacles matrix for jth FoV.
    VOBS = VObs{j};

    % visibile cells 
    [ visDom ] = fVisibileDomain( map_coverage,map_env,sensorDom,obstacDom,VOBS );

    % Update Visibility matrix.
    V(visDom,i) = 1;  
    
end

% sensing config can not visualize own cell
for i = 1:numel(fixedConf)
    V(fixedConf_cell_ind(i),i) = 0;
end

oV = V; % oV := visibility with obstacles

% VISIBILITY MATRIX CORRESPOND TO THE HOTSPOT ONLY:
% ------------------------------------------------------------------------------

% ---- cells:

% all cells covered by candidate conf.
coveredCell_ind = sum(oV,2);
coveredCell_num = find(coveredCell_ind);
coveredCell_vec = zeros(1,size(oV,1));
coveredCell_vec(coveredCell_num) = 1;

% cells need to be covered -- according to the coverage map
needToCoveCell_ind = find(map_coverage);
needToCoveCell_vec = zeros(1,size(oV,1));
% size(coverage_cell_vector)
needToCoveCell_vec(needToCoveCell_ind) = 1;
% size(coverage_cell_vector)

% unoccupied cells
unoccupiedCell_vec = ones(1,size(oV,1));
unoccupiedCell_vec(OBS.ind) = 0;

cellToKeep_vec = coveredCell_vec.*needToCoveCell_vec.*unoccupiedCell_vec;
cellKept       = find(cellToKeep_vec);
cellRemoved    = find(cellToKeep_vec==0);
V              = oV(cellKept,:);
Vvalidity      = oV(subHot,:);
Vindex         = sum(V,1);



end



function [ V,...
           oV,...
           constV,...
           confKept,...
           confRemoved,...
           cellKept,...
           cellRemoved,...
           FoVfaces ] = fV( map_candConfFusion_1st,...
                            map_coverage,...
                            map_env,...
                            n,...
                            f,...
                            c,...
                            FoV,...
                            FoVfaces,...
                            SensorRange,...
                            OBS,...
                            subHotspots,...
                            para_ )
%fV generates Visibility matrix V.

% --- ANGLES:
FoVfaces.lgt = [1/FoVfaces.num:1/FoVfaces.num:1-(1/FoVfaces.num) 0]'*360;
FoVfaces.alf = FoVfaces.lgt-(FoV/2);
FoVfaces.bay = FoVfaces.lgt+(FoV/2);
% compensate negative alf
FoVfaces.alf(FoVfaces.alf<0) = 360+FoVfaces.alf(FoVfaces.alf<0);
FoVfaces.ang = [FoVfaces.alf FoVfaces.bay];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SensorDomain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('.. Sensor Domain')

switch para_.SenDom
    
    case 'ComputeNew'
        
        % For each cell of the map.
        sD = cell(1,f); % sensor domain (visible).
        oD = cell(1,f); % sensor's osbtacle domain (used to check the 
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
                sD{j} = visibleRange;
                
                
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
                oD{j} = obsRange;
                
            end
        end

        %--------------------------------------
        % EFFECTED VISIBILITY BY EACH OBSTACLE
        %--------------------------------------
        
        disp('.. Visibility vs Obstacles')
        VObs = cell(1,f); % Visibility vs Obstacle
        
        for i = 1
            sensor = [0,0];
            % For each configuration of a cell.
            for j = 1:f
                
                % let take sensor's mask
                sensorDom  = sD{j};
                % and obstacle domain
                obstacDom  = oD{j};
                
                VObs{j} = zeros(size(obstacDom,1),size(sensorDom,1)); % no more a square matrix
                % rows    := obstacle cell,
                % columns := visibility effect.
                
                %for k = 1:size(sensorDom,1)
                for k = 1:size(obstacDom,1)
                    
                    % kth cell is an obstacle
                    %obstacDomOBS(k,3)  = 0; 
                    obs = obstacDom(k,:);
                    
                    % visibility effect on all the cells
                    [ vDom ] = fVisSDvObs( sensor,sensorDom,obs );
                    
                    % update visibility against kth obstacle
                    VObs{j}(k,:) = vDom;
                   
                end
               
            end
        end
        
        % save the mask and visibility info to a file.
        save(['SenDom_Rng',num2str(SensorRange.cell),'FoV',num2str(FoV)],...
            'sD','oD','VObs');
        
    case 'UseEarlier'
        
        %load(['SenDom_Rng',num2str(SensorRange.cell),'FoV',num2str(FoV)]);
        load(['preprocess/SenDom_Rng',num2str(SensorRange.cell),'FoV',num2str(FoV)]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('.. Visibility Matrix')
% n
% c

V = sparse(n,c); % initializing
% V = zeros(n,c); % initializing

% For each cell of the map.
for i = 1:n
    if map_candConfFusion_1st(i) % if its not occupied
        [rr,cc] = ind2sub(size(map_candConfFusion_1st),i);
        sensor = [rr,cc];
        
        % For each configuration of a cell.
        for j = 1:f
            
            % sensorDom := sensor domain translation w.r.t sensor position.
            sensorDom = [sD{j}(:,1)+sensor(1) sD{j}(:,2)+sensor(2)];
            
            % obstacDom := obstalce domain translated w.r.t sensor's
            % position.
            obstacDom = [oD{j}(:,1)+sensor(1) oD{j}(:,2)+sensor(2)];
            
            % VOBS := Visibility vs Obstacles matrix for jth FoV.
            VOBS = VObs{j};
            
            % visibile cells 
            [ visDom ] = fVisibileDomain( map_coverage,map_env,sensorDom,obstacDom,VOBS );
            
            % Update Visibility matrix.
            V(visDom,((i-1)*f)+j) = 1;
        end
    end
end



% sensing config can not visualize own cell
for i = 1:numel(map_candConfFusion_1st)    
    V(i,(i-1)*f+1:(i-1)*f+f) = 0;
end

oV = V; % oV := visibility with obstacles

% V = sparse(V);

%{
% OBSTACLE FREE VISIBILITY:
% ------------------------------------------------------------------------------
% -- removing cells
V(OBS.ind,:) = [];

% -- removing conf.
confOBS = zeros(numel(OBS.ind)*f,1);
for i = 1:numel(OBS.ind)
    % removing sensing configurations.
    %V(:,((OBS.ind(i)-1)*f+1)-(f*(i-1)):((OBS.ind(i)-1)*f+f)-(f*(i-1))) = [];
    confOBS(((i-1)*f)+1:((i-1)*f)+f) = (((OBS.ind(i)-1)*f)+1):(((OBS.ind(i)-1)*f)+f);
end
V(:,confOBS) = [];
%}

% VISIBILITY MATRIX CORRESPOND TO THE HOTSPOT ONLY:
% ------------------------------------------------------------------------------

% ---- cells:

% all cells covered by candidate conf.
coveredCell_ind = sum(oV,2);
coveredCell_num = find(coveredCell_ind);
coveredCell_vec = zeros(1,size(oV,1));
coveredCell_vec(coveredCell_num) = 1;

% cells need to be covered -- according to the coverage map
needToCoveCell_ind = find(map_coverage);
needToCoveCell_vec = zeros(1,size(oV,1));
% size(coverage_cell_vector)
needToCoveCell_vec(needToCoveCell_ind) = 1;
% size(coverage_cell_vector)

% unoccupied cells
unoccupiedCell_vec = ones(1,size(oV,1));
unoccupiedCell_vec(OBS.ind) = 0;

% size(coverage_cell_vector)

% pause
% --- selecting cells that are covered and need to be covered and unoccupied.
% size(covered_cell_vector)
% size(coverage_cell_vector)
% size(free_cell_vector)
cellToKeep_vec = coveredCell_vec.*needToCoveCell_vec.*unoccupiedCell_vec;
cellKept       = find(cellToKeep_vec);
cellRemoved    = find(cellToKeep_vec==0);
V              = oV(cellKept,:);
constV         = oV(subHotspots.ind,:);

% ---- configurations:

% --- selecting conf number corresponding to candidate conf map only
% cells for the candidate conf.
candConfCells = find(map_candConfFusion_1st); 
% conf numbers
confCand_vec = zeros(1,size(oV,2));
for i = 1:numel(candConfCells)
    confCand_vec( (candConfCells(i)-1)*f+1:(candConfCells(i)-1)*f+f ) = 1;
end

% --- conf that cover hotspot ---------------
% tmpV = oV;

% subjectHotspots_ind = find(subjectHotspots)

% numel(subjectHotspots_ind)
% 
% hotspotSEQ(subjectHotspots_ind)
% 
% subHotspots_ind


% tmpV = zeros(numel(subjectHotspots_ind),size(oV,2));
% for i = 1:numel(subjectHotspots_ind)
%     tmpV(i,:) = oV(hotspotSEQ(subjectHotspots_ind(i)),:);
% end
% sumV_hotspots = ones(1,numel(subjectHotspots_ind))*tmpV;


sumV_hotspots = ones(1,numel(subHotspots.ind))*oV(subHotspots.ind,:);

% conf_cover_all_hotspots_num = find(sumV_hotspots==numel(subjectHotspots_ind));
ConfCoverAllHotspot_ind = find(sumV_hotspots>=numel(subHotspots.ind));
ConfCoverAnyHotspot_ind = find(sumV_hotspots>=1);

confCoverHotspot_vec = zeros(1,size(oV,2));
% confCoverHotspot_vec(ConfCoverAllHotspot_ind) = 1;
confCoverHotspot_vec(ConfCoverAnyHotspot_ind) = 1;
% ---------------------------------------------

% cand conf belong to unoccupied cells
confUnoccupiedCell_vec = ones(1,size(oV,2));
for i = 1:numel(OBS.ind)
    confUnoccupiedCell_vec( (OBS.ind(i)-1)*f+1:(OBS.ind(i)-1)*f+f ) = 0;
end

% --- keeping candidate conf in the list that are from candidate conf map
% AND cover the hotspot AND are placed over unoccupied cells 

% size(confCand_vec)
% size(confCoverHotspot_vec)
% size(confUnoccupiedCell_vec)

validConf_vec = confCand_vec.*confCoverHotspot_vec.*confUnoccupiedCell_vec;
confKept      = find(validConf_vec);
confRemoved   = find(validConf_vec==0);
V             = V(:,confKept);
V             = full(V);
constV        = constV(:,confKept); % contraint V: hotspots cells times conf kept
constV        = full(constV);


end

function [ visibleRange ] = fSenDomain( sensor,SensorRange,zawia )
%fSenDomain returns visible cells (domain/mask) for a given sensing
%configuration. 
%

% INTERPOLATING INTERMEDIATE ANGLES:
% -------------------------------------
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

% initialize the obstacle range.
obsRange = []; 

% find the list of cells (for obstacles) in the sensing range.
for i = 1:numel(angl)
    
    % angle.
    ALF = angl(i);
    
    % Finding cells that are intersected along the ray.
    [ obsCells ] = fRayTracOD(sensor,SensorRange,ALF);
    
    % Updating the list of cells to be checked for obstacles.
    obsRange = unique([obsRange; obsCells],'rows');

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

function [ vDom ] = fVisSDvObs( sensor,sensorDom,obs )
%fVisibileArea returns visible cells for a given sensing configuration.
%

% initialize the list of visible cells.
vDom = zeros(size(sensorDom,1),1);

% for all cells in the sensor domain, check visibility against obs.
for i = 1:size(sensorDom,1)
    %i
    % cell to be validated for visibility
    subCell = [sensorDom(i,1) sensorDom(i,2)];

    % validatation for the visibility
    %[ cellValid ] = fCellValidSD( sensor,subCell,obstacDomOBS );
    [ validation ] = fCellValidSD( sensor,subCell,obs );

    % update visDom
    vDom(i) = validation;
    %pause
end

end


function [ validation ] = fCellValidSD( sensor,subCell,obs )
% 
% 

% angle to the center of the cell
ang = atan2(((subCell(2))-(sensor(2))),((subCell(1))-(sensor(1))))*(180/pi);
ang = ang+(360*(ang<0));

% conditions to project the line along (+/-x,+/-y).
cond_x1 = ang>=0&&ang<=45;
cond_x2 = ang>=135&&ang<=225;
cond_x3 = ang>=315&&ang<=360;

cond_y1 = ang>=45&&ang<=135;
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
cell_xy = unique(cell_xy,'rows');

% validation for the visibility
ind  = find(cell_xy(:,1)==obs(1)&cell_xy(:,2)==obs(2));
% if size(ind,1) > 0
%     validation = 0;
% end
validation = (size(ind,1) < 1);

end

function [ visDom ] = fVisibileDomain( map_coverage,map_env,sensorDom,obstacDom,VOBS )
%fVisibileDomain returns visible cells for a given sensing configuration.
%It takes the mask (sensor domain), finds the non-visible cells due to
%obstacles in the map and remove them from the mask.
%
% INPUTS:
% ------
% map       := [x,y] map
% sensorDom := [n,2] mask for the sensing configuration in question.
% obstacDom := [m,2] mask for the sensing configuration in question (used
%                    to check visibility against obstacles).
% VOBS      := [m,n] effect on visibility (columns) against each obstacle 
%                    (rows)
% OUTPUTS:
% ------
% visDom    := [k] visibile cells for the configuration in question 
%              (subscripts)

% rows and columns for the map
[map_row,map_col] = size(map_coverage);

% temporary domain/masks to handle intermediate computation quickly.
visDomX = sensorDom;
obsDomX = obstacDom;

% preRemoveList:= cells outside the map (sensor domain).
% along x-axis:
   i = find(visDomX(:,1)<=0);      % less than row number
   j = find(visDomX(:,1)>map_row); % greater than row number
   preRemoveList = [i;j];
% along y-axis:
   i = find(visDomX(:,2)<=0);      % less than col number
   j = find(visDomX(:,2)>map_col); % greater than col number
   preRemoveList = [preRemoveList;i;j];

% remove preRemoveList from temporary visible domain/mask.
% visDomX(preRemoveList,:) = [];

% preRemoveListOD:= cells outside the map (obstacle domain).
% along x-axis:
   i = find(obsDomX(:,1)<=0);      % less than row number
   j = find(obsDomX(:,1)>map_row); % greater than row number
   preRemoveListOD = [i;j];
% along y-axis:
   i = find(obsDomX(:,2)<=0);      % less than col number
   j = find(obsDomX(:,2)>map_col); % greater than col number
   preRemoveListOD = [preRemoveListOD;i;j];

% remove preRemoveList from temporary visible domain/mask.
obsDomX(preRemoveListOD,:) = [];


% convert visDomX from ind to sub.
% visDomX_ind = sub2ind(size(map),visDomX(:,1),visDomX(:,2));
obsDomX_ind = sub2ind(size(map_coverage),obsDomX(:,1),obsDomX(:,2));

% find the list of obstacles in visDomX;
% obs_list = find(map(visDomX_ind)==0);
obs_listOD = find(map_env(obsDomX_ind)==0);

% VOBSx:= VOBS with outside cells are removed.
% VOBSx = VOBS; 
% VOBSx(preRemoveList,:) = [];

VOBSx = VOBS; 
VOBSx(preRemoveListOD,:) = [];

% non-visible cells due to obstacles (logical values).
% nVisCells = (VOBSx(obs_list',:)==0);
nVisCells = (VOBSx(obs_listOD',:)==0);

% convert logical values to a list.
list = find(nVisCells'==1);

% translate the list to the column numbers of VOBS, which is remove list.
[removeList,~] = ind2sub(size(VOBSx'),list);

% update the RemoveList with preRemoveList and removeList.
RemoveList = [preRemoveList;removeList];

% finally, we know non-visible cells in the domain/mask.
visDom = sensorDom;
visDom(RemoveList,:) = [];

% convert indices into subscripts and sort them.
if visDom  
    visDom = sort(sub2ind(size(map_coverage),visDom(:,1),visDom(:,2)));
end

end

% ---------------------------------------------------- End of the document -------------------------------------------------------
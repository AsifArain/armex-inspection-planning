function [ SensorRange,...
           OBSenv,...
           OBSconf,...
           V,...
           oV,...
           FoVfaces,...
           confKept,...
           confRemoved,...
           cellKept,...
           cellRemoved,...
           T,...
           E,...
           eE,...
           A,...
           oA,...
           Ct,...
           o,...
           P,...
           tPreprocess ] = fPreprocessGD( FoV,...
                                          SensorRange,...
                                          para_,...
                                          dir_,...
                                          visualize )
% fPreProcess_R4 genereates required parameters for solving a global coverage problem.
% Date: 2014-01-19, Rev 4
% 
% n  := total number of cells in the map.
% nf := number of free cells in the map.
% f  := total number of sensing configurations per cell.
% 
% OUTPUTS:
% -----------------------------------------------------------------------------------
% cnt: [scalar][ ]   4 or 8 connectivity graph.
% o  : [scalar][ ]   Number of outgoing edges for each conf.
% E  : [nf*f*o][ ]   Edge connectivity. Index number is ID of outgoing edge
%                    and array value is ID of connected incoming edge.%                   
% vE : [------][ ]   E (as above) with reduced variables using Corner Reduction Method.
% eE : [nf*f*o][ ] 
% evE: [------][ ]   eE (as above) with reduced variables using Corner Reduction Method.
% T  : [nf*f,o][sec] Travelling cost. Rows are conf number and columns are outgoing edges.
% vT : [------][sec] T (as above) with reduced variables using Corner Reduction Method.
% V  : [nf,nf*f][binary] Visibility matrix, indicates visiblity status (1
%                        for visible, 0 otherwise) of nf cells (rows) from each conf (columns).
% vV : [-------][binary] V (as above) with reduced variables using Corner Reduction Method.
% vkS: [][]          .cell := selected special vertex cell number.
%                    .conf := configurations that cover selected special vertex cell.
% CRL: [][]          Corner Reduction List, list of configuration selected
%                    to be removed using Corner Reduction technique. 
% SensorRange: [scalar][integer] .cell sensor range in number of cells.
% map: [l,m][binary] Grid map matrix, 0 for occupied cells and 1 for
%                    unoccupied cells. l and m are any positive integer numbers. 
% FoVfaces [][]: .alf first/start angle for each FoV orientation.
%                .bay second/final angle for each FoV orientation.
%                .ang set of first/start and second/final angles of each FoV orientation.
%                .num total number of FoV faces/orientations for a sensing position.
% OBS: [][]      .ind Indices of occupied cells.
%                .sub Subscripts (row,col) of occupied cells.
% tPreprocess: [scalar][sec] Total computation time.
% 
% INPUTS:
% -----------------------------------------------------------------------------------
% FoV          : [][integer] Field of view.
% SensorRange  : [scalar][integer] .m sensor range in meters.
% paramPreprocess     : [][] Solution parameters.
%                     .CRL      := compute list of corner reduction variables (yes/no).
%                     .vkS      := selection of special vertex (yes/no).
%                     .Ct       := compute travelling cost matrix (yes/no).
%                     .GraphCnt := graph connectivity (4/8).
%                     .Res      := required cell size in meter.
%                     .TravCost := travelling cost (sec) for one meter travel.
%                     .PixSize  := original (map) pixel size (meters).
%                     .numFoVs  := number of elementray sensing action per cell.
%                     .A        := generate connectivity matrix A (yes/no).
%                     .AWM      := artificial world map (yes/no).
%                     .RWM      := real world map (yes/no).
%                     .MapFig   := map file e.g. 'frieburg.png'.
%                     .MapTxt   := map file e.g. 'kjula.txt'.
%                     .SenDom   := computing sensor domain. 'ComputeNew' to
%                                  compute from scratch, and 'UseEarlier' to use 
%                                  earlier computed sensor domain parameters.
% map          : [l,m][binary] Grid map matrix, 0 for occupied cells and 1
%                              for unoccupied cells. l and m are any positive 
%                              integer numbers. 
%
% NOTES/UPDATES:
% -----------------------------------------------------------------------------------
% R4 : Visibility Method VIS4 (2013-12-13) is used.
% R3 : Visibility Method VIS3 (2013-12-04) is used.
% R2 : ... 
% R1 : ...
% R0 : ...

tPreprocessi = tic; % initialize pre computational time.

% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             ENVIRONMENT MAP:
% 
% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
[ map_env,...
  map_conf,...
  origin_env,...
  origin_conf,...
  cellsize_env,...
  cellsize_conf,...
  cellsize_org ] = fMapStuff( para_,dir_ );



% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             SENSING PARAMETERS AND CONSTANTS:
% 
% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
FoVfaces.num = para_.numFoVs;   % number of conf per cell.
n            = numel(map_env);  % number of cells in the map_env.
f            = FoVfaces.num;    % for simplification.
c            = n*f;             % total conf.
cnt          = para_.GraphCnt;  % 4 or 8 connected neighbours.

% Sensor range in number of cells.
SensorRange.cell = (SensorRange.m/cellsize_env); % sensor range in cells



% o := number of outgoing edges for each configuration.
if cnt == 4
    if f == 1
        o = 4;
        oO = [1 1 1 1;...    % T - 1st 1
              0 0 0 0;...    % T - 2nd 1.4
              0 0 0 0;...    % T - 3rd inf
              0 0 0 0;...    % R - 1st 0.1
              0 0 0 0];      % R - 2nd 0.2          
         eO = oO;
          
    elseif f == 4
        o = 3;
        oO = [1 0 0;...    % T - 1st 1
              0 0 0;...    % T - 2nd 1.4
              0 0 0;...    % T - 3rd inf
              0 0 0;...    % R - 1st 0.1
              0 1 1];      % R - 2nd 0.2
        eO = oO;        
    end   
    
elseif cnt == 8
    if f == 1
        o = 8;
        %     1 2 3 4 5 6 7 8 -     IDs of outgoing edges
        oO = [1 0 1 0 1 0 1 0;...    % T - 1st 1
              0 1 0 1 0 1 0 1;...    % T - 2nd 1.4
              0 0 0 0 0 0 0 0;...    % T - 3rd inf
              0 0 0 0 0 0 0 0;...    % R - 1st 0.1
              0 0 0 0 0 0 0 0];      % R - 2nd 0.2
        eO = oO;
        
    elseif f == 8
        
        o = 3;
        oO = [1 0 0;...    % T - 1st 1
              0 0 0;...    % T - 2nd 1.4
              0 0 0;...    % T - 3rd inf
              0 1 1;...    % R - 1st 0.1
              0 0 0];      % R - 2nd 0.2
          
        eO = [0 0 0;...    % T - 1st 1
              1 0 0;...    % T - 2nd 1.4
              0 0 0;...    % T - 3rd inf
              0 1 1;...    % R - 1st 0.1
              0 0 0];      % R - 2nd 0.2
    end
end
e = c*o; % # of edges


% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             OBSTACLE LIST:
% 
% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% Finding subscripts of the obstacles.
OBSenv.ind            = find(~map_env);     % linear indices of obstacles.
[OBSenv_r,OBSenv_c]   = find(~map_env);     % rows and columns of obstacles.
OBSenv.sub            = [OBSenv_r OBSenv_c];  % subscripts of obstacles.

OBSconf.ind           = find(~map_conf);     % linear indices of obstacles.
[OBSconf_r,OBSconf_c] = find(~map_conf);     % rows and columns of obstacles.
OBSconf.sub           = [OBSconf_r OBSconf_c];  % subscripts of obstacles.


% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                             EDGES AND COST RELATED MATRICS:
% 
% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
disp('ETAC...')

[  T,...
   E,...
   eE,...
   A,...
   oA,...
   Ct,...
   P ] = fETAC( map_env,...
                cnt,...
                n,...
                f,...
                c,...
                o,...
                e,...
                oO,...
                eO,...
                cellsize_env,...
                OBSenv,...
                para_);


% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
% 
%                                   VISIBILITY MATRIX:
% 
% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
disp('V...')
[  V,...
   oV,...
   FoVfaces,...
   confKept,...
   confRemoved,...
   cellKept,...
   cellRemoved ] = fV( map_env,...
                       map_conf,...
                       n,...
                       f,...
                       c,...
                       FoV,...
                       FoVfaces,...
                       SensorRange,...
                       OBSenv,...
                       OBSconf,...
                       para_ );

% size(V)
% numel(find(map_env))
% numel(find(map_conf))
% pause
% -----------------------------------------------------------------------------------
% Total computation time.
tPreprocess =  1e-4*(round(toc(tPreprocessi)*1e4));

disp(['tPreprocess: ',num2str(tPreprocess),' sec']);

end


function [ T,...
           E,...
           eE,...
           A,...
           obsA,...
           Ct,...
           P ] = fETAC( map_env,...
                        cnt,...
                        n,...
                        f,...
                        c,...
                        o,...
                        e,...
                        oO,...
                        eO,...
                        cell_length_m,...
                        OBS,...
                        para_)
                        
%fTt returns traversing time matrix Ct and predecessor matrix P (which is not
%used so far).
%

% CONSTANTS:
% ----------
[map_row,map_col] = size(map_env);

ID1oO = find(oO(1,:)); % ID 1
ID2oO = find(oO(2,:)); % ID 2
ID3oO = find(oO(3,:)); % ID 3
ID4oO = find(oO(4,:)); % ID 4
ID5oO = find(oO(5,:)); % ID 5
if ~all(size(ID1oO))
    ID1oO = 0;
end
if ~all(size(ID2oO))
    ID2oO = 0;
end
if ~all(size(ID3oO))
    ID3oO = 0;
end
if ~all(size(ID4oO))
    ID4oO = 0;
end
if ~all(size(ID5oO))
    ID5oO = 0;
end
ID1eO = find(eO(1,:)); % ID 1
ID2eO = find(eO(2,:)); % ID 2
ID3eO = find(eO(3,:)); % ID 3
ID4eO = find(eO(4,:)); % ID 4
ID5eO = find(eO(5,:)); % ID 5
if ~all(size(ID1eO))
    ID1eO = 0;
end
if ~all(size(ID2eO))
    ID2eO = 0;
end
if ~all(size(ID3eO))
    ID3eO = 0;
end
if ~all(size(ID4eO))
    ID4eO = 0;
end
if ~all(size(ID5eO))
    ID5eO = 0;
end

TstepS = para_.TravCost*cell_length_m*1.000; % translational step - straight
TstepD = para_.TravCost*cell_length_m*1.400; % Translational step - diagonal
Rstep1 = para_.TravCost*cell_length_m*0.125; % rotational step - 1st deg
Rstep2 = para_.TravCost*cell_length_m*0.250; % rotatioanl step - 2nd deg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize T:
T = zeros(c,o);

for i = 1:n
    if map_env(i)
        [rr,col] = ind2sub(size(map_env),i); % index number            
        
        if cnt == 4 && f~=8
            % -- List of all neighbors - straight
            N = [rr+0 col+1;...
                 rr-1 col+0;...
                 rr+0 col-1;...
                 rr+1 col+0];
        elseif cnt == 8 || f == 8
            N = [rr+0 col+1;...
                 rr-1 col+1;...
                 rr-1 col+0;...
                 rr-1 col-1;...
                 rr+0 col-1;...
                 rr+1 col-1;...
                 rr+1 col+0;...
                 rr+1 col+1];
        end
        
        for j = 1:f
            oCONFn = (i-1)*f+j; % conf number (from/outgoing)

            for k = 1:o
                l = (j*(f~=1))+(k*(f==1)); % selecting neighbor cell.
                   
                if (any(ID1oO==k)) || (any(ID2oO==k)) || (ID2eO==k)
                    
                    if ~(N(l,1) <= 0 || N(l,1) > map_row ||...
                           N(l,2) <= 0 || N(l,2) > map_col ||...
                           (~~ID3eO*(~mod(l,2))) ||...
                           map_env(N(l,1),N(l,2)) == 0)
                       
                       T(oCONFn,k) = TstepS;
                    else
                       T(oCONFn,k) = 1e300;
                    end
                end
                
                if ((any(ID2oO==k))*(mod(l,2)))||((any(ID2eO==k))*(~mod(l,2)))
                    
                    if ~(N(l,1) <= 0 || N(l,1) > map_row ||...
                           N(l,2) <= 0 || N(l,2) > map_col ||...
                           map_env(N(l,1),N(l,2)) == 0)
                       
                       T(oCONFn,k) = TstepD;
                    else
                       T(oCONFn,k) = 1e300;
                    end
                end
                
                if (any(ID5oO==k)||any(ID4oO==k)) && (k==2||k==3)
                    T(oCONFn,k) = (Rstep1*(any(ID4oO==k)))+(Rstep2*(any(ID5oO==k)));
                end
            end
        end
    end
end

% OBSTACLE FREE:
% Removing obstacle indices:
for i = 1:numel(OBS.ind)
    T(((OBS.ind(i)-1)*f+1)-(f*(i-1)):((OBS.ind(i)-1)*f+f)-(f*(i-1)),:)  = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize E:
E = zeros(e,1);

for i = 1:n
    if map_env(i)
        [rr,col] = ind2sub(size(map_env),i); % index number            
        
        if cnt == 4 && f~=8
            % -- List of all neighbors - straight
            N = [rr+0 col+1;...
                 rr-1 col+0;...
                 rr+0 col-1;...
                 rr+1 col+0];
        elseif cnt == 8 || f == 8
            N = [rr+0 col+1;...
                 rr-1 col+1;...
                 rr-1 col+0;...
                 rr-1 col-1;...
                 rr+0 col-1;...
                 rr+1 col-1;...
                 rr+1 col+0;...
                 rr+1 col+1];
        end

        for j = 1:f            
            oCONFn = (i-1)*f+j; % conf number (from/outgoing)             

            for k = 1:o
                oEDGEn = (oCONFn-1)*o+k;   % edge num (from/outgoing)
                l = (j*(f~=1))+(k*(f==1)); % selecting neighbor cell.

                if (any(ID1oO==k)) || (any(ID2oO==k)) || (ID2eO==k)
                    
                    if ~(N(l,1) <= 0 || N(l,1) > map_row ||...
                           N(l,2) <= 0 || N(l,2) > map_col ||...
                           map_env(N(l,1),N(l,2)) == 0)
                       
                       Ni = sub2ind(size(map_env),N(l,1),N(l,2));
                       iCONFn = (Ni-1)*f+j;         % conf num (to/incoming)
                       iEDGEn = (iCONFn-1)*o+k;     % edge num (to/incoming)
                       E(oEDGEn) = iEDGEn;
                       
                    else
                       iCONFn = oCONFn+(((f/2)-((j/f>0.5)*f))*~mod(f,2));
                       iEDGEn = ((oEDGEn+(((o/2)-((k/o>0.5)*o))*~mod(o,2)))*(iCONFn==oCONFn))+...
                           (((iCONFn-1)*o+k)*(iCONFn~=oCONFn));
                       E(oEDGEn) = iEDGEn;
                    end                    
                end

                if (any(ID5oO==k)||any(ID4oO==k)) && (k==2)
                    % neighbor to the 2nd edge
                    iCONFn = oCONFn+((1-((~mod(j,f))*f))*~mod(f,2));  % incoming conf#                    
                    iEDGEn = (iCONFn-1)*o+k;
                    E(oEDGEn) = iEDGEn;
                end

                if (any(ID5oO==k)||any(ID4oO==k)) && (k==3)                    
                    % neighbor to the 3rd edge
                    iCONFn = oCONFn-1+((((mod(j,f)==1)*f))*~mod(f,2)); % incoming conf#
                    iEDGEn = (iCONFn-1)*o+k;
                    E(oEDGEn) = iEDGEn;                     
                end
            end
        end
    end
end


% Rearranging IDs of connected neighbours for E:
for i = numel(OBS.ind):-1:1
    E(E>(OBS.ind(i)*(f*o)))= E(E>(OBS.ind(i)*(f*o)))-(f*o);
end

% OBSTACLE FREE:
% Removing obstacle indices:
for i = 1:numel(OBS.ind)
    E(((OBS.ind(i)-1)*(f*o)+1)-((f*o)*(i-1)):((OBS.ind(i)-1)*(f*o)+(f*o))-((f*o)*(i-1)),:) = [];
end

eE = sub2ind([numel(find(map_env))*f o],fix(E/o)+(~~mod(E,o)),mod(E,o)+(~mod(E,o))*o);
% eE = sub2ind([c o],fix(E/o)+(~~mod(E,o)),mod(E,o)+(~mod(E,o))*o);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A    = []; % initialize A.
obsA = []; % initialize A with obstacles.

if para_.A
    
    % Initialize A:
    A = inf(length(T(:,1))); % for all obstacle free conf
    for i = 1:length(T(:,1))
        a = find((T(i,:))<=TstepD);
        for j = 1:numel(a)        
            oEDGE = (i-1)*o+a(j);
            iEDGE = E(oEDGE);
            iCONF = fix(iEDGE/o)+(~~mod(iEDGE,o));
            A(i,iCONF) = T(i,a(j));
        end
    end
    obsA = A; % A with obstacles.
    for i = 1:numel(OBS.ind)
        obsA = [obsA(:,1:(OBS.ind(i)-1)*f) inf(length(obsA(:,1)),f) ...
            obsA(:,(OBS.ind(i)-1)*f+1:end)];
        obsA = [obsA(1:(OBS.ind(i)-1)*f,:); inf(f,length(obsA(1,:))); ...
            obsA((OBS.ind(i)-1)*f+1:end,:)];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Ct and P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ct = []; P = [];

if para_.Ct

    % GENERATING MATRIX Ct and P:
    % ---------------------------
    Ct = inf(c); % initiate Ct 
    P  = inf(c); % initiate P

    originalA = A; % save A
    %[A.rp,A.ci,A.ai]=sparse_to_csr(obsA);
    [rp,ci,ai]=sparse_to_csr(obsA);
    A = [];
    A.rp = rp; A.ci = ci; A.ai = ai;

    for i = 1:c
        %[d, pred] = dijkstra(obsA,i);
        [d, pred] = dijkstra(A,i);
        Ct(i,:)   = d;
        P(i,:)    = pred;
    end
    A = [];
    A = originalA; % recover original A.

    % Removing OBS
    % ------------
    for i = 1:numel(OBS.ind)    
        % Ct
        Ct(:,((OBS.ind(i)-1)*f+1)-(f*(i-1)):((OBS.ind(i)-1)*f+f)-(f*(i-1))) = []; % from
        Ct(((OBS.ind(i)-1)*f+1)-(f*(i-1)):((OBS.ind(i)-1)*f+f)-(f*(i-1)),:) = []; % to
        % P
        P(:,((OBS.ind(i)-1)*f+1)-(f*(i-1)):((OBS.ind(i)-1)*f+f)-(f*(i-1))) = []; % from
        P(((OBS.ind(i)-1)*f+1)-(f*(i-1)):((OBS.ind(i)-1)*f+f)-(f*(i-1)),:) = []; % to
    end
end

% % Number of edges in the map_env
% nEDGES = numel(find(A~=inf & A~=0));


end


function [ V,...
           oV,...
           FoVfaces,...
           confKept,...
           confRemoved,...
           cellKept,...
           cellRemoved ] = fV( map_env,...
                               map_conf,...
                               n,...
                               f,...
                               c,...
                               FoV,...
                               FoVfaces,...
                               SensorRange,...
                               OBSenv,...
                               OBSconf,...
                               para_ )
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

FoVfaces.lgt
FoVfaces.ang

%********************************************
%         SENSOR DOMAIN TEMPLATE
%********************************************
disp('.. Sensor Domain')
% fSensorDomainTemplate(FoVfaces,SensorRange,f,para_);

%********************************************
%            VISIBILITY MATRIX
%********************************************
disp('.. Visibility Matrix')
[V] = fVisbilityMatrix(n,f,c,FoVfaces,SensorRange,FoV,map_conf,map_env);

oV = V;
% V = sparse(V);


% unoccupied cells
free_cell_vector = ones(1,size(oV,1));
free_cell_vector(OBSenv.ind) = 0;

cellKept     = find(free_cell_vector);
cellRemoved  = find(free_cell_vector==0);
V            = V(cellKept,:);


% ---- configurations:

% --- selecting conf number corresponding to candidate conf map only
% cells for the candidate conf.
candConfCells = find(map_conf); 
% conf numbers
% candConfNum = zeros(numel(map_candConf)*f,1);
candConfNum = zeros(1,size(oV,2));
for i = 1:numel(candConfCells)
    candConfNum( (candConfCells(i)-1)*f+1:(candConfCells(i)-1)*f+f ) = 1;
end

% conf belong to unoccupied cells
conf_free = ones(1,size(oV,2));
for i = 1:numel(OBSconf.ind)
    conf_free( (OBSconf.ind(i)-1)*f+1:(OBSconf.ind(i)-1)*f+f ) = 0;
end


% disp('here')
% find(candConfNum~=conf_free)
% disp('and here')
% pause

% --- keeping candidate conf in the list that are from candidate conf map
% AND cover the hotspot AND are placed over unoccupied cells 
conf_vector  = candConfNum.*conf_free;
confKept     = find(conf_vector);
confRemoved  = find(conf_vector==0);
V            = V(:,confKept);
% V            = full(V);



%{
% OBSTACLE FREE VISIBILITY:
% ------------------------------------------------------------------------------
% for i = 1:numel(OBS.ind)
%     V(OBS.ind(i)-(i-1),:) = []; % removing cells.
%     % removing sensing configurations.
%     V(:,((OBS.ind(i)-1)*f+1)-(f*(i-1)):((OBS.ind(i)-1)*f+f)-(f*(i-1))) = [];
% end

V(OBS.ind,:) = []; % removing cells
% removing conf.
confOBS = zeros(numel(OBS.ind)*f,1);
for i = 1:numel(OBS.ind)
    % removing sensing configurations.
    %V(:,((OBS.ind(i)-1)*f+1)-(f*(i-1)):((OBS.ind(i)-1)*f+f)-(f*(i-1))) = [];
    confOBS(((i-1)*f)+1:((i-1)*f)+f) = (((OBS.ind(i)-1)*f)+1):(((OBS.ind(i)-1)*f)+f);
end
V(:,confOBS) = [];

% V = sparse(V);

%}


end




function [ CRL ] = fCRL( o,f,A,E,V,Ct,map_env,OBS )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Ei
Ei    = zeros(numel(E),1);
Ei(E) = 1:numel(E);

% Set of vertices belong to the boundry cells in the map_env.
erodeScale        = 1; % scale to erode.
mapBndry          = map_env-bwmorph(map_env,'erode',erodeScale); % map with boundry cells.
mapBndry(OBS.ind) = []; % remove obstacle cells.
candCell          = find(mapBndry); % Candidate cells for CRL.

candCONF = zeros(numel(candCell)*f,1); % initialize candidate conf for CRL.
for i = 1:numel(candCell)
    candCONF(((i-1)*f)+(1:f)) = (((candCell(i)-1)*f)+(1:f));
end

% --- CORNER REDUCTION LIST:
CRL       = []; % initialize corner reduction list (CRL).
closeList = []; % initialize close list.

while numel(candCONF) % run till candidate conf is null.
    
    i = candCONF(1);                % pick the first candidate conf 'Ci'.
    candCONF(candCONF==i) = [];     % remove it from the list.
    closeList = union(closeList,i); % update close list.
    
    % Edge numbers of the candidate conf.
    oEDGES = ((i-1)*o)+(1:o);
    iEDGES = oEDGES;
    
    % List of incoming vertices/conf to the subject conf.
    iCONF = fix((Ei(iEDGES))/o)+(~~mod((Ei(iEDGES)),o));
    iCONF = unique(iCONF);
    iCONF(iCONF==i) = [];
    
    % List of outgoing vertices/conf from the subject conf.
    oCONF = fix((E(oEDGES))/o)+(~~mod((E(oEDGES)),o));
    oCONF = unique(oCONF);
    oCONF(oCONF==i) = [];
    
    % Pairs of incoming and outgoing vertices/conf of the subject conf
    % 'Ci'.
    [p,q] = meshgrid(iCONF,oCONF);
    pairs = [p(:) q(:)];
    pairs(eq(pairs(:,1),pairs(:,2)),:) = []; % remove pairs of same vertices.
    
    % coverage set of the subject conf 'Ci'.
    Oi = find(V(:,i));
    
    % initiate condition array for all pairs.
    Cond = zeros(1,length(pairs(:,1)));

    for l = 1:length(pairs(:,1)) % for all pairs.

        Oj = find(V(:,pairs(l,1))); % set of visible cells from Cj.
        Ok = find(V(:,pairs(l,2))); % set of visible cells from Ck.
        
        if all(ismember(Oi,union(Oj,Ok))) % Oi is in OjUOk (first condition).
            
            % --- Distance before removing Ci.
            bTjk = Ct(pairs(l,1),pairs(l,2));

            % --- Distance after removing Ci.
            Ax          = A;    % some temporary A matrix
            Ax(i,:)     = inf;  % inf weight from Ci to all
            Ax(:,i)     = inf;  % inf weight from all to Ci
            [d,~]       = dijkstra(Ax,pairs(l,1)); % distance from Cj to all
            aTjk        = d(pairs(l,2)); % distance from Cj to Ck

            % --- Update condition array
            Cond(l) = bTjk==aTjk;

            if ~(bTjk==aTjk) % if second condition does'nt hold.
                break % no need to check rest of the pairs
            end            
        end
    end
    if all(Cond) % if second cond hols for all pairs.
        CRL     = union(CRL,i);     % add Ci in the list.
        
        % ---- Update candidate conf list by adding neighbour conf.s.        
        %oEDGES  = ((i-1)*o)+(1:o);  % outgoing edges of i.
        iEDGES  = E(oEDGES); % incoming edges.
        % Incoming vertices/conf.
        iCONF   = fix(iEDGES/o)+(~~mod(iEDGES,o));
        iCONF   = unique(iCONF);
        iCONF(iCONF==i) = [];
        % insert neighbor conf to the candidate conf list if they are not
        % already in the closed list.
        candCONF = union(candCONF,iCONF(~ismember(iCONF,closeList)));
    end
end


end


function [ vV,vT,vE,evE ] = fvVTE( f,o,V,T,E,CRL,map_env )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% vV
vV = V;
vV(:,CRL) = [];

%%                                vT & vE
% ----------------------------------------------------------------------------------------

% eCRL = zeros(numel(CRL)*o,1);
eCRL = [];
for i = 1:numel(CRL)
    eCRL = [eCRL (CRL(i)-1)*o+(1:o)];
end

vT = T;
vE = E;
for i = 1:numel(eCRL)
    a = find((E(:)==eCRL(i)));
    if ~any(a==eCRL)
        vE(a) = vE(vE(a));
        if o == 4
            vE(a) = vE(vE(a));
        end
    end
    vT(fix(a/o)+(~~mod(a,o)),mod(a,o)+(~mod(a,o)*o)) = 1e300;
end
vE(eCRL)  = [];
vT(CRL,:) = [];

% Rearranging IDs of connected neighbours for vE:
for i = numel(eCRL):-1:1
    vE(vE>(eCRL(i)))= vE(vE>(eCRL(i)))-1;
end

evE = sub2ind([(numel(find(map_env)))*f o],...
    fix(vE/o)+(~~mod(vE,o)),mod(vE,o)+(~mod(vE,o))*o);


end


function [ vkS ] = fvkS( V )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Select visibility matrix
% if DIM.CnrRed
%     Visibility = V;%vV;
% elseif ~DIM.CnrRed
%     Visibility = V;
% end
Visibility = V;

% Number of conf for each cell
Vk = sum(Visibility,2);

% Pick a cell observed by min conf and list those conf for special vertices
vkS.cell = find(Vk==min(Vk));
vkS.cell = vkS.cell(1);
vkS.conf = find(Visibility(vkS.cell,:));

end


% ------------------------- End of the document -------------------------------------

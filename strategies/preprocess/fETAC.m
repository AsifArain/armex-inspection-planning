function [ T,E,eE,A,oA,Ct,P ] = fETAC( map_env,...
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
                                       para_ )
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


%%                                    T
% ----------------------------------------------------------------------------------------
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

%%                                    E
% ----------------------------------------------------------------------------------------

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

%%                                    A
% ----------------------------------------------------------------------------------------
A  = []; % initialize A.
oA = []; % initialize A with obstacles.

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
    oA = A; % A with obstacles.
    for i = 1:numel(OBS.ind)
        oA = [oA(:,1:(OBS.ind(i)-1)*f) inf(length(oA(:,1)),f) ...
            oA(:,(OBS.ind(i)-1)*f+1:end)];
        oA = [oA(1:(OBS.ind(i)-1)*f,:); inf(f,length(oA(1,:))); ...
            oA((OBS.ind(i)-1)*f+1:end,:)];
    end
end

%%                             Ct and P
% ----------------------------------------------------------------------------------------
Ct = []; P = [];

if para_.Ct

    % GENERATING MATRIX Ct and P:
    % ---------------------------
    Ct = inf(c); % initiate Ct 
    P  = inf(c); % initiate P

    originalA = A; % save A
    %[A.rp,A.ci,A.ai]=sparse_to_csr(obsA);
    [rp,ci,ai]=sparse_to_csr(oA);
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

% % Number of edges in the map
% nEDGES = numel(find(A~=inf & A~=0));


end

% ---- lets keep new fETAC
%{
function [ T,E,eE,A,oA,D,P ] = fETAC( map_env,...
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
                                      startConf,...
                                      confKept,...
                                      paramPreprocess )
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

TstepS = paramPreprocess.TravCost*cell_length_m*1.000; % translational step - straight
TstepD = paramPreprocess.TravCost*cell_length_m*1.400; % Translational step - diagonal
Rstep1 = paramPreprocess.TravCost*cell_length_m*0.125; % rotational step - 1st deg
Rstep2 = paramPreprocess.TravCost*cell_length_m*0.250; % rotatioanl step - 2nd deg


%%                                    T
% ----------------------------------------------------------------------------------------
% Initialize T:
T = zeros(c,o);
size(T)
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

oT = T; % keep T with obs
size(oT)

% OBSTACLE FREE:
% Removing obstacle indices:
for i = 1:numel(OBS.ind)
    T(((OBS.ind(i)-1)*f+1)-(f*(i-1)):((OBS.ind(i)-1)*f+f)-(f*(i-1)),:)  = [];
end

size(T)

%%                                    E
% ----------------------------------------------------------------------------------------

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

%%                                    A
% ----------------------------------------------------------------------------------------
A  = []; % initialize A.
oA = []; % initialize A with obstacles.

if paramPreprocess.A
    
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
    oA = A; % A with obstacles.
    for i = 1:numel(OBS.ind)
        oA = [oA(:,1:(OBS.ind(i)-1)*f) inf(length(oA(:,1)),f) ...
            oA(:,(OBS.ind(i)-1)*f+1:end)];
        oA = [oA(1:(OBS.ind(i)-1)*f,:); inf(f,length(oA(1,:))); ...
            oA((OBS.ind(i)-1)*f+1:end,:)];
    end

end

%%                             Distance (D) and Predecessor (P) matrix
% ----------------------------------------------------------------------------------------
D = []; P = [];

if paramPreprocess.D

    % GENERATING MATRIX Ct and P:
    % ---------------------------
    %Ct = inf(c); % initiate Ct 
    %P  = inf(c); % initiate P

    originalA = A; % save A
    %[A.rp,A.ci,A.ai]=sparse_to_csr(obsA);
    [rp,ci,ai]=sparse_to_csr(oA);
    A = [];
    A.rp = rp; A.ci = ci; A.ai = ai;

    %for i = 1:c
    %    %[d, pred] = dijkstra(obsA,i);
    %    [d, pred] = dijkstra(A,i);
    %    Ct(i,:)   = d;
    %    P(i,:)    = pred;
    %end
    oA
    size(oA)
    [D,P] = dijkstra(oA,startConf);
    
    A = [];
    A = originalA; % recover original A.

    % Removing OBS
    % ------------
    %{
    for i = 1:numel(OBS.ind)    
        % Ct
        Ct(:,((OBS.ind(i)-1)*f+1)-(f*(i-1)):((OBS.ind(i)-1)*f+f)-(f*(i-1))) = []; % from
        Ct(((OBS.ind(i)-1)*f+1)-(f*(i-1)):((OBS.ind(i)-1)*f+f)-(f*(i-1)),:) = []; % to
        % P
        P(:,((OBS.ind(i)-1)*f+1)-(f*(i-1)):((OBS.ind(i)-1)*f+f)-(f*(i-1))) = []; % from
        P(((OBS.ind(i)-1)*f+1)-(f*(i-1)):((OBS.ind(i)-1)*f+f)-(f*(i-1)),:) = []; % to
    end
    %}
    
    % update distance and predecessor matrix for the candidate conf only.
    D = D(confKept);
    P = P(confKept);
    
    
end

% % Number of edges in the map
% nEDGES = numel(find(A~=inf & A~=0));


end
%}


% --------------------------- End of the document -----------------------------------

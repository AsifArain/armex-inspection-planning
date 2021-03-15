function [ gCost,gPred,tAstar ] = fAstar( n,f,o,sConf,gConf,E,Eo,T,map,OBS,ParamAstar )
%fAstar is A* algorithm.
% Dated: 2014-01-17, Rev 1
% 
% --------------------------------------------------------------------------------------------------------------------------------
% OUTPUTS:
% --------------------------------------------------------------------------------------------------------------------------------
% gCost  : [scalar][sec] Cost to travel from start conf to goal conf.
% gPred  : [1,n][#] Predecessor to move from goal conf to start conf.
% tAstar : [scalar][sec] Computation time for fAstar.
% 
% --------------------------------------------------------------------------------------------------------------------------------
% INPUTS:
% --------------------------------------------------------------------------------------------------------------------------------
% n            : [ ][ ] 
% f            : [ ][ ] 
% o            : [ ][ ] 
% sConf        : [ ][ ] 
% gConf        : [ ][ ] 
% E            : [ ][ ] 
% T            : [ ][ ] 
% map          : [ ][ ] 
% OBS          : [ ][ ]
% ParamAstar   : [][] A star parameters.
%                     .TravCost := travelling cost (sec) for one meter travel.
%                     .ConfOBS  := start and goal conf. are with/without obstacle indices (yes/no).
% 
% --------------------------------------------------------------------------------------------------------------------------------
% NOTES/UPDATES:
% --------------------------------------------------------------------------------------------------------------------------------
% R1     :
%  

% Queue        := queue list.
% minConf      := conf with minimum cost value.
% gConf        := goal conf.
% sConf        := start conf.
% Cst          := cost list for all conf.
% Key          := key list for all conf.

tAstari         = tic;
if ~ParamAstar.ConfOBS
    for i = 1:numel(OBS.ind)
        sConf(sConf>((OBS.ind(i)-1)*f))= sConf(sConf>((OBS.ind(i)-1)*f))+f;
        gConf(gConf>((OBS.ind(i)-1)*f))= gConf(gConf>((OBS.ind(i)-1)*f))+f;
    end
end

Queue           = sConf; % sConf := start conf
minConf         = sConf;

% COST:
Cst             = inf(1,n*f);
Cst(sConf)      = 0;

% KEY:
Key             = inf(1,n*f);
d               = ParamAstar.TravCost;                                           % manhattan distance between two cells.
sCell           = fix(sConf/f)+(~~mod(sConf,f));
gCell           = fix(gConf/f)+(~~mod(gConf,f));
[sCellR,sCellC] = ind2sub(size(map),sCell);
[gCellR,gCellC] = ind2sub(size(map),gCell);
Key(sConf)      = Cst(sConf)+(d*(abs(sCellR-gCellR)+abs(sCellC-gCellC)));
                             % (d*(abs(sCellR-gCellR)+abs(sCellC-gCellC)));      % Manhattan Distance
                             % sqrt(((gCellR-sCellR)^2)+((gCellC-sCellC)^2));    % Euclidean Distance

% PREDECESSOR:
Pred            = inf(1,n*f);
Pred(sConf)     = 0;

% CLOSED LIST:
ClosedList      = [];

% If Eo and T are not already computed.
if ~ParamAstar.EoT
    % Eo & T:
    Eo              = zeros(size(T'));
    i               = 1:length(E);
    Eo(i)           = fix(E(i)/o)+(~~mod(E(i),o));
    Eo              = Eo';
    % Fixing Eo and T for not connected Conf(s).
    i = 1:length(Eo(:,1));
    T((fix(Eo(i,1)/f)+(~~mod(Eo(i,1),f)))'==fix(i/f)+(~~mod(i,f)),1)  = inf;
    Eo((fix(Eo(i,1)/f)+(~~mod(Eo(i,1),f)))'==fix(i/f)+(~~mod(i,f)),1) = inf;

    % Insert obstacles.
    for i = 1:numel(OBS.ind)
        T = [T(1:(OBS.ind(i)-1)*f,:); inf(f,length(T(1,:))); T((OBS.ind(i)-1)*f+1:end,:)];
        Eo = [Eo(1:(OBS.ind(i)-1)*f,:); inf(f,length(Eo(1,:))); Eo((OBS.ind(i)-1)*f+1:end,:)];    
    end
    % Rearrange IDs after inserting Obstacles.
    for i = 1:numel(OBS.ind)
        Eo(Eo>((OBS.ind(i)-1)*f))= Eo(Eo>((OBS.ind(i)-1)*f))+f;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A-STAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while minConf~=gConf
    
    [~,minInd]          = min(Key(Queue));
    %minInd              = argmin(Key(Queue));
    minConf             = Queue(minInd);
    %ClosedList          = union(ClosedList,minConf);
    ClosedList          = union_sorted(ClosedList,minConf);
    nConf               = Eo(minConf,:);
    nConf(nConf==inf)   = [];
    
    for i = 1:numel(nConf)
        
        %if ~ismember(nConf(i),ClosedList)
        if ~ismember_sorted(nConf(i),ClosedList)
            testCost = Cst(minConf)+T(minConf,Eo(minConf,:)==nConf(i));
            
            if Cst(nConf(i)) > testCost
                
                Cst(nConf(i))   = testCost;
                %Queue           = union(Queue,nConf(i));
                Queue           = union_sorted(Queue,nConf(i));
                nCell           = fix(nConf(i)/f)+(~~mod(nConf(i),f));
                [nCellR,nCellC] = ind2sub(size(map),nCell);
                
                Key(nConf(i))   = Cst(nConf(i))+(d*(abs(nCellR-gCellR)+abs(nCellC-gCellC)));
                                  %sqrt(((gCellR-nCellR)^2)+((gCellC-nCellC)^2)); % Euclidean Distance
                                  %(d*(abs(nCellR-gCellR)+abs(nCellC-gCellC)));   % Manhattan Distance
                Pred(nConf(i))  = minConf;
            end
        end
    end
    [~,~,IQ]  = intersect(minConf,Queue);
    %IQ  = intersect_sorted(minConf,Queue);
    Queue(IQ) = []; % remove minimum cost conf from Queue.
end

% numel(ClosedList)

gCost = Cst(gConf);
gPred = Pred(gConf);

if gCost
    while gPred(end) ~= sConf    
        gPred = [gPred Pred(gPred(end))];
    end
end

tAstar = 1e-4*(round(toc(tAstari)*1e4));

end

function [c,a_match,b_match] = union_sorted(a,b)
%UNION_SORTED  Set union of sorted sets.
% UNION_SORTED(A,B) where A and B are vectors returns the combined values
% from A and B with no repetitions.  A (and B) must be sorted and unique, and 
% the result will be sorted and unique.
%
% [C,A_MATCH,B_MATCH] = UNION_SORTED(A,B) also returns
%   A_MATCH = MATCH_SORTED(A,C)
%   B_MATCH = MATCH_SORTED(B,C)
%
% Examples:
%   union_sorted([20 30 40], [10 20 30])
%   [c,a_match,b_match] = union_sorted([20 30 40], [10 20 30])

% Written by Tom Minka
% (c) Microsoft Corporation. All rights reserved.

% instead of a full sort, you could do a merge of the two sorted lists.
if nargout <= 1
  c = sort([a(~ismember_sorted(a,b)) b]);
else
  [tf,loc] = ismember_sorted(a,b);
  [c,i] = sort([a(~tf) b]);
  % c = [a(~tf) b](i)
  nc = length(c);
  nb = length(b);
  na = length(a);
  index = zeros(1,nc);
  index(i) = 1:nc;
  % c(index) = [a(~tf) b]
  b_match = index((nc-nb+1):nc);
  a_match = zeros(1,na);
  a_match(~tf) = index(1:(nc-nb));
  a_match(tf) = b_match(loc(tf));
end
end

function [tf,loc] = ismember_sorted(a,s)
%ISMEMBER_SORTED   True for member of sorted set.
% ISMEMBER_SORTED(A,S) for the vector A returns an array of the same size as A
% containing 1 where the elements of A are in the set S and 0 otherwise.
% A and S must be sorted and cannot contain NaN.
%
% [TF,LOC] = ISMEMBER_SORTED(A,S) also returns an index array LOC where
% LOC(i) is the index in S which matches A(i) (highest if there are ties)
% or 0 if there is no such index.
%
% See also ISMEMBER, MATCH_SORTED, INTERSECT_SORTED, SETDIFF_SORTED, UNION_SORTED.

% Written by Tom Minka
% (c) Microsoft Corporation. All rights reserved.

% The internal function ismembc comes from ismember.m
% It requires non-sparse arrays.
a = full(a);
s = full(s);
if nargout < 2
  tf = ismembc(a,s);
else
  loc = ismembc2(a,s);
  tf = (loc > 0);
end
end

%% TRASH - fAstar.m
%{
function [ gCost,gPred ] = fAstar( n,f,o,sConf,gConf,E,T,map,OBS )
%fAstar is A* algorithm.
% 

% Queue     := queue list.
% minConf   := conf with minimum cost value.
% gConf     := goal conf.
% sConf     := start conf.
% Cst       := cost list for all conf.
% Key       := key list for all conf.

%  -----

% for i = 1:numel(OBS.ind)
%     % Configuration:
%     T = [T(1:(OBS.ind(i)-1)*f,:); zeros(f,length(T(1,:))); ...
%         T((OBS.ind(i)-1)*f+1:end,:)];
%     
%     E = [E(1:(OBS.ind(i)-1)*f,:); zeros(f,length(E(1,:))); ...
%         E((OBS.ind(i)-1)*f+1:end,:)];    
% end

% for i = 1:numel(OBS.ind)
%     E(E>((OBS.ind(i)-1)*f))= E(E>((OBS.ind(i)-1)*f))+f;
%     sConf(sConf>((OBS.ind(i)-1)*f))= sConf(sConf>((OBS.ind(i)-1)*f))+f;
%     gConf(gConf>((OBS.ind(i)-1)*f))= gConf(gConf>((OBS.ind(i)-1)*f))+f;
% end



figure;
[map_row, map_col] = size(map);
for i = 1:map_row
    for j = 1:map_col
        % if obstacle/free
        if ~map(i,j) % occupied
            col = [0.1 0.1 0.1];
        elseif map(i,j) % free
            col = [0.9 0.9 0.9];
        end
        plot(i-0.5,j-0.5,'S', 'color',col,'MarkerFaceColor',col,'MarkerSize',50); hold on;
    end
end
axis equal
% ----


Queue           = sConf; % sConf := start conf
minConf         = sConf;





% COST:
Cst             = inf(1,n*f);
Cst(sConf)      = 0;

% KEY:
Key             = inf(1,n*f);
d               = 10; % manhattan distance between two cells.
sCell           = fix(sConf/f)+(~~mod(sConf,f));
gCell           = fix(gConf/f)+(~~mod(gConf,f))
[sCellR,sCellC] = ind2sub(size(map),sCell);
[gCellR,gCellC] = ind2sub(size(map),gCell)
Key(sConf)      = Cst(sConf)+(d*(abs(sCellR-gCellR)+abs(sCellC-gCellC)));
% (d*(abs(sCellR-gCellR)+abs(sCellC-gCellC)));      % Manhattan Distance
% sqrt(((gCellR-sCellR)^2)+((gCellC-sCellC)^2));    % Euclidean Distance

% PREDECESSOR:
Pred            = inf(1,n*f);
Pred(sConf)     = 0;

% CLOSED LIST:
ClosedList      = [];

% Eo & T:
Eo              = zeros(size(T'));
i               = 1:length(E);
Eo(i)           = fix(E(i)/o)+(~~mod(E(i),o));
Eo              = Eo';
% Fixing Eo and T for not connected Conf(s).
i = 1:length(Eo(:,1));
T((fix(Eo(i,1)/f)+(~~mod(Eo(i,1),f)))'==fix(i/f)+(~~mod(i,f)),1) = inf;
Eo((fix(Eo(i,1)/f)+(~~mod(Eo(i,1),f)))'==fix(i/f)+(~~mod(i,f)),1) = inf;
for i = 1:numel(OBS.ind)
    % Configuration:
    T = [T(1:(OBS.ind(i)-1)*f,:); inf(f,length(T(1,:))); ...
        T((OBS.ind(i)-1)*f+1:end,:)];
    
    Eo = [Eo(1:(OBS.ind(i)-1)*f,:); inf(f,length(Eo(1,:))); ...
        Eo((OBS.ind(i)-1)*f+1:end,:)];    
end
for i = 1:numel(OBS.ind)
    Eo(Eo>((OBS.ind(i)-1)*f))= Eo(Eo>((OBS.ind(i)-1)*f))+f;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A-STAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while minConf~=gConf
    
    [~,minInd]          = min(Key(Queue));
    minConf             = Queue(minInd);
    disp('------------')
    Queue
    CCost = Cst(Queue)
    KKey  = Key(Queue)
    minConf
    % -----
    mConf = minConf;
    %for i = 1:numel(OBS.ind)
    %    mConf(mConf>((OBS.ind(i)-1)*f))= mConf(mConf>((OBS.ind(i)-1)*f))+f;
    %end
    mCell = fix(mConf/f)+(~~mod(mConf,f));
    [mCellR,mCellC] = ind2sub(size(map),mCell);
    if mod(minConf,f) == 1
        mCellC = mCellC+0.25;
    elseif mod(minConf,f) == 2
        mCellR = mCellR-0.25;
    elseif mod(minConf,f) == 3
        mCellC = mCellC-0.25;
    elseif mod(minConf,f) == 0
        mCellR = mCellR+0.25;
    end
    h1(1) = plot(mCellR-0.5,mCellC-0.5,'sg','MarkerSize',15);
    % ----
    %pause
    
    
    ClosedList          = union(ClosedList,minConf);
    % -------
    for c = 1:numel(ClosedList)
        cConf = ClosedList(c);
        %for i = 1:numel(OBS.ind)
        %    cConf(cConf>((OBS.ind(i)-1)*f))= cConf(cConf>((OBS.ind(i)-1)*f))+f;
        %end
        cCell = fix(cConf/f)+(~~mod(cConf,f));
        [cCellR,cCellC] = ind2sub(size(map),cCell);
        if mod(cConf,f) == 1
            cCellC = cCellC+0.25;
        elseif mod(cConf,f) == 2
            cCellR = cCellR-0.25;
        elseif mod(cConf,f) == 3
            cCellC = cCellC-0.25;
        elseif mod(cConf,f) == 0
            cCellR = cCellR+0.25;
        end    
        CLOS(c,:) = [cCellR cCellC];
    end
    h3(1) = plot(CLOS(:,1)-0.5,CLOS(:,2)-0.5,'ok','MarkerFaceColor','k','MarkerSize',7);
    % ------
    %pause
    
    nConf               = Eo(minConf,:);
    nConf(nConf==inf)   = [];
    
    
    
    for i = 1:numel(nConf)
        
        if ~ismember(nConf(i),ClosedList)
            testCost = Cst(minConf)+T(minConf,Eo(minConf,:)==nConf(i));
            
            if Cst(nConf(i)) > testCost
                
                % -----
                % -----
                oConf = nConf(i);
                %for j = 1:numel(OBS.ind)
                %    oConf(oConf>((OBS.ind(j)-1)*f))= oConf(oConf>((OBS.ind(j)-1)*f))+f;
                %end
                oCell = fix(oConf/f)+(~~mod(oConf,f));
                [oCellR,oCellC] = ind2sub(size(map),oCell);
                if mod(nConf(i),f) == 1
                    oCellC = oCellC+0.25;
                elseif mod(nConf(i),f) == 2
                    oCellR = oCellR-0.25;
                elseif mod(nConf(i),f) == 3
                    oCellC = oCellC-0.25;
                elseif mod(nConf(i),f) == 0
                    oCellR = oCellR+0.25;
                end
                xx = [mCellR oCellR]; yy = [mCellC oCellC];
                h4(1) = plot(xx-0.5,yy-0.5,'-k','MarkerSize',2);
                % ----
                % -----
                
                nnConf = nConf(i)
                Cst(nConf(i))   = testCost;
                Queue           = union(Queue,nConf(i));
                nCell           = fix(nConf(i)/f)+(~~mod(nConf(i),f));
                [nCellR,nCellC] = ind2sub(size(map),nCell);
                
                Key(nConf(i))   = Cst(nConf(i))+...
                    (d*(abs(nCellR-gCellR)+abs(nCellC-gCellC)));
                                  %sqrt(((gCellR-nCellR)^2)+((gCellC-nCellC)^2));
                                  %(d*(abs(nCellR-gCellR)+abs(nCellC-gCellC)));
                nnManDist = (d*(abs(nCellR-gCellR)+abs(nCellC-gCellC)))
                Pred(nConf(i))  = minConf;                
            end
        end
    end
    [~,~,IQ]  = intersect(minConf,Queue);
    Queue(IQ) = []; % remove minimum cost conf from Queue.
     for q = 1:numel(Queue)
        qConf = Queue(q);
        %for i = 1:numel(OBS.ind)            
        %    qConf(qConf>((OBS.ind(i)-1)*f))= qConf(qConf>((OBS.ind(i)-1)*f))+f;
        %end
        qCell = fix(qConf/f)+(~~mod(qConf,f));
        [qCellR,qCellC] = ind2sub(size(map),qCell);
        if mod(qConf,f) == 1
            qCellC = qCellC+0.25;
        elseif mod(qConf,f) == 2
            qCellR = qCellR-0.25;
        elseif mod(qConf,f) == 3
            qCellC = qCellC-0.25;
        elseif mod(qConf,f) == 0
            qCellR = qCellR+0.25;
        end    
        Q(q,:) = [qCellR qCellC];
    end
    h2(1) = plot(Q(:,1)-0.5,Q(:,2)-0.5,'ob','MarkerSize',10);
    %pause
    delete(h1(1))
    %delete(h2(1))
    %delete(h3(1))    
    % ----- 
end

gCost = Cst(gConf);
gPred = Pred(gConf);

if gCost
    while gPred(end) ~= sConf    
        gPred = [gPred Pred(gPred(end))];
    end
end


end
%}

%% TRASH - test_Eo.m
%{
clc
Eo = zeros(size(T'));
% for i = 1:length(E)    
%     Eo(i) = fix(E(i)/o)+(~~mod(E(i),o));
% end
i = 1:length(E)    
Eo(i) = fix(E(i)/o)+(~~mod(E(i),o));
Eo = Eo'

i = 1:length(Eo(:,1));    
%     fix(Eo(i,1)/f)+(~~mod(Eo(i,1),f)) == fix(i/f)+(~~mod(i,f))

T((fix(Eo(i,1)/f)+(~~mod(Eo(i,1),f)))'==fix(i/f)+(~~mod(i,f)),1) = inf
Eo((fix(Eo(i,1)/f)+(~~mod(Eo(i,1),f)))'==fix(i/f)+(~~mod(i,f)),1) = inf
%}
function [ C,...
           C0,...
           C1,...
           Cu,...
           combC,...
           lBound,...
           uBound,...
           logs_rwl1,...
           tSPP_rwl1,...
           tSPP_Comb,...
           tSPP ] = fPlan_GD_FullCoverage( V,para_ ) 
       
%fPlan_GD_FullCoverage is Sensor Placement Problem solution based on
%Reweighted l1 minimization. conv-SPP
% Date: 2014-10-17, Rev 5
% 
% n  := total number of cells in the map.
% nf := number of free cells in the map.
% f  := total number of sensing configurations per cell.
% 
% ----------------------------------------------------------------------------------------
% OUTPUTS:
% ----------------------------------------------------------------------------------------
% C         : [nf*f,1][binary] Include/exclude status of all sensing conf.
% 1 := selected, 0 := not selected. 
% C0        : [---][integer] IDs of zero conf using re-weighted l1 minimization.
% C1        : [---][integer] IDs of selected conf using re-weighted l1 minimization.
% Cu        : [---][integer] IDs of uncertain conf using re-weighted l1 minimization.
% Cu        : [---][integer] IDs of selected conf using combinatorial optimization.
% lBound    : [scalar][real] Lower bound of the solution.
% uBound    : [scalar][integer] Upper bound of the solution.
% logs_rwl1  : [k,5][] Logs for each (k) iterations of re-weightes l1
% minimization. Loged parameters are: iterat := iteration 
%                     number, epsilon := epsilon to find next weights,
%                     improvCu := difference of uncertain conf from the last 
%                     iterations, numCu := number of uncertain conf, and
%                     tItr := computation time for the iteration. 
% tSPP_rwl1 : [scalar][sec] computation time for re-weighted l1 minimization.
% tSPP_Comb : [scalar][sec] computation time for combinatorial optimization solution�.
% tSPP      : [scalar][sec] total computation time to solve AG problem.
% 
% ----------------------------------------------------------------------------------------
% INPUTS:
% ----------------------------------------------------------------------------------------
% V         : [nf,nf*f][binary] Visibility matrix. 1 := visibile, 0 := not
% visibile. rows are cells and colums are conf. 
% FoVfaces  : .num [][] number of FoVs for each cell.
% map       : [l,m][binary] l times m map. 1 := free cell, 0 := occupied cell.
% paramSPP  : [][] .numHistIR := number of last iterations to keep in
% history table of improvement rate for uncertain conf. 
%                  .maxIterat := maximum number of iterations for rwl1.
%                  .numCu     := maximum number of uncertain conf to
%                  terminate rwl1 search. 
%                  .Prog      := program script to be used. 'CVX' for CVX,
%                  'GRB' for gurobi, and 'MOT' for matlab 
%                                optimization toolbox.
% 
% ----------------------------------------------------------------------------------------
% NOTES/UPDATES:
% ----------------------------------------------------------------------------------------
% R5 2014-10-17: AGP -> SPP
% R4 2013-12-16: Evaluation of 'epsilon' is changed as epsilon = epsilon^(1+((iterat-1)*1e-1)).
% R3 2013-12-13: Evaluation of 'epsilon' is changed as epsilon = (epsilon/(iterat)^0.25).
% R2           : 
% R1           :
% R0           :
%  

% Change solver to Gurobi/MOSEK
% cvx_solver Gurobi
% cvx_solver

%% DIMENSIONS, CONSTANTS, INPUT MATRIX ETC.:

n       = size(V,1); %numel(find(map));   % number of unoccupied cells.
% f       = FoVfaces.num;       % rename for simplification.
c       = size(V,2); %n*f;                % total number of conf.
cutC    = 1e-2;               % cutoff/threshold value to select/reject a conf.

%% SPP based on RE-WEIGHTED L1 MINIMIZATION:

tSPPi        = tic;                          % timer for sensor placement problem.
tSPP_rwl1i   = tic;                          % timer for convex relaxation solution.
W            = ones(c,1);                    % initialize weight vector.
iterat       = 0;                            % iterations count.
cond         = 0;                            % condition to go for next iteration.
numCu_prev   = inf;                          % initialize num of uncertain conf.
HistIR       = inf(para_.numHistIR,1);    % history array of uncertain conf improvement rate for n iterations.
logs_rwl1    = zeros(para_.maxIterat,5);  % data logging for each iteration of important parameters.

CIterativeValues = [];

while ~cond
    
    tItri = tic;
    switch para_.Prog
        
        case 'CVX'
            %cvx_solver Gurobi
            cvx_begin %quiet
            
                % --- VARIABLES ---
                variable C(c) nonnegative

                % --- OBJECTIVE FUNCTION ---
                minimize (W'*C)

                subject to
                    % --- COVERAGE CONSTRAINT:
                    % Each unoccupied cell is observed by at least one
                    % sensing configuration. 
                    V*C >= 1
            cvx_end
            
        case 'GRB'
            try
                clear model;
                model.A          = V;
                model.obj        = W;
                model.rhs        = ones(n,1);
                model.sense      = '>';
                model.vtype      = 'C';
                model.modelsense = 'min';

                clear params;
                params.outputflag = 0;
                resultrwl1 = gurobi(model,params);

            catch gurobiError
                fprintf('Error reported\n');
                gurobiError
                
                fprintf('%s\n',gurobiError);
            end
            C = resultrwl1.x;
            
        case 'MOT'
            f        = W;
            A        = -V;
            b        = -1*ones(n,1);
            Aeq      = [];
            beq      = [];
            lb       = zeros(c,1);
            options  = optimset('LargeScale','on');
            C        = linprog(f,A,b,Aeq,beq,lb,[],[],options);
    end
    
        
    % ----------------- Update Weights ---------------------------------------------------
    iterat  = iterat+1;                      % Iteration count.
    epsilon = 1/(exp(1)-1);                  % threshold/epsilon.
    epsilon = epsilon^(1+((iterat-1)*1e-1)); % Update threshold/epsilon for the iteration.
    W       = epsilon./((C+epsilon));        % update weight vector.
    % ------------------------------------------------------------------------------------
      
    % ------------ Condition for the next iteration --------------------------------------
    tItr       = 1e-4*(round(toc(tItri)*1e4));  % time taken for the current iteration.
    Cu         = find(C>=cutC & C<=1-cutC);     % Uncertain Configurations.
    numCu      = numel(Cu);                     % num of uncertain conf.
    improvCu   = numCu_prev-numCu;              % difference of uncertain conf from the previous iteration.
    numCu_prev = numCu;                         % num of (previous) uncertain conf for the next iteration.
    HistIR     = circshift(HistIR,1);           % circle rotation of history of improved uncertain conf for the last n iterations.
    HistIR(1)  = improvCu;                      % update current improvement history of uncertain conf list.
    cond       = (iterat>=para_.maxIterat|numCu<=para_.numCu|~any(HistIR));  % condition for the next iteration.
    %cond       = (iterat>=para_.maxIterat|numCu<=5|~any(HistIR))  % condition for the next iteration.
    %iterat>=para_.maxIterat
    %numCu<=2
    %~any(HistIR)
    
    %figure; stem(C)
    %CIterativeValues(:,iterat) = C;
    % ------------------------------------------------------------------------------------
    %disp(['Iteration: ',num2str(iterat),...
    %    ', Num of uncertain conf.: ',num2str(numCu),...
    %    ', Improvement rate: ',num2str(improvCu),...
    %    ', Time: ',num2str(tItr)]);
    disp(['it: ',num2str(iterat),...
        ', uncertain confs.: ',num2str(numCu),...
        ', improvement rate: ',num2str(improvCu)]);
    
    % -------------------- Lower Bound ---------------------------------------------------
    if (iterat == 1)
        lBound = sum(C);
    end
    % ------------------------------------------------------------------------------------
    
    % --------------------------------------- Data Logging -------------------------------
    logs_rwl1(iterat,:) = [iterat,epsilon,improvCu,numCu,tItr];
    
    %pause
end

logs_rwl1(iterat+1:end,:)  = [];                                  % remove unused rows.
tSPP_rwl1                  = 1e-4*(round(toc(tSPP_rwl1i)*1e4));   % computation time for convex relaxation solution.
C0                         = find(C <= cutC );                    % not selected/zero conf.
C1                         = find(C >= 1-cutC );                  % selected/one conf.
C                          = (C >= 1-cutC );                      % selected/one conf.
% ----------------------------------------------------------------------------------------
disp(['tSPP_rwl1: ',num2str(tSPP_rwl1)]);


%% COMBINATORIAL OPTIMIZATION FOR UNCERTAIN CONFIGURATIONS:

tSPP_Combi  = tic;               % initialize timer for comb solution.
Vu          = V(:,Cu);           % visibility matrix for uncertain conf.
V1          = sum(V(:,C1),2);    % visibility status for the selected/one conf from the previous solution.
W           = ones(numel(Cu),1); % unity weights.

switch para_.Prog
    case 'CVX'
        cvx_begin %quiet
            % --- VARIABLES ---
            variable xC(numel(Cu)) binary
            
            % --- OBJECTIVE FUNCTION ---
            minimize (W'*xC)

            subject to
                % --- COVERAGE CONSTRAINT:
                % Each unoccupied cell is observed by at least one sensing configuration.
                Vu*xC + V1 >= 1
        cvx_end
        xC = round(xC);
        
    case 'GRB'
        try
            clear model;
            model.A          = sparse(Vu);
            model.obj        = W;
            model.rhs        = ones(n,1)-V1;
            model.sense      = '>';
            model.vtype      = 'B';
            model.modelsense = 'min';

            clear params;
            params.outputflag = 0;
            resultL0 = gurobi(model, params);

        catch gurobiError
            fprintf('Error reported\n');
        end
        xC = round(resultL0.x);
        
    case 'MOT'
        f         = W;
        A         = -Vu;
        b         = -1*(ones(n,1)-V1);
        Aeq       = [];
        beq       = [];
        options   = optimset('LargeScale','on');
        xC        = bintprog(f,A,b,Aeq,beq,[],options);
        xC        = round(xC);
end

combC     = Cu(~~xC);                           % IDs of selected conf out of uncertain conf using combinatorial opt.
C(combC)  = 1;                                  % update C with combC.
uBound    = sum(xC)+numel(C1);                  % upperbound is sum of rwl1 and Comb. optimization results.
tSPP_Comb = 1e-4*(round(toc(tSPP_Combi)*1e4));  % computation time for comb optimization.
tSPP      = 1e-4*(round(toc(tSPPi)*1e4));       % total computation time.


% CIterativeValues(:,iterat+1) = C
% figure; stem(C)

save('CIterativeValues.mat','CIterativeValues');

% ----------------------------------------------------------------------------------------
disp(['Num of zero conf.     (rwl1): ',num2str(numel(C0))])
disp(['Num of selected conf. (rwl1): ',num2str(numel(C1))])
disp(['Num of selected conf. (Comb): ',num2str(numel(combC))])
disp(['Total conf.                 : ',num2str(numel(find(C)))]);
disp(['lower bound.                : ',num2str(lBound)]);

disp(['tSPP_Comb: ',num2str(tSPP_Comb),' sec']);
disp(['tSPP     : ',num2str(tSPP),' sec']);

end


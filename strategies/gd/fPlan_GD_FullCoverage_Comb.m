function [ C,statCVX,tSPP ] = fPlan_GD_FullCoverage_Comb( V,FoVfaces,map,paramSPP )
%fAGP_CombR1_CVX is an Art Gallery Problem solution based on combinatorial optimization.
% Date: 2013-12-20, Rev 1
% 
% n  := total number of cells in the map.
% nf := number of free cells in the map.
% f  := total number of sensing configurations per cell.
% 
% --------------------------------------------------------------------------------------------------------------------------------
% OUTPUTS:
% --------------------------------------------------------------------------------------------------------------------------------
% C    : [nf*f,1][binary] Include/exclude status of all sensing conf. 1 := selected, 0 := not selected.
% tSPP : [scalar][sec] total computation time to solve AG problem.
% 
% --------------------------------------------------------------------------------------------------------------------------------
% INPUTS:
% --------------------------------------------------------------------------------------------------------------------------------
% V        : [nf,nf*f][binary] Visibility matrix. 1 := visibile, 0 := not visibile. rows are cells and colums are conf.
% FoVfaces : .num [][] number of FoVs for each cell.
% map      : [l,m][binary] l times m map. 1 := free cell, 0 := occupied cell.
% paramSPP : [][] .Prog := program script to be used. 'CVX' for CVX, 'GRB' for gurobi, and 'MOT' for matlab optimization
%                          toolbox.
% 
% --------------------------------------------------------------------------------------------------------------------------------
% NOTES/UPDATES:
% --------------------------------------------------------------------------------------------------------------------------------
%  

% Change solver to Gurobi/MOSEK
% cvx_solver Gurobi
% 
% cvx_solver

% Time Limit - Gurobi
% cvx_solver_settings( 'TimeLimit', 3600 );

%% DIMENSIONS, CONSTANTS, INPUT MATRIX ETC.:

n = numel(find(map));   % number of unoccupied cells.
f = FoVfaces.num;       % rename for simplification.
c = n*f;                % total number of conf.

statCVX = [];

%% COMBINATORIAL OPTIMIZATION SOLUTION:

tSPPi = tic;         % timer for sensor placement problem.
W     = ones(c,1);   % initialize weight vector.

switch paramSPP.Prog
    
    case 'CVX'
        cvx_solver Gurobi
        cvx_solver
        % Time Limit - Gurobi
        cvx_solver_settings( 'TimeLimit', 3600 );
        
        cvx_begin quiet
            % --- VARIABLES ---
            variable C(c) binary
            
            % --- OBJECTIVE FUNCTION ---
            minimize (W'*C)

            subject to
                % --- COVERAGE CONSTRAINT:
                % Each unoccupied cell is observed by at least one sensing configuration.
                V*C >= 1
        cvx_end
        
        statCVX.cvx_optval = cvx_optval;
        statCVX.cvx_status = cvx_status;
        C = round(C);
        
    case 'GRB'
        try
            clear model;
            model.A          = V;
            model.obj        = W;
            model.rhs        = ones(n,1);
            model.sense      = '>';
            model.vtype      = 'B';
            model.modelsense = 'min';

            clear params;
            params.outputflag = 0;
            result = gurobi(model,params);

        catch gurobiError
            fprintf('Error reported\n');
        end
        C = round(result.x);
        
    case 'MOT'
        f         = W;
        A         = -V;
        b         = -1*ones(n,1);
        Aeq       = [];
        beq       = [];
        options   = optimset('LargeScale','on');
        C         = bintprog(f,A,b,Aeq,beq,[],options);
        C         = round(C);
end

% total computation time.
tSPP = 1e-4*(round(toc(tSPPi)*1e4)); 


% --------------------------------------------------------------------------------------------------------------------------------
disp(['Total conf.: ',num2str(numel(find(C)))]);
disp(['tSPP       : ',num2str(tSPP),' sec']);

end

function [mean_map_vector] = CalculateGDM(A,ppmm,lambda)

% --- original function
%[mean_map_vector, variance_map_vector, var_process_noise] = CalculateGDM(A,ppmm,lambda)

    %MEAN MAP
    %[mean_map_vector, lm] = SolveQP(A,ppmm,lambda);
    [mean_map_vector, ~] = SolveQP(A,ppmm,lambda); % -- tempo
    
    %invalid_variance_entry = logical(lm > 1e-10); % --- temporarily commited ---
    
    %PROCESS NOISE
    %[num_measurements, num_cells] = size(A); % --- temporarily commited ---
    %res = ppmm - A*mean_map_vector; % --- temporarily commited ---
    %k = rank(full(A));
    %num_measurements
    %num_cells
    %var_process_noise = res'*res/(num_measurements - num_cells); % -- tempo
    
    %VARIANCE MAP
    % --- temporarily commited ---
    %inv_des = inv(A'*A + lambda*eye(num_cells));
    %variance_map_vector = var_process_noise./diag(inv_des);
    
end

function [solution, lagrange_multipliers] = SolveQP(A,ppmm,lambda)
    %build the quadratic term
    H = A'*A + lambda*eye(size(A,2));
    %build the linear term
    q = full(-2*ppmm'*A);

    clear model;
    model.A = sparse(eye(size(A,2)));
    model.obj = q;
    model.sense = '>';
    model.rhs = zeros(size(A,2),1);
    %model.lb = zeros(size(A,2),1);
    model.Q = sparse(H);

    clear params;
    params.method = -1;
    params.OutputFlag = 0;

    result = gurobi(model, params);
    solution = result.x;
    lagrange_multipliers = result.pi;
end
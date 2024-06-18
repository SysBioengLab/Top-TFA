function feasibilityProblem = looplessCheck(model,flag)
% Extract relevant information
Nint  = null(full(model.S(:,model.internal)),'r')';
[m,n] = size(Nint);
feasibilityProblem.obj = zeros(n,1);
feasibilityProblem.modelsense = 'max';

% Model definition (stays the same)
feasibilityProblem.A = sparse(Nint); % Loop condition

% Sign assignation (stays the same)
feasibilityProblem.sense = '=';

% Assignation of variable types (stays the same)
feasibilityProblem.vtype = 'C';

% Right-hand side definition
feasibilityProblem.rhs = zeros(m,1);

% Write problem for debugging (optional)
if nargin > 1 && flag
    gurobi_write(feasibilityProblem,'check.lp');
end
function [sample,warmupPoints] = warmupFcTMFA(sample,vTol)
% Generates a set of samples. Based on warmupLooplessACHRB.

% Initialize main variables
if (nargin<2); vTol = 1e-8; end

% Define number of warmup points
warmupPoints = zeros(sample.numVars,2*sample.numVars);
LPproblem.A       = sample.S;
LPproblem.b       = sample.b;
LPproblem.lb      = sample.lb;
LPproblem.ub      = sample.ub;
LPproblem.c       = zeros(size(LPproblem.ub));
LPproblem.csense  = repmat('E',size(LPproblem.A,1),1);

% Define parameters for the warmup process.
lbRef = LPproblem.lb;
ubRef = LPproblem.ub;
alpha = 0.95;

% Main loop
for counter = 1:sample.numVars

    % Change objective
    LPproblem.c(counter)  = 1;
    
    % Solve maximization problem
    LPproblem.osense      = -1;    
    LPproblem.ub(counter) = lbRef(counter) + alpha*(ubRef(counter)-lbRef(counter));
    solMax                = solveCobraLP(LPproblem);
    solMax.full(abs(solMax.full)<vTol) = 0;
    LPproblem.ub(counter) = ubRef(counter);
    warmupPoints(:,2*counter-1) = solMax.full;

    % Solve minimization problem
    LPproblem.osense      = 1;
    LPproblem.lb(counter) = lbRef(counter) + (1-alpha)*(ubRef(counter)-lbRef(counter));
    solMin                = solveCobraLP(LPproblem);
    solMin.full(abs(solMin.full)<vTol) = 0;
    LPproblem.lb(counter) = lbRef(counter);
    warmupPoints(:,2*counter) = solMin.full;
    
    % Reset original obj. values
    LPproblem.c(counter) = 0;

end

end
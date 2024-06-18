function modGurobi = looplessStructureMILP(model,flag)
% Builds loopless structure problem
% Inputs:  model structure, flag (optional)
% Outputs: loopless model structure
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2015 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[l,n] = size(model.S);

% Useful information for later
model = parseInternalRxns(model,n);

% Reorganize the stoichiometric matrix as [ìnternalRxns|exchangeRxns]
modGurobi.numRxns  = n;
modGurobi.internal = model.internal;
modGurobi.exchange = model.exchange;
modGurobi.S        = model.S;
modGurobi.rxns     = model.rxns;
modGurobi.rxnNames = model.rxnNames;
modGurobi.mets     = model.mets;
modGurobi.metNames = model.metNames;
modGurobi.description = model.description;

% Calculate null basis of Sint
try
     modGurobi.Nint = (fastSNP(model,'gurobi'))';
     if isempty(modGurobi.Nint)
         modGurobi.Nint = null(full(model.S(:,modGurobi.internal)),'r')';
     end
catch
    modGurobi.Nint = null(full(model.S(:,modGurobi.internal)),'r')';
end
[p,m] = size(modGurobi.Nint);

% Determine constrained internal reactions
M = eye(n);
M(modGurobi.exchange,:) = [];

% Define binary constant, model objective and sense
K = 1e4;
modGurobi.obj        = zeros(n+2*m,1);
modGurobi.modelsense = 'max';

% Constraints definition
modGurobi.A = sparse([model.S,zeros(l,2*m);...
    zeros(p,n),modGurobi.Nint,zeros(p,m);...
    zeros(m,n),eye(m),(1+K)*eye(m);...
    zeros(m,n),-eye(m),-(1+K)*eye(m);...
    M,zeros(m),-diag(model.ub(modGurobi.internal));...
    -M,zeros(m),-diag(model.lb(modGurobi.internal))]);

% RHS model
modGurobi.rhs = [zeros(l+p,1);K*ones(m,1);-ones(m,1);zeros(m,1);-model.lb(modGurobi.internal)];

% Sign assignation
modGurobi.sense = blanks(l+p+4*m);
for i = 1:l+p+4*m
    if i <= l+p
        modGurobi.sense(i) = '=';
    else
        modGurobi.sense(i) = '<';
    end
end

% Assignation of variable types
modGurobi.vtype = blanks(n+2*m);
for i = 1:n+2*m
    if i <= n+m
        modGurobi.vtype(i) = 'C';
    else
        modGurobi.vtype(i) = 'B';
    end
end

% Bounds assignation
modGurobi.lb = [model.lb;-K*ones(m,1);zeros(m,1)];
modGurobi.ub = [model.ub;K*ones(m,1);ones(m,1)];

% Write problem for debugging
if flag
    gurobi_write(modGurobi,'looplessMILP.lp');
end
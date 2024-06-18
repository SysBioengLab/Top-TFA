function modGurobi = looplessStructureLP(model)
% Builds loopless structure problem
% Inputs:  model structure
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
modGurobi.Nint = null(full(modGurobi.S(:,modGurobi.internal)),'r')';
[p,m] = size(modGurobi.Nint);

% Define model objective and sense
K = 1e4;
modGurobi.obj        = zeros(n+m,1);
modGurobi.modelsense = 'min';

% Constraints definition
modGurobi.A = sparse([model.S,zeros(l,m);zeros(p,n),modGurobi.Nint]);

% RHS model
modGurobi.rhs = zeros(l+p,1);

% Sign assignation
modGurobi.sense = '=';

% Assignation of variable types
modGurobi.vtype = 'C';

% Bounds assignation
modGurobi.lb  = [model.lb;-K*ones(m,1)];
modGurobi.ub  = [model.ub;K*ones(m,1)];
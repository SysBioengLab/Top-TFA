function modelLP = structureLP(model)
% Builds loopless structure problem
% Inputs:  model structure
% Outputs: loopless model structure
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2015 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[l,n] = size(model.S);

% Useful information for later
model = parseInternalRxns(model,n);

% Reorganize the stoichiometric matrix as [ìnternalRxns|exchangeRxns]
modelLP.numRxns  = n;
modelLP.internal = model.internal;
modelLP.exchange = model.exchange;
modelLP.S        = model.S;
modelLP.rxns     = model.rxns;
modelLP.rxnNames = model.rxnNames;
modelLP.mets     = model.mets;
modelLP.metNames = model.metNames;
modelLP.description = model.description;

% Define model objective and sense
modelLP.obj        = zeros(n,1);
modelLP.modelsense = 'min';

% Constraints definition
modelLP.A = sparse([model.S]);

% RHS model
modelLP.rhs = zeros(l,1);

% Sign assignation
modelLP.sense = '=';

% Assignation of variable types
modelLP.vtype = 'C';

% Bounds assignation
modelLP.lb  = model.lb;
modelLP.ub  = model.ub;
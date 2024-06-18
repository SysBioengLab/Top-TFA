function model = looplessFVA(model)
% Performs loopless-FVA
% Inputs:  model structure
% Outputs: model structure with new bounds and directionalities
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2015 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.outputflag = 0;
tol = 1e-9;

% Define and run FVA problem
f = model.obj;
FVA = zeros(model.numRxns,2);

% FVA cycle
for i = 1:model.numRxns
    model.obj = f;

    % Mimization problem
    model.obj(i) = -1;
    sol = gurobi(model,params);
    FVA(i,1) = sol.x(i);
    
    % Maximization problem;
    model.obj(i) = 1;
    sol = gurobi(model,params);
    FVA(i,2) = sol.x(i);
end

% Find zero reactions
zeroRxns = find((abs(FVA(:,1))<tol)+(abs(FVA(:,2))<tol)==2);

if isempty(zeroRxns)   
    % Reset objective function
    model.obj = f;
    
    % Redefine model boundaries
    model.lb(1:model.numRxns) = FVA(:,1).*(abs(FVA(:,1))>tol);
    model.ub(1:model.numRxns) = FVA(:,2).*(abs(FVA(:,2))>tol);
    
    % Re-format the problem
    model = looplessStructureMILP(model,0);
else
    % Remove zero rxns    
    modelTemp.S = model.S;    
    modelTemp.S(:,zeroRxns') = [];
    modelTemp.c = zeros(size(modelTemp.S,2),1);
        
    % Find orphan metabolites
    orphanMets = find(sum(full(abs(modelTemp.S)),2)==0);
    modelTemp.S(orphanMets,:) = [];
        
    % Re-define other quantities
    modelTemp.b  = zeros(size(modelTemp.S,1),1);
    modelTemp.lb = FVA(:,1).*(abs(FVA(:,1))>tol);
    modelTemp.lb(zeroRxns) = [];
    modelTemp.ub = FVA(:,2).*(abs(FVA(:,2))>tol);
    modelTemp.ub(zeroRxns) = [];
    modelTemp.rxns = model.rxns;
    modelTemp.rxns(zeroRxns) = [];
    modelTemp.rxnNames = model.rxnNames;
    modelTemp.rxnNames(zeroRxns) = [];
    modelTemp.mets = model.mets;
    modelTemp.mets(orphanMets) = [];
    modelTemp.metNames = model.metNames;
    modelTemp.metNames(orphanMets) = [];
    modelTemp.description = model.description;
    
    % Re-format the problem
    model = looplessStructureMILP(modelTemp,0);
end

% Assign reversibilities to thermodynamically allowable rxns
model.rev = (model.lb(1:model.numRxns)<-tol).*(model.ub(1:model.numRxns)>tol);
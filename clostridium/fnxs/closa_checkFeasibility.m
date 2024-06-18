function [feasible,bounds,FVAbounds,dGrBounds] = closa_checkFeasibility(model,tfm)
% Fix reversible rxns directions by tfm and run MVA. Checks if the
% topologies are feasible and return the bounds given by MVA.

nTopologies = size(tfm,1);
feasible = false(nTopologies,1);
nVars  = size(model.lb,1) + size(model.lbc,1) + size(model.dGftLB,1); % With proton concentration
bounds = zeros(nVars,2,nTopologies);
FVAbounds = zeros(size(model.lb,1),2,nTopologies);
flags  = zeros(nVars,2,nTopologies);
exchange  = findExRxns(model);
dGrBounds = zeros(sum(~exchange),2,nTopologies);

% Main loop
for i=1:nTopologies
    fixedModel = fixDrxns(model,tfm(i,:));
    newModel = closa_newMVA(fixedModel);
    bounds(:,:,i) = [newModel.tlb newModel.tub];
    flags(:,:,i)  = [newModel.flagLB newModel.flagUB];
    dGrBounds(:,:,i) = newModel.dGrBounds;
    if any(flags(:,:,i)<0,'all')
        feasible(i) = 0;
    else
        feasible(i) = 1;
    end
    fixedModel.lb(fixedModel.lb == 0) =  0.001;
    fixedModel.ub(fixedModel.ub == 0) = -0.001;
    [FVAbounds(:,1,i),FVAbounds(:,2,i)] = fluxVariability(fixedModel,0);
end
end
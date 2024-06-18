function [feasible,bounds,FVAbounds,fluxSolutions,dGrBounds,meanTime] = checkFeasibility(model,tModel,tfm,relax)
% Fix reversible rxns directions by tfm and run MVA. Checks if the
% topologies are feasible and return the bounds given by MVA.

nTopologies = size(tfm,1);
feasible = false(nTopologies,1);
% nVars  = size(tModel.lb,1) + size(tModel.lbc,1) -2; % Without proton concentration
nVars     = size(tModel.lb,1) + size(tModel.lbc,1) + size(tModel.dGftLB,1); % With proton concentration
bounds    = zeros(nVars,2,nTopologies);
FVAbounds = zeros(size(tModel.lb,1),2,nTopologies);
flags     = zeros(nVars,2,nTopologies);
optTimes  = zeros(1,nTopologies);
exchange  = findExRxns(tModel);
fluxSolutions = zeros(nVars,2*size(tModel.lb,1),nTopologies);
dGrBounds = zeros(sum(~exchange),2,nTopologies);

% Main loop
for i=1:nTopologies
    fixedModel = fixDrxns(model,tfm(i,:));
    fixedTherModel = fixDrxns(tModel,tfm(i,:));
    optStart = tic;
    newModel = newMVA(fixedModel,fixedTherModel,relax);
    optEnd = toc(optStart);
    optTimes(i) = optEnd;
    bounds(:,:,i) = [newModel.tlb newModel.tub];
    flags(:,:,i)  = [newModel.flagLB newModel.flagUB];
    fluxSolutions(:,:,i) = newModel.X(:,1:2*size(tModel.lb,1));
    dGrBounds(:,:,i) = newModel.dGrBounds;
    if any(flags(:,:,i)<0,'all')
        feasible(i) = 0;
    else
        feasible(i) = 1;
    end
    if relax == 1 % Arbitrary: To avoid doing FVA 2 times
        fixedModel.lb(fixedModel.lb == 0) =  0.001;
        fixedModel.ub(fixedModel.ub == 0) = -0.001;
        [FVAbounds(:,1,i),FVAbounds(:,2,i)] = fluxVariability(fixedModel,0);
    end
end

meanTime = mean(optTimes(feasible));

end
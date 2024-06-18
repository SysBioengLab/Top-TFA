function [fixedModel] = fixDrxns(model,dirVec)
% Fixes directions of reversible reactions given by dirVec.
% Raises a Warning if the rxn does not match the bounds

completeDir = zeros(size(model.rxns));
model.rev = (model.lb<-1e-8 & model.ub>1e-8);
completeDir(model.rev) = dirVec;            % Must have equal size
fixedModel = model;
fixedModel.ub(completeDir==-1) = 0;
fixedModel.lb(completeDir==1)  = 0;

end
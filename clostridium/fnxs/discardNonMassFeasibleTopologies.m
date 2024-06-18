function [tfm, FVAbounds] = discardNonMassFeasibleTopologies(model,tfm)

nTopologies = size(tfm,1);
massFeasible = false(nTopologies,1);
FVAbounds = zeros(size(model.lb,1),2,nTopologies);

for i=1:nTopologies
    fixedModel = fixDrxns(model,tfm(i,:));
    try
        [FVAbounds(:,1,i),FVAbounds(:,2,i)] = fluxVariability(fixedModel,1);
        massFeasible(i) = 1;
    catch
        continue
    end
end

tfm(~massFeasible,:) = [];
FVAbounds(:,:,~massFeasible) = [];
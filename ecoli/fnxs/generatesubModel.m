function [subModel,subTherModel] = generatesubModel(model,rxnList,unbalancedMets)
% Generates a copy of model but only includes the rxns in rxnList. It also
% removes the unbalancedMets from the S matrix. Compress the model and
% checks if FVA is working properly.

% Generate vector of metabolites to remove
rmvdMets = false(size(model.mets));
for i=1:length(unbalancedMets)
    rmvdMets = rmvdMets + ismember(model.mets,unbalancedMets{i}); % Adds exact match
end
rmvdMets = logical(rmvdMets);

% Generate vector of reactions to remove
rmvdRxns = false(size(model.rxns));
for i=1:length(rxnList)
    rmvdRxns = rmvdRxns + ismember(model.rxns,rxnList{i}); % Adds exact match
end
rmvdRxns = ~rmvdRxns;

% Make copy of model
subModel = model;
subModel.description = "mass";

% Remove reactions from other fields
subModel.rxns(rmvdRxns)           = [];
subModel.rxnNames(rmvdRxns)       = [];
subModel.subSystems(rmvdRxns)     = [];
subModel.S(:,rmvdRxns)            = [];
subModel.T(:,rmvdRxns)            = [];
subModel.charVec(rmvdRxns)        = [];
subModel.lb(rmvdRxns)             = [];
subModel.ub(rmvdRxns)             = [];
subModel.c(rmvdRxns)              = [];
subModel.rev(rmvdRxns)            = [];
subModel.rxnGeneMat(rmvdRxns,:)   = [];
subModel.grRules(rmvdRxns)        = [];

% If there are unbalanced mets, generate subTherModel to save thermodynamic
% information
if sum(rmvdMets) > 0
    subTherModel = subModel;
    subTherModel.description = "thermo";
    subTherModel.unbalancedMets = rmvdMets;
    subTherModel = compressModel(subTherModel);
end

% Remove unbalanced mets from fields
subModel.mets(rmvdMets)        = [];
subModel.metNames(rmvdMets)    = [];
subModel.metFormulas(rmvdMets) = [];
subModel.metCharge(rmvdMets)   = [];
subModel.S(rmvdMets,:)         = [];
subModel.T(rmvdMets,:)         = [];
subModel.b(rmvdMets)           = [];
subModel.lbc(rmvdMets)         = [];
subModel.ubc(rmvdMets)         = [];
subModel.dGft_mu(rmvdMets)     = [];
subModel.dGft_cov(rmvdMets,:)  = [];
subModel.dGft_cov(:,rmvdMets)  = [];
subModel.dGft_std(rmvdMets)    = [];
subModel.dGftLB(rmvdMets)      = [];
subModel.dGftUB(rmvdMets)      = [];
subModel.NH(rmvdMets)          = [];

% FVA run and compression of model
[subModel.lb, subModel.ub] = fluxVariability(subModel,5);
subModel = compressModel(subModel);
subModel.rev(subModel.lb<0 & subModel.ub>0) = 1;     % Update rev rxns
subModel.rev(~(subModel.lb<0 & subModel.ub>0)) = 0;  % Update rev rxns

% Correction of subTherModel bounds and rxns
if sum(rmvdMets) > 0
    subTherModel = sincronizeModels(subModel,subTherModel);
    subTherModel = compressModel(subTherModel);
    subTherModel.lb = subModel.lb; subTherModel.ub = subModel.ub;
    subTherModel.rev = subModel.rev;  % Update rev rxns
else
    subTherModel = subModel;
end

end
function [myTherModel] = closa_rmvLoopyRxns(myTherModel,model)

br = myTherModel.rxns;
br(ismember(br,model.rxns)) = [];
bR = ismember(myTherModel.rxns,br);

if sum(bR)>0
    disp("The following loopy reactions have been removed:")
    disp(myTherModel.rxns(bR))
    bR = ismember(myTherModel.rxns,br);
    myTherModel.rxnGeneMat(bR,:) = [];
    myTherModel.grRules(bR)      = [];
    myTherModel.rxns(bR)         = [];
    myTherModel.rxnNames(bR)     = [];
    myTherModel.subSystems(bR)   = [];
    myTherModel.S(:,bR)          = [];
    myTherModel.T(:,bR)          = [];
    myTherModel.charVec(bR)      = [];
    myTherModel.lb(bR)           = [];
    myTherModel.ub(bR)           = [];
    myTherModel.c(bR)            = [];
    myTherModel.rev(bR)          = [];
    
    % Run FVA to update boundsfixedModel
    [myTherModel.lb, myTherModel.ub] = fluxVariability(myTherModel,5);
    myTherModel = compressModel(myTherModel);
    myTherModel.rev(myTherModel.lb<0 & myTherModel.ub>0) = 1;     % Update rev rxns
    myTherModel.rev(~(myTherModel.lb<0 & myTherModel.ub>0)) = 0;  % Update rev rxns
else
    disp("No loopy rxns found")
end

end
function [myModel,myTherModel] = rmvLoopyRxns(myModel,myTherModel,model)

br = myModel.rxns;
br(ismember(br,model.rxns)) = [];
bR = ismember(myModel.rxns,br);

if sum(bR)>0
    disp("The following loopy reactions have been removed:")
    disp(myModel.rxns(bR))
end

myModel.rxnGeneMat(bR,:) = [];
myModel.grRules(bR)      = [];
myModel.rxns(bR)         = [];
myModel.rxnNames(bR)     = [];
myModel.subSystems(bR)   = [];
myModel.S(:,bR)          = [];
myModel.T(:,bR)          = [];
myModel.charVec(bR)      = [];
myModel.lb(bR)           = [];
myModel.ub(bR)           = [];
myModel.c(bR)            = [];
myModel.rev(bR)          = [];

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
[myModel.lb, myModel.ub] = fluxVariability(myModel,5);
myModel = compressModel(myModel);
myModel.rev(myModel.lb<0 & myModel.ub>0) = 1;     % Update rev rxns
myModel.rev(~(myModel.lb<0 & myModel.ub>0)) = 0;  % Update rev rxns
% Update TherModel
myTherModel = sincronizeModels(myModel,myTherModel);
myTherModel = compressModel(myTherModel);
myTherModel.lb = myModel.lb; myTherModel.ub = myModel.ub;
myTherModel.rev = myModel.rev;  % Update rev rxns
% myTherModel.lb = myModel.lb; myTherModel.ub = myModel.ub;
% myTherModel.rev = myModel.rev;

end
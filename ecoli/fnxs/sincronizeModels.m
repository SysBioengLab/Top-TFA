function [myTherModel] = sincronizeModels(myModel,myTherModel)

br = myTherModel.rxns;
br(ismember(br,myModel.rxns)) = [];
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

end
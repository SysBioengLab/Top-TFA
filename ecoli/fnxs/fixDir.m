function newModel = fixDir(model,sol1)
% Fixes the direction of reversible rxns.
% Temporalmente se incluyen las reacciones de exchange
newModel = model;

j = 1;
for i=newModel.rev
    if sol1(j) == 1
        newModel.lb(i) = 0;
    else
        newModel.ub(i) = 0;
    end
    j = j+1;
end
end
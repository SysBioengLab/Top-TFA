function [model1, model2] = buildComplementaryModels(model)
% Creates two models
% Fields S, lb and ub are required.
model = removeBlockedRxns(model);
model1.S = model.S;
model2.S = model.S;
lb = model.lb;
ub = model.ub;
model1.lb = lb;
model2.lb = lb;
model1.ub = ub;
model2.ub = ub;

rev = zeros(size(lb));
rev(sign(lb) .* sign(ub) == -1) = 1;
rev = logical(rev);
model1.lb(rev) = 0;
model2.ub(rev) = 0;
model1.rev = rev;
model2.rev = rev;

end
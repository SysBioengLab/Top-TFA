function [minFlux,maxFlux] = FVA(model)
% Uses linprog to run Flux Variability Analysis

if ~isfield(model,'FVAoptions')
    model.FVAoptions = optimoptions('linprog','Display','none');
end

% FVA Parameters
nr_mets = size(model.S,1);
nr_rxns = size(model.S,2);
Aineq = [];
bineq = [];
Aeq = model.S;
beq = zeros(nr_mets,1);
lb = model.lb;
ub = model.ub;
minFlux = zeros(nr_rxns,1);
maxFlux = zeros(nr_rxns,1);
FVAoptions = model.FVAoptions;

%% FVA run
for i=1:1:nr_rxns
    f = zeros(1,nr_rxns);
    f(i) = 1;
    x = linprog(f,Aineq,bineq,Aeq,beq,lb,ub,FVAoptions);
    minFlux(i) = x(i);
    f(i) = -1;
    x = linprog(f,Aineq,bineq,Aeq,beq,lb,ub,FVAoptions);
    maxFlux(i) = x(i);
end

end
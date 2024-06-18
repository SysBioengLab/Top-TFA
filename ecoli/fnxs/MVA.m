function [model] = MVA(model)
% Performs Metabolite Variability Analysis. Updates model bounds and adds
% thermodinamic restrictions. Returns the updated model and 2*(n°
% variables) feasible points.

% Build Thermodynamic Model
model = buildTherModels(model);
% Run MVA
MVAoptions = optimoptions('intlinprog','Display','none');
nVar       = size(model.tlb,1);
nOptvar    = size(model.S,1)+size(model.S,2);
X          = zeros(nVar,2*nOptvar);           % Hacer + eficiente (no minimizar para las reversibles)
flagLB     = zeros(nOptvar,1);
flagUB     = zeros(nOptvar,1);
newLB      = model.tlb;
newUB      = model.tub;
f          = zeros(1,nVar);

for i=1:1:nOptvar
    f(i) = 1;
    [x, ~, flag] = intlinprog(f,model.intcon,model.Aineq,model.bineq, ...
        model.Aeq,model.beq,model.tlb,model.tub,[],MVAoptions);
    flagLB(i) = flag;
    if flag > 0
        newLB(i) = x(i);
        X(:,2*i-1) = x;
    end
    f(i) = -1;
    [x, ~, flag] = intlinprog(f,model.intcon,model.Aineq,model.bineq, ...
        model.Aeq,model.beq,model.tlb,model.tub,[],MVAoptions);
    flagUB(i) = flag;
    if flag > 0
        newUB(i) = x(i);
        X(:,2*i) = x;
    end
    f(i) = 0;
end
disp("MVA Complete!")

% Update model
model.tlb    = newLB;
model.tub    = newUB;
model.lbc    = model.tlb(size(model.S,2)+1:nOptvar);
model.ubc    = model.tub(size(model.S,2)+1:nOptvar);
model.trxns  = model.rxns;
model.flagLB = flagLB;
model.flagUB = flagUB;
model.X      = X;
model.tS     = model.S;
model.S      = model.oldS;
model.lb     = model.oldLB;
model.ub     = model.oldUB;
model.rxns   = model.oldrxns;

r = 0;
for i=1:size(model.oldLB,1)
    if ismember(i,model.rev)
        if model.tub(i+r)>0 && model.tub(i+r+1)>0 % rev
            model.ub(i) =  model.tub(i+r);
            model.lb(i) = -model.tub(i+r+1);
        elseif model.tub(i+r)==0 && model.tub(i+r+1)==0 % blocked
            model.lb(i) = 0;
            model.ub(i) = 0;
        elseif model.tub(i+r)==0 % backwards
            model.lb(i) = -model.tub(i+r+1);
            model.ub(i) = -model.tlb(i+r+1);
        else
            model.lb(i) = model.tlb(i+r);
            model.ub(i) = model.tub(i+r);
        end
        r = r+1;
    else
        model.lb(i) = model.tlb(i+r);
        model.ub(i) = model.tub(i+r);
    end
end
end
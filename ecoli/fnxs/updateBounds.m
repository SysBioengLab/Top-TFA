function model = updateBounds(model)
% Uses linprog to run Flux Variability Analysis and update model bounds.
% It updates 'rev' field.
% INPUTS:
%              model (structure):    (the following fields are required - others can be supplied)
%                                    * S  - `m x 1` Stoichiometric matrix
%                                    * lb - `n x 1` Lower bounds
%                                    * ub - `n x 1` Upper bounds
%                                    * b  - `m x 1` Solution vector of the system
% -------------------------------------------------------------------------

% FVA Parameters
nVars = size(model.S,2);
Aineq = [];
bineq = [];
Aeq = model.S;
beq = model.b;
lb  = model.lb;
ub  = model.ub;
f   = zeros(1,nVars);
newLB = zeros(nVars,1);
newUB = zeros(nVars,1);
FVAoptions = optimoptions('linprog','Display','none');

% FVA run
for j=1:nVars
    f(j) = 1;
    x = linprog(f,Aineq,bineq,Aeq,beq,lb,ub,FVAoptions);
    newLB(j) = x(j);
    f(j) = -1;
    x = linprog(f,Aineq,bineq,Aeq,beq,lb,ub,FVAoptions);
    newUB(j) = x(j);
    f(j) = 0;
end

% Update bounds
newLB(abs(newLB)<1e-8) = 0;
newUB(abs(newUB)<1e-8) = 0;
model.lb = newLB;
model.ub = newUB;
model.rev = (model.lb<-1e-8 & model.ub>1e-8);

end
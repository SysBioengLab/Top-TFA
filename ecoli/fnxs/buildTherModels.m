function model = buildTherModels(model)
% Add the necessary fields to solve MVA.

% Division of reversible fluxes
rev = false(size(model.S,2),1);
rev(model.lb<0 & model.ub>0) = 1;
newNrxns = size(model.S,2)+sum(rev);
newS = zeros(size(model.S,1),newNrxns);
newLB = zeros(newNrxns,1);
newUB = zeros(newNrxns,1);
newrxns = cell(size(model.rxns));
rev = find(rev);
model.rev = rev;

i = 1;
j = 1;
while i <= newNrxns
    newS(:,i) = model.S(:,j);
    if ismember(j,rev)
        newLB(i)  = 0;
        newUB(i)  = model.ub(j);
        newrxns(i) = model.rxns(j);
        i = i+1;
        newS(:,i) = -model.S(:,j);
        newLB(i)  = 0;
        newUB(i)  = -model.lb(j);
        newrxns(i) = strcat(model.rxns(j),'_r');
    else
        newLB(i) = model.lb(j);
        newUB(i) = model.ub(j);
        newrxns(i) = model.rxns(j);
    end
    i = i+1;
    j = j+1;
end

model.oldS  = model.S;
model.oldLB = model.lb;
model.oldUB = model.ub;
model.oldrxns = model.rxns;
model.S     = newS;
model.lb    = newLB;
model.ub    = newUB;
model.rxns  = newrxns;
model.exchange_rxns = findExcRxns(model);
model.thermo_rxns = ~model.exchange_rxns;
ext = logical(sum(abs(model.S) .* model.comp,1));
int = logical(sum(abs(model.S) .* ~model.comp,1));
model.transport_rxns = transpose(logical(int & ext)); % Rxns that involve intra and extra celular mets
charge = chargeMatrix(model);

% model       = parseInternalRxns(model); 

% Parameters definition
nMets    = size(model.S,1);
nRxns    = size(model.S,2);
nIntRxns = nRxns - sum(model.exchange_rxns);
slope    = ones(nIntRxns,1) .* 4.3478;                                     % slope that minimizes error
std_u    = model.dGf;                                                      % formation energies (KJ/mol)
R        = 8.3144626 * 1e-3;                                               % (KJ/mol K)
T        = 298.15;                                                         % (K) | 25 °C
F        = 96500;                                                          % Faraday constant
potential = -0.14;                                                         % -140 mV
% charge   = model.charge;                                                   % charge of mets
K        = 1e5;

% Equality constraints
model.Aeq = [model.S zeros(nMets,nMets+nIntRxns)];
model.beq = zeros(nMets,1);
% Inequality constraints
I1 = eye(nRxns);
I1 = I1(logical(model.thermo_rxns),:);                                     % thermo rxns
I       = eye(nIntRxns,nIntRxns);
vMaxint = model.ub(logical(model.thermo_rxns));
SintT   = transpose(model.S(:,logical(model.thermo_rxns)));
charge_term = F * potential * sum(transpose(charge(:,logical(model.thermo_rxns))) .* SintT,2);
cons1   = [I1 zeros(nIntRxns,nMets) -vMaxint.*I];
cons2   = [zeros(nIntRxns,nRxns) R*T*SintT K*I];
cons3   = [I1 R*T*(vMaxint./slope).*SintT K*I];
c1      = zeros(nIntRxns,1);
c2      = K*ones(nIntRxns,1)-SintT*std_u-charge_term-1e-5*ones(nIntRxns,1); % Strict inequality
c3      = K*ones(nIntRxns,1)-(vMaxint./slope).*SintT*std_u-(vMaxint./slope).*charge_term;
model.Aineq  = [cons1; cons2; cons3];
model.bineq  = [c1; c2; c3];
model.tlb    = [model.lb; model.lbc; zeros(nIntRxns,1)];
model.tub    = [model.ub; model.ubc; ones(nIntRxns,1)];
model.intcon = nRxns+nMets+1:size(model.tlb,1);

end
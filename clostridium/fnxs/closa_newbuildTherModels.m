function LPproblem = closa_newbuildTherModels(model)
% Creates the LPproblem that has to be solved in order to run MVA

tModel = model; % With closa they are the same. But recieve different treatment.
exchange = findExRxns(tModel);
tModel = closa_removeProtons(tModel); % Removing protons
fluxDir = ones(size(tModel.rxns));
fluxDir(tModel.ub<=1e-8) = -1;
Vmax = tModel.vmax;
Vmax(fluxDir==-1) = tModel.vmin(fluxDir==-1);
Vmax(exchange) = [];
fluxDir(exchange) = [];
St = transpose(tModel.S(:,~exchange));
Tt = transpose(tModel.T(:,~exchange));

% Parameters
Vmax      = Vmax ./ 4.3478;                                                % slope that minimizes error
R         = 8.3144626 * 1e-3;                                              % (KJ/mol K)
T         = 310.15;                                                        % (K) | 37 °C
F         = 96485;                                                         % Faraday constant C/mol
Hint     = exp(tModel.lbc(find(strcmp(tModel.mets,'c242[c]'))));           % Internal proton concentration
Hext     = exp(tModel.lbc(find(strcmp(tModel.mets,'c242[e]'))));           % External proton concentration
vpH = zeros(size(tModel.mets));
vpH(contains(tModel.mets,'[c]')) = R*T*log(Hint);
vpH(contains(tModel.mets,'[e]')) = R*T*log(Hext);
dpH = abs(-log10(Hint) - -log10(Hext));                                    % pH difference
potential = 1e-3 * (33.33 * dpH - 143.33);                                 % V  = J/C
potential = 1e-3 * potential;                                              % kV = kJ/C

% Auxiliar variables
nMetsMass    = size(model.S,1);
nRxnsMass    = size(model.S,2);
nMetsThermo  = size(St,2);

% Forcing flux condition (non-blocked rxns)
fflb = tModel.lb;
ffub = tModel.ub;
fflb(tModel.lb == 0) =  0.001;
ffub(tModel.ub == 0) = -0.001;

% Equality constraints
LPproblem.Aeq = [model.S zeros(nMetsMass,2*nMetsThermo)];
LPproblem.beq = zeros(nMetsMass,1);

% Inequality constraints
I = eye(nRxnsMass);
I(exchange,:) = [];
LPproblem.Aineq = [fluxDir.*I R*T*Vmax.*St Vmax.*St];
LPproblem.bineq = -Vmax.*(Tt*(tModel.dGf0_H + (tModel.NH .* vpH))...
    + F*potential*tModel.charVec(~exchange));

% dGr matrix
LPproblem.dGrMatrix = [zeros(size(St,1),nRxnsMass) R*T*St St];
LPproblem.dGrVector = Tt*(tModel.dGf0_H + (tModel.NH .* vpH)) + F*potential*tModel.charVec(~exchange); % Constant vector

% Removing H2Ot from thermodynamical analysis
rxnNames = tModel.rxnNames(~exchange);
LPproblem.Aineq(contains(rxnNames,'H2Ot5'),:) = zeros(1,size(LPproblem.Aineq,2));
LPproblem.bineq(contains(rxnNames,'H2Ot5'))   = 0;

% Adding reduction constraints
[redCons, redB] = closa_addReductionConstraints(tModel);
LPproblem.Aineq = [LPproblem.Aineq; zeros(size(redCons,1),nRxnsMass) redCons zeros(size(redCons,1),nMetsThermo)];
LPproblem.bineq = [LPproblem.bineq; redB];

% Bounds & Options
LPproblem.lb = [fflb; tModel.lbc; tModel.dGftLB];
LPproblem.ub = [ffub; tModel.ubc; tModel.dGftUB];
LPproblem.options = optimoptions('linprog','Display','none');
LPproblem.solver  = 'linprog';

end
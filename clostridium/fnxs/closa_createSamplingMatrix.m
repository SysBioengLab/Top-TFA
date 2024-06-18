function matrix_struct = closa_createSamplingMatrix(model,bounds)
% Construct the matrix to be sampled. Define de required fields for execute
% sampling.

%% Define models, bounds and fixed variables
tModel = model;  % With closa they are the same. But recieve different treatment.
tModel.lb  = bounds(1:size(tModel.S,2),1); % This defines the rxn direction
tModel.ub  = bounds(1:size(tModel.S,2),2); % This defines the rxn direction
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

% Equality constraints
Aeq = [model.S zeros(nMetsMass,2*nMetsThermo)];
beq = zeros(nMetsMass,1);

% Inequality constraints
I = eye(nRxnsMass);
I(exchange,:) = [];
Aineq = [fluxDir.*I R*T*Vmax.*St Vmax.*St];
bineq = -Vmax.*(Tt*(tModel.dGf0_H + (tModel.NH .* vpH)) + F*potential*tModel.charVec(~exchange));

% Save dGr matrix 4 later
matrix_struct.dGrMatrix = [R*T*St St Tt*(tModel.dGf0_H + (tModel.NH .* vpH)) + F*potential*tModel.charVec(~exchange)]; % Last column is constant
matrix_struct.dGrRxns = model.rxns(~exchange);

% Removing H2Ot from thermodynamical analysis
rxnNames = tModel.rxnNames(~exchange);
Aineq(contains(rxnNames,'H2Ot5'),:) = zeros(1,size(Aineq,2));
bineq(contains(rxnNames,'H2Ot5'))   = 0;

% Adding reduction constraints
[redCons, redB] = closa_addReductionConstraints(tModel);
Aineq = [Aineq; zeros(size(redCons,1),nRxnsMass) redCons zeros(size(redCons,1),nMetsThermo)];
bineq = [bineq; redB];

%% Add slack variables
nSlack = size(Aineq,1);
Aeq   = [Aeq zeros(size(Aeq,1),nSlack)];
Aineq = [Aineq eye(nSlack)];
matrix_struct.S  = [Aeq; Aineq];
matrix_struct.lb = [bounds(:,1); zeros(nSlack,1)];
matrix_struct.ub = [bounds(:,2); ones(nSlack,1)*1e6];
matrix_struct.b  = [beq; bineq];

sl1 = cell(nSlack,1);
sl1(:) = {'slack_'};
sl2 = cellstr(string(1:nSlack)');
sl3 = strcat(sl1,sl2);

% dGft
dG1 = cell(nMetsThermo,1);
dG1(:) = {'dGft_'};
dG2 = strcat(dG1,tModel.mets);

matrix_struct.vars = [tModel.rxns; tModel.mets; dG2; sl3];

end
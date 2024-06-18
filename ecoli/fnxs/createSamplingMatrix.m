function matrix_struct = createSamplingMatrix(bigModel,bounds,TMFA)
% Construct the matrix to be sampled. Define de required fields for execute
% sampling.

%% Define models, bounds and fixed variables
model  = bigModel.mass;
tModel = bigModel.thermo;
tModel.lb  = bounds(1:size(tModel.S,2),1);
tModel.ub  = bounds(1:size(tModel.S,2),2);
exchange = findExRxns(tModel);
tModel = removeProtons(tModel); % Removing protons
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
vpH = zeros(size(tModel.mets));
vpH(contains(tModel.mets,'_c')) = R*T*log(10^(-7.0));
vpH(contains(tModel.mets,'_e')) = R*T*log(10^(-6.3));
dpH = 0.7;
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
if TMFA
    Z = zeros(sum(~exchange),nRxnsMass);
    Aineq = [Z R*T*fluxDir.*St fluxDir.*St];
    bineq = -fluxDir.*(Tt*(tModel.dGf0_H + (tModel.NH .* vpH)) + F*potential*tModel.charVec(~exchange)) - 0.00001;
else
    I = eye(nRxnsMass);
    I(exchange,:) = [];
    Aineq = [fluxDir.*I R*T*Vmax.*St Vmax.*St];
    bineq = -Vmax.*(Tt*(tModel.dGf0_H + (tModel.NH .* vpH)) + F*potential*tModel.charVec(~exchange));
end

% Save dGr matrix 4 later
matrix_struct.dGrMatrix = [R*T*St St Tt*(tModel.dGf0_H + (tModel.NH .* vpH)) + F*potential*tModel.charVec(~exchange)]; % Lats column is constant
matrix_struct.dGrRxns = model.rxns(~exchange);

% Removing H2Ot from thermodynamical analysis
rxnNames = tModel.rxns(~exchange);
if sum(contains(rxnNames,'H2Ot')) > 0
    Aineq(contains(rxnNames,'H2Ot'),:) = zeros(1,size(Aineq,2));
    bineq(contains(rxnNames,'H2Ot'))   = 0;
end

% Adding reduction constraints
[redCons, redB] = addReductionConstraints(tModel);
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
function closa_addThermoFields(model,meditions)
% Adds dGf and concentration limits for all metabolites. Also adds flux
% data.

%% Flux data
% Biological replicate 1 (High BC) (SIM 4) is used.
model.lb(find(strcmp(model.rxns,'c1311_e_b'))) = -13.099;   % H2 uptake
model.ub(find(strcmp(model.rxns,'c1311_e_b'))) = -12.873;   % H2 uptake
model.lb(find(strcmp(model.rxns,'EX_c1405[e]'))) = -31.56;  % CO uptake
model.ub(find(strcmp(model.rxns,'EX_c1405[e]'))) = -29.57;  % CO uptake
model.lb(find(strcmp(model.rxns,'c1298_e_b'))) =  0;        % Ethanol production
model.ub(find(strcmp(model.rxns,'c1298_e_b'))) =  3.82;     % Ethanol production
model.lb(find(strcmp(model.rxns,'c241_e_b')))  =  3.87;     % Acetate production
model.ub(find(strcmp(model.rxns,'c241_e_b')))  =  10.14115; % Acetate production
model.lb(find(strcmp(model.rxns,'c1305_e_b'))) =  0;        % 2,3-BDO production
model.ub(find(strcmp(model.rxns,'c1305_e_b'))) =  0.15686;  % 2,3-BDO production
model.lb(find(strcmp(model.rxns,'c1466_e_b'))) =  8.35851;  % CO2 production
model.ub(find(strcmp(model.rxns,'c1466_e_b'))) =  12.91;    % CO2 production

%% Concentration data
% Concentration units is molar (M)
% Loose bounds
model.lbc         = log(1e-7) .* ones(size(model.S,1),1);          % 10^-7 M
model.ubc         = log(0.01) .* ones(size(model.S,1),1);          % 0.01 M
% External pH
model.lbc(find(strcmp(model.mets,'c242[e]'))) = log(10^(-5));      % pH 5
model.ubc(find(strcmp(model.mets,'c242[e]'))) = log(10^(-5));      % pH 5
% Internal pH
model.lbc(find(strcmp(model.mets,'c242[c]'))) = log(1e-6);         % pH 6
model.ubc(find(strcmp(model.mets,'c242[c]'))) = log(1e-6);         % pH 6
% Water
model.lbc(find(strcmp(model.mets,'c1392[c]'))) = log(1);           % 1 M
model.ubc(find(strcmp(model.mets,'c1392[c]'))) = log(1);           % 1 M
model.lbc(find(strcmp(model.mets,'c1392[e]'))) = log(1);           % 1 M
model.ubc(find(strcmp(model.mets,'c1392[e]'))) = log(1);           % 1 M

% Acetolactate (fixed)
model.lbc(find(strcmp(model.mets,'c1432[c]'))) = log(2.15977e-06);
model.ubc(find(strcmp(model.mets,'c1432[c]'))) = log(2.15977e-06);

% Meditions
% Extracted from (Valgepea, 2017) DOI:https://doi.org/10.1016/j.cels.2017.04.008
if meditions
    % NADH
    model.lbc(find(strcmp(model.mets,'c638[c]'))) = log(7.09823e-07);
    model.ubc(find(strcmp(model.mets,'c638[c]'))) = log(1.66075e-05);
    % NAD
    model.lbc(find(strcmp(model.mets,'c1164[c]'))) = log(0.001037366);
    model.ubc(find(strcmp(model.mets,'c1164[c]'))) = log(0.00149531);
    % NADPH
    model.lbc(find(strcmp(model.mets,'c1130[c]'))) = log(3.42402e-05);
    model.ubc(find(strcmp(model.mets,'c1130[c]'))) = log(4.56127e-05);
    % NADP
    model.lbc(find(strcmp(model.mets,'c624[c]'))) = log(0.000507405);
    model.ubc(find(strcmp(model.mets,'c624[c]'))) = log(0.000536135);
    % Pyruvate
    model.lbc(find(strcmp(model.mets,'c1337[c]'))) = log(3.55259e-05);
    model.ubc(find(strcmp(model.mets,'c1337[c]'))) = log(7.83589e-05);
    % Acetyl-CoA
    model.lbc(find(strcmp(model.mets,'c947[c]'))) = log(1.33153e-05);
    model.ubc(find(strcmp(model.mets,'c947[c]'))) = log(0.000146784);
    % Acetylphosphate
    model.lbc(find(strcmp(model.mets,'c17[c]'))) = log(8.01318e-05);
    model.ubc(find(strcmp(model.mets,'c17[c]'))) = log(0.000247642);
    
    % Mean and std of mets (In order)
    model.concMu     = [8.65866e-06 0.001266338 3.99265e-05 0.00052177 5.69424e-05 8.00495e-05 0.000163887]';
    model.concStd    = [3.97442e-06 0.000114486 2.84314e-06 7.18246e-6 1.07083e-05 3.33671e-05 4.18776e-05]';
    model.concNormalMets = {'c638[c]','c1164[c]','c1130[c]','c624[c]',...
        'c1337[c]','c947[c]','c17[c]'}; % Leave empty if no meditions available
else
    model.concNormalMets = {}; % Leave empty if no meditions available
end

%% dG incorporation
% Using data from eQuilibrator API
load("ClosA_Thermo_Info.mat","dGft_mu","dGft_cov")

dGft_mu = dGft_mu';
banMets = isnan(dGft_mu);
dGft_mu(banMets) = 0;

dGft_cov(banMets,:) = zeros(sum(banMets),size(dGft_cov,2));
tempCov = dGft_cov;
banIndexes = find(banMets);
banIndexes = banIndexes';
for i=banIndexes
    tempCov = [tempCov(:,1:i-1) zeros(size(tempCov,1),1) tempCov(:,i:end)];
    tempCov(i,i) = 1e6;
end

tempCov = nearestSPD(tempCov);   % Treatment of small differences. John D'Errico (2023)

% Saving variables
model.dGft_mu  = dGft_mu;
model.dGft_cov = tempCov;
model.dGft_std = sqrt(diag(model.dGft_cov));

% Add H+ data
R         = 8.3144626 * 1e-3;
T         = 310.15;
model.dGf0_H  = 0;
model.dGft_Hc = model.dGf0_H - 1*(model.dGf0_H + R*T*model.lbc(find(strcmp(model.mets,'c242[c]'))));
model.dGft_He = model.dGf0_H - 1*(model.dGf0_H + R*T*model.lbc(find(strcmp(model.mets,'c242[e]'))));

% dGft bounds
model.dGftLB = [model.dGft_He; model.dGft_Hc; model.dGft_mu - 2*model.dGft_std];
model.dGftUB = [model.dGft_He; model.dGft_Hc; model.dGft_mu + 2*model.dGft_std];

% Add transport matrix and charge vector
model = closa_completeTherModel(model);

% Add nH of each met
model = generateNH(model);

% Define reversible rxns
model.rev = (model.lb<0 & model.ub>0);

% Name change
thermoRCA = model;

save(strcat('models/thermo_rca.mat'),"thermoRCA")
clear;

end
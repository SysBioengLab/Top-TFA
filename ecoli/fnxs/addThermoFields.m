function addThermoFields(vX)
% Function that generates thermoVX. Adds dGf and concentration
% limits for all metabolites.

load("models/e_coli_core.mat","e_coli_core");
model = bounds1000to100(e_coli_core);
% model = e_coli_core;
% Aditional: ATPS4r fixation
model.lb(contains(model.rxns,'ATPS4r')) = 0;

% Concentration units is mol/L (M)
% Data extracted from TMFA (2007)
model.lbc         = log(1e-7) .* ones(size(model.S,1),1);        % 10^-7 M
model.ubc         = log(0.01) .* ones(size(model.S,1),1);        % 0.01 M
% Water
model.lbc(find(strcmp(model.mets,'h2o_c'))) = log(1);            % 1 M
model.ubc(find(strcmp(model.mets,'h2o_c'))) = log(1);            % 1 M
model.lbc(find(strcmp(model.mets,'h2o_e'))) = log(1);            % 1 M
model.ubc(find(strcmp(model.mets,'h2o_e'))) = log(1);            % 1 M

% Internal pH
model.lbc(find(strcmp(model.mets,'h_c'))) = log(10^(-7.0));      % pH 7
model.ubc(find(strcmp(model.mets,'h_c'))) = log(10^(-7.0));      % pH 7
% External pH
model.lbc(find(strcmp(model.mets,'h_e'))) = log(10^(-6.3));      % pH 6.3
model.ubc(find(strcmp(model.mets,'h_e'))) = log(10^(-6.3));      % pH 6.3

% Concentrations of nutrients in media (TMFA)
model.lbc(find(strcmp(model.mets,'pi_e'))) = log(.056);          % Phospate
model.ubc(find(strcmp(model.mets,'pi_e'))) = log(.056);          % Phospate
model.lbc(find(strcmp(model.mets,'nh4_e'))) = log(.019);         % Ammonium
model.ubc(find(strcmp(model.mets,'nh4_e'))) = log(.019);         % Ammonium
model.lbc(find(strcmp(model.mets,'co2_e'))) = log(.0001);        % CO2
model.ubc(find(strcmp(model.mets,'co2_e'))) = log(.0001);        % CO2
model.lbc(find(strcmp(model.mets,'o2_e'))) = log(8.2e-6);        % O2
model.ubc(find(strcmp(model.mets,'o2_e'))) = log(8.2e-6);        % O2
model.lbc(find(strcmp(model.mets,'glc__D_e'))) = log(.0222);     % Glucose
model.ubc(find(strcmp(model.mets,'glc__D_e'))) = log(.0222);     % Glucose

% Add H+ data
R         = 8.3144626 * 1e-3;
T         = 310.15;
model.dGf0_H  = 0;
model.dGft_Hc = model.dGf0_H - 1*(model.dGf0_H + R*T*log(10^(-7.0)));
model.dGft_He = model.dGf0_H - 1*(model.dGf0_H + R*T*log(10^(-6.3)));

%% dG incorporation
% Using data from eQuilibrator API
load("Ecore_Thermo_Info_pHc7.mat","dGft_mu","dGft_cov")

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
model = generateNH(model);

% dGft bounds
model.dGftLB = model.dGft_mu - 2*model.dGft_std;
model.dGftUB = model.dGft_mu + 2*model.dGft_std;
hcindx = find(strcmp(model.mets,'h_c'));
model.dGftLB = [model.dGftLB(1:hcindx-1); model.dGft_Hc; model.dGft_He; model.dGftLB(hcindx:end)];
model.dGftUB = [model.dGftUB(1:hcindx-1); model.dGft_Hc; model.dGft_He; model.dGftUB(hcindx:end)];

model.description = strcat('COREv',num2str(vX));

if vX >= 2
    % Add concentration data
    model = addConcentrationData(model);
else
    model.concNormalMets = {}; % No meditions
end

% Add transport matrix and charge vector
model = completeTherModel(model);

save(strcat('models/COREv',num2str(vX),'.mat'),"model")
clear;

end
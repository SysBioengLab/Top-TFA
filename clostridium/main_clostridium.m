function main_clostridium(medition,options)
% main_clostridium_autoethanogenum

% In this script we are going to apply Top-TMFA to the model of Clostridum
% autoethanogenum.
% clear; clc
% addpath models\
addpath files\
addpath fnxs\
addpath clusterFiles\
addpath '..\DFS\'
addpath '..\DFS\fxns'

load("models\rcaModel.mat","rcaModel") % Load model

% Generate thermodynamic model
closa_addThermoFields(rcaModel,medition)
clear rcaModel

% FVA run
load("models\thermo_rca.mat","thermoRCA")
thermoRCA = updateBounds(thermoRCA);

% Topology search
clusterID = 'thermoRCA_cluster';
findTopologiesNewClusters3(thermoRCA,clusterID)
load(strcat('tfm_cluster',clusterID,'.mat'),"tfm","model")
% Remove rxns that form a loop
thermoRCA = closa_rmvLoopyRxns(thermoRCA,model);
tfm = reorderTFM(tfm,model,thermoRCA);
% Add Vmax and Vmin fields
thermoRCA.vmax = thermoRCA.ub;
thermoRCA.vmin = thermoRCA.lb;

% Discard non mass-feasible topologies
[tfm, ~] = discardNonMassFeasibleTopologies(thermoRCA,tfm);

% Run Top-TMFA
[feasible,boundsTcTMFA,boundsFVA,dGrbounds] = closa_checkFeasibility(thermoRCA,tfm);
disp(feasible)

% Save data
if medition
    save(strcat("results\DATA_v2.mat"),"thermoRCA","tfm","boundsFVA", ...
        "feasible","boundsTcTMFA","dGrbounds")
else
    save(strcat("results\DATA_v1.mat"),"thermoRCA","tfm","boundsFVA", ...
        "feasible","boundsTcTMFA","dGrbounds")
end

% Run if two versions are updated:
% closa_data2Excel()

%% Sampling
addpath '..\looplessFxns\'
% if ~medition
%     name1 = "test1_ADSB_NU";
% else
%     name1 = options.name;
% end
name1 = options.name;
nTop = sum(feasible);
volProp = ones(nTop,1) / nTop;
bounds = boundsTcTMFA(:,:,feasible);
options.numSamples = 1e5;
options.diagnostics = 0;
options.stepsPerPoint = 10; % Thinning factor
options.algorithm = 'ADSB'; % CHRR ADSB
options.numDiscarded = 1e3;
options.uniform = 0;
options.discreteFactor = 10000;
options.parallelFlag = 1;
options.loopless = 0;
options.populationScale = 15;
closa_tcTMFAsampler(thermoRCA,bounds,volProp,name1,options)

end
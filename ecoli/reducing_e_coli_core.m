% Reducing e_coli_core
% The first steps are the same to avoid repeated work

clear;
addpath 'D:\UC\Magíster\MATLAB\current\models'            % load model
addpath 'D:\UC\Magíster\MATLAB\current\fnxs'              % load required functions
% initCobraToolbox()

%% Fixation of concentration limits
% Se revisará el caso de crecimiento aeróbico con glucosa como fuente de
% carbono
load("e_coli_core.mat");
model = e_coli_core;
[minFlux, maxFlux] = fluxVariability(model);
% La unidad de las concentraciones será molar (M)
% Se ocupan datos del paper de TMFA (2007)
model.lbc         = log(1e-5) .* ones(size(model.S,1),1);        % 10^-5 M
model.ubc         = log(0.02) .* ones(size(model.S,1),1);        % 0.02 M
% Water
model.lbc(find(strcmp(model.mets,'h2o_c'))) = log(32);           % 32 M
model.ubc(find(strcmp(model.mets,'h2o_c'))) = log(32);           % 32 M
model.lbc(find(strcmp(model.mets,'h2o_e'))) = log(32);           % 32 M
model.ubc(find(strcmp(model.mets,'h2o_e'))) = log(32);           % 32 M
% Internal pH
model.lbc(find(strcmp(model.mets,'h_c'))) = log(1e-7);           % pH 7
model.ubc(find(strcmp(model.mets,'h_c'))) = log(1e-7);           % pH 7
% External pH
model.lbc(find(strcmp(model.mets,'h_e'))) = log(1e-11);          % pH 11
model.ubc(find(strcmp(model.mets,'h_e'))) = log(1e-4);           % pH 4
% Internal o2
model.lbc(find(strcmp(model.mets,'o2_c'))) = log(1e-7);          % 10^-7 M
model.ubc(find(strcmp(model.mets,'o2_c'))) = log(8.2e-6);        % 8.2e-6 M
% Concentrations of nutrients in media
model.lbc(find(strcmp(model.mets,'pi_e'))) = log(.056);          % Phospate
model.ubc(find(strcmp(model.mets,'pi_e'))) = log(.056);          % Phospate
model.lbc(find(strcmp(model.mets,'nh4_e'))) = log(.019);         % Ammonium
model.ubc(find(strcmp(model.mets,'nh4_e'))) = log(.019);         % Ammonium
model.lbc(find(strcmp(model.mets,'co2_e'))) = log(.0001);        % CO2
model.ubc(find(strcmp(model.mets,'co2_e'))) = log(.0001);        % CO2
model.lbc(find(strcmp(model.mets,'o2_e'))) = log(8.2e-6);        % O2
model.ubc(find(strcmp(model.mets,'o2_e'))) = log(8.2e-6);        % O2
model.lbc(find(strcmp(model.mets,'glc__D_e'))) = log(.02);       % Glucose
model.ubc(find(strcmp(model.mets,'glc__D_e'))) = log(.02);       % Glucose
% Lista de compuestos no encontrados en modelo: sulfato, sodio, potasio y
% Fe+2
% ¿Se deberían fijar el resto de los mets externos en cero?

%% dG incorporation
% Cargar datos de eQuilibrator API
load("Ecore_Thermo_Info.mat")
comp = double(comp);
comp = comp';
dGf = dGf';
% Corregir dGf de compuestos problematicos
dGf(13)   = -1945.427899; % adp
dGf(16)   = -1033.960494; % amp
dGf(17)   = -2811.578332; % atp
dGf(30)   = -915.5829776; % d-fructose
dGf(45)   = -1155.773217; % iso-citrate
model.dGf = dGf;
% Agregar cargas
% Copy-paste de hoja de excel
model.charge = [-4 -3 -3 -3 -2 -1 -1 0 0 -4 -3 -2 -2 -2 -2 -2 -3 -3 0 0 -4 -2 -2 0 0 -2 -4 -1 -1 0 -2 -2 -2 -2 0 0 0 -1 -1 -1 0 0 1 1 -3 -1 -1 -2 -2 -1 -1 -3 -3 1 1 0 0 -2 -3 -2 -2 -1 -1 0 0 -2 -2 -2 -2 -2 -5 -2]';
%% Define transport reactions
model.comp = comp;
ext = logical(sum(abs(model.S) .* comp,1));
int = logical(sum(abs(model.S) .* ~comp,1));
model.transport_rxns = logical(int & ext);                                 % Rxns that involve intra and extra celular metabolites
model.exchange_rxns  = findExcRxns(model);

%% Reduce model (new step)
% En primera instancia queremos quedarnos con una vía pequeña, que ojalá
% tenga reacciones de transporte para poner a prueba las restricciones
% nuevas.

% Incorporar manualmente las rxns de interés
rmvdRxns = false(size(model.rxns));
rmvdRxns = rmvdRxns + contains(model.rxns,'EX_glc__D_e');
rmvdRxns = rmvdRxns + contains(model.rxns,'GLCpts');
rmvdRxns = rmvdRxns + contains(model.rxns,'PGI');
rmvdRxns = rmvdRxns + contains(model.rxns,'PFK');
rmvdRxns = rmvdRxns + contains(model.rxns,'FBP');
rmvdRxns = rmvdRxns + contains(model.rxns,'FBA');
rmvdRxns = rmvdRxns + contains(model.rxns,'TPI');
rmvdRxns = rmvdRxns + contains(model.rxns,'GAPD');
rmvdRxns = rmvdRxns + contains(model.rxns,'PGK');
rmvdRxns = rmvdRxns + contains(model.rxns,'PGM');
rmvdRxns = rmvdRxns + contains(model.rxns,'ENO');
rmvdRxns = rmvdRxns + contains(model.rxns,'PYK');
rmvdRxns = rmvdRxns + contains(model.rxns,'PPS');
rmvdRxns = rmvdRxns + contains(model.rxns,'PYRt2');
rmvdRxns = rmvdRxns + contains(model.rxns,'EX_pyr_e');
rmvdRxns = rmvdRxns + contains(model.rxns,'THD2');
rmvdRxns = rmvdRxns + contains(model.rxns,'NADTRHD');
rmvdRxns = rmvdRxns + contains(model.rxns,'EX_h_e');
rmvdRxns = rmvdRxns + contains(model.rxns,'ADK1');
rmvdRxns = rmvdRxns + contains(model.rxns,'ATPM');
rmvdRxns = rmvdRxns + contains(model.rxns,'ATPS4r');
rmvdRxns = rmvdRxns + contains(model.rxns,'EX_h2o_e');
rmvdRxns = rmvdRxns + contains(model.rxns,'H2Ot');
rmvdRxns = rmvdRxns + contains(model.rxns,'EX_pi_e');
rmvdRxns = rmvdRxns + contains(model.rxns,'PIt2r');
rmvdRxns = ~rmvdRxns;
% Hacer copia del modelo
model2 = model;
% Remover otros subsistemas
model2.rxns(rmvdRxns)           = [];
model2.rxnNames(rmvdRxns)       = [];
model2.subSystems(rmvdRxns)     = [];
model2.S(:,rmvdRxns)            = [];
model2.lb(rmvdRxns)             = [];
model2.ub(rmvdRxns)             = [];
model2.c(rmvdRxns)              = [];
model2.rev(rmvdRxns)            = [];
model2.transport_rxns(rmvdRxns) = [];
model2.exchange_rxns(rmvdRxns)  = [];
model2.rxnGeneMat(rmvdRxns,:)   = [];
model2.grRules(rmvdRxns)        = [];
% Agregar exchanges (virtuales) de metabolitos necesarios
newExmets = {'pep_c','nad_c','nadh_c'}';
model2 = addExchangeRxn(model2,newExmets);
% Optimización
[minFlux1, maxFlux1] = fluxVariability(model);
[minFlux2, maxFlux2] = fluxVariability(model2);

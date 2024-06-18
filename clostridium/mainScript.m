% mainScript

addpath '..\looplessFxns\'
load("results\DATA_v2.mat","boundsTcTMFA","thermoRCA")

% nReps = 3;

% Common for all cases
nTop = 1;
volProp = 1;
bounds = boundsTcTMFA;
options.numSamples = 5e4;
options.diagnostics = 0;
options.stepsPerPoint = 3000; % Thinning factor
options.numDiscarded = 1e3; % 1e4
options.parallelFlag = 1;
options.uniform = 0;
options.algorithm = 'ADSB';
options.loopless = 0;
options.discreteFactor = 100;

options.name = "non_uni_4";
closa_tcTMFAsampler(thermoRCA,bounds,volProp,options.name,options)

options.name = "non_uni_5";
closa_tcTMFAsampler(thermoRCA,bounds,volProp,options.name,options)

% 2 repeticiones del uniforme con estos parámetros
options.uniform = 1;

options.name = "uni_4";
closa_tcTMFAsampler(thermoRCA,bounds,volProp,options.name,options)

options.name = "uni_5";
closa_tcTMFAsampler(thermoRCA,bounds,volProp,options.name,options)

% for j=1:nReps
%     options.name = strcat("test_noKWB_",num2str(j));
%     closa_tcTMFAsampler(thermoRCA,bounds,volProp,options.name,options)
% end




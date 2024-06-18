%%%%%%%%%%%%%%%%%%%%%%%%% Part 1: Model curation %%%%%%%%%%%%%%%%%%%%%%%%%%

% Adding required functions and models to the path
clear; clc;
% initCobraToolbox(false)
addpath models\
addpath fnxs\
addpath ..\DFS
addpath ..\DFS\fxns\

%% Add thermodynamic data
% Additional data as dGf (obtained from eQuilibrator API) and concentration
% bounds is added to build new models (COREv1, COREv2).

addThermoFields(1) % vX=1
addThermoFields(2) % vX=2

%% Model reduction
% As first try, the model is reduced in order to test all functions.
% Two models are created (ecModel1/2/3 and ecTherModel1/2/3).
ecM = 1:3;
version = 1:2;
for N=ecM
    for vX=version
        modelReduction(N,vX)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%% Part 2: Pre-processing %%%%%%%%%%%%%%%%%%%%%%%%%%
% Run of optimization problem with TMFA and Top-TMFA, for all variables,
% all models and both datasets.

ecM = 1:3;
version = 1:2;
for N=ecM
    disp(["ecModel: ",num2str(N)])
    my_dir = strcat('files\ec',num2str(N),"\");
    for vX=version
        tStart = tic;
        disp(["Version: ",num2str(vX)])
        load(strcat("models\subModels\ec",num2str(N),"v",num2str(vX),".mat"))
        if N == 1
            bigModel.mass   = ecModel1;
            bigModel.thermo = ecTherModel1;
        elseif N == 2
            bigModel.mass   = ecModel2;
            bigModel.thermo = ecTherModel2;
        else
            bigModel.mass   = ecModel3;
            bigModel.thermo = ecTherModel3;
        end
        % Original
        ogBounds = [bigModel.thermo.lb bigModel.thermo.ub; bigModel.thermo.lbc bigModel.thermo.ubc];
        % Topology search
        clusterID = strcat('_ec',num2str(N),'v',num2str(vX));
        findTopologiesNewClusters3(bigModel.mass,clusterID)
        load(strcat('clusterFiles\tfm_cluster',clusterID,'.mat'),"tfm","model")
        % Remove rxns that form a loop
        [bigModel.mass,bigModel.thermo] = rmvLoopyRxns(bigModel.mass,bigModel.thermo,model);
        tfm = reorderTFM(tfm,model,bigModel.mass);
        % Add Vmax and Vmin fields
        [bigModel.thermo.lb,bigModel.thermo.ub] = FVA(bigModel.mass);
        bigModel.mass.lb = bigModel.thermo.lb;
        bigModel.mass.ub = bigModel.thermo.ub;
        bigModel.thermo.vmax = bigModel.thermo.ub;
        bigModel.thermo.vmin = bigModel.thermo.lb;
        % Run TMFA
        relax = true;
        [feasibleTMFA,boundsTMFA,boundsFVA,fluxSolutionsTMFA,dGrboundsTMFA,timeTMFA] = checkFeasibility(bigModel.mass,bigModel.thermo,tfm,relax);
        % Run Top-TMFA
        relax = false;
        [feasibleFcTMFA,boundsFcTMFA,~,fluxSolutionsFcTMFA,dGrboundsFcTMFA,timeFcTMFA] = checkFeasibility(bigModel.mass,bigModel.thermo,tfm,relax);
        preProcessingTime = toc(tStart);
        % Save data
        save(strcat(my_dir,"DATAv",num2str(vX),".mat"),"bigModel","tfm","feasibleTMFA", ...
            "feasibleFcTMFA","ogBounds","boundsFVA","boundsTMFA","boundsFcTMFA",...
            "dGrboundsTMFA","dGrboundsFcTMFA","fluxSolutionsFcTMFA", ...
            "fluxSolutionsTMFA","preProcessingTime","timeTMFA","timeFcTMFA")
    end
end

data2Excel(ecM,version) % Generate tables on excel

%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 3: Sampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Test sampling (Top-TMFA)
load("files\ec3\DATAv2.mat")
nTop = 1;
volProp = 1;
bounds = boundsFcTMFA(:,:,1);
name = 'topTFA_ADS_pro';
options.numSamples = 1e5;
options.diagnostics = 0;
options.algorithm = 'ADSB'; % CHRR ADSB
options.parallelFlag = 1;
options.stepsPerPoint = 1e3;
topTFAsampler(bigModel,bounds,volProp,name,options)
% Test with TMFA
name = 'TMFA_ADS_pro';
options.TMFA = 1;
bounds = boundsTMFA(:,:,1);
topTFAsampler(bigModel,bounds,volProp,name,options)

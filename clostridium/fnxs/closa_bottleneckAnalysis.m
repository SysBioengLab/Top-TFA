function bottleOutputs = closa_bottleneckAnalysis(path)

% Load data
% path = "results\DATA_v2.mat";
load(path,"thermoRCA","boundsTcTMFA","feasible","boundsFVA","dGrbounds")

% Define variables
tModel = thermoRCA;
exchange = findExRxns(tModel);
dGrRxns = tModel.rxns(~exchange);
ndGr = size(dGrRxns,1);
nTop = sum(feasible);
nRxns = size(tModel.rxns,1);
cleanBoundsFVA = boundsFVA;
cleanBoundsFVA(exchange,:,:) = [];
cleanBoundsTcTMFA = boundsTcTMFA(1:nRxns,:,:);
cleanBoundsTcTMFA(exchange,:,:) = [];
clasification = zeros(ndGr,nTop); % 0:NA 1:TB 2:TR 3:TA
maxBoundsFVA = zeros(ndGr,nTop);
maxBoundsTcTMFA = zeros(ndGr,nTop);
maxdGr = zeros(ndGr,nTop);
dGrFVAMatrix = zeros(ndGr,nTop);
dGrTcTMFAMatrix = zeros(ndGr,nTop);

% Main loop (por topología)
for i=find(feasible')
    % Vmax calculation
    fluxDir = ones(size(tModel.rxns));
    fluxDir(exchange) = [];
    fluxDir(cleanBoundsTcTMFA(:,2,i)<=1e-8) = -1;
    Vmax = tModel.vmax; Vmax(exchange) = [];
    Vmin = tModel.vmin; Vmin(exchange) = [];
    Vmax(fluxDir==-1) = Vmin(fluxDir==-1);
    Vmax = Vmax ./ 4.3478;
    % Check bounds
    maxBoundsFVA(fluxDir==1,i) = cleanBoundsFVA(fluxDir==1,2,i);
    maxBoundsFVA(fluxDir==-1,i) = cleanBoundsFVA(fluxDir==-1,1,i);
    maxBoundsTcTMFA(fluxDir==1,i) = cleanBoundsTcTMFA(fluxDir==1,2,i);
    maxBoundsTcTMFA(fluxDir==-1,i) = cleanBoundsTcTMFA(fluxDir==-1,1,i);
    NA = abs(maxBoundsFVA(:,i) - maxBoundsTcTMFA(:,i)) < 1e-6;
    % dGrFVA and dGrFcTMFA calculus
    dGrFVA = -fluxDir.*(maxBoundsFVA(:,i) ./ Vmax);
    dGrTcTMFA = -fluxDir.*(maxBoundsTcTMFA(:,i) ./ Vmax);
    dGrFVAMatrix(:,i) = dGrFVA;
    dGrTcTMFAMatrix(:,i) = dGrTcTMFA;
    % Optimization
    tempModel = tModel;
    tempModel.lb = boundsFVA(:,1,i);
    tempModel.ub = boundsFVA(:,2,i);
    LPproblem = closa_newbuildTherModels(tempModel);
    dGrMatrix = LPproblem.dGrMatrix;
    dGrVector = LPproblem.dGrVector;
    LPproblem = rmfield(LPproblem,{'dGrMatrix' 'dGrVector'});
    count = 0;
    for j=find((~exchange)')
        count = count + 1;
        backUpLB = LPproblem.lb(j); backUpUB = LPproblem.ub(j);
        if boundsTcTMFA(j,2,i)<=1e-8 % Negative reaction, dGr >0
            fixValue = boundsTcTMFA(j,1,i);
            optsense = -1;
            LPproblem.lb(j) = fixValue;
            LPproblem.ub(j) = fixValue;
        else
            fixValue = boundsTcTMFA(j,2,i);
            optsense = 1;
            LPproblem.lb(j) = fixValue;
            LPproblem.ub(j) = fixValue;
        end
        LPproblem.f = optsense*dGrMatrix(count,:);
        [~, fval, flag] = linprog(LPproblem);
        if flag > 0
            maxdGr(count,i) = optsense*fval + dGrVector(count);
        else
            disp("Problem with dGr optimization")
            break
        end
        LPproblem.lb(j) = backUpLB; LPproblem.ub(j) = backUpUB;
    end
    % Clasification
    clasification(abs(dGrTcTMFA-maxdGr(:,i))<=1e-6 & ~NA,i) = 1; % Thermodynamic bottlenecks
    clasification(abs(maxdGr(:,i))-abs(dGrTcTMFA)>1e-6 & abs(dGrFVA) > abs(maxdGr(:,i)) & ~NA,i) = 2; % Thermodynamic restricted
    clasification(abs(dGrFVA) <= abs(maxdGr(:,i)) & ~NA,i) = 3; % Thermodynamic affected
end

% Thermodynamic sensitivity
Tsen = false(ndGr,nTop);
for i=1:nTop
    Tsen(:,i) = min(abs(dGrbounds(:,:,i)),[],2) < 4.3478;
end

bottleOutputs.clasification   = clasification;
bottleOutputs.maxdGr          = maxdGr;
bottleOutputs.dGrFVAMatrix    = dGrFVAMatrix;
bottleOutputs.dGrTcTMFAMatrix = dGrTcTMFAMatrix;
bottleOutputs.dGrRxns         = dGrRxns;
bottleOutputs.Tsen            = Tsen;

end
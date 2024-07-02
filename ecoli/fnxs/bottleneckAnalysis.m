% bottleneckAnalysis

% Load data
load("files\ec3\DATAv2.mat","bigModel","boundsFcTMFA","feasibleFcTMFA", ...
    "dGrboundsFcTMFA","boundsFVA","fluxSolutionsFcTMFA")

% Define variables
model = bigModel.mass;
tModel = bigModel.thermo;
exchange = findExRxns(tModel);
dGrRxns = tModel.rxns(~exchange);
ndGr = size(dGrRxns,1);
nTop = sum(feasibleFcTMFA);
nRxns = size(tModel.rxns,1);
cleanBoundsFVA = boundsFVA(1:nRxns,:,:);
cleanBoundsFVA(exchange,:,:) = [];
cleanBoundsFcTMFA = boundsFcTMFA(1:nRxns,:,:);
cleanBoundsFcTMFA(exchange,:,:) = [];
clasification = zeros(ndGr,nTop); % 0:NA 1:TB 2:TR 3:TA
maxBoundsFVA = zeros(ndGr,nTop);
maxBoundsFcTMFA = zeros(ndGr,nTop);
maxdGr = zeros(ndGr,nTop);
X = zeros(nRxns+2*size(tModel.mets,1),ndGr,nTop);

% Main loop (por topología)
for i=find(feasibleFcTMFA')
    % Vmax calculation
    fluxDir = ones(size(tModel.rxns));
    fluxDir(exchange) = [];
    fluxDir(cleanBoundsFcTMFA(:,2,i)<=1e-8) = -1;
    Vmax = tModel.vmax; Vmax(exchange) = [];
    Vmin = tModel.vmin; Vmin(exchange) = [];
    Vmax(fluxDir==-1) = Vmin(fluxDir==-1);
    Vmax = Vmax ./ 4.3478;
    % Check bounds
    maxBoundsFVA(fluxDir==1,i) = cleanBoundsFVA(fluxDir==1,2,i);
    maxBoundsFVA(fluxDir==-1,i) = cleanBoundsFVA(fluxDir==-1,1,i);
    maxBoundsFcTMFA(fluxDir==1,i) = cleanBoundsFcTMFA(fluxDir==1,2,i);
    maxBoundsFcTMFA(fluxDir==-1,i) = cleanBoundsFcTMFA(fluxDir==-1,1,i);
    NA = abs(maxBoundsFVA(:,i) - maxBoundsFcTMFA(:,i)) < 1e-6;
    % dGrFVA and dGrFcTMFA calculus
    dGrFVA = -fluxDir.*(maxBoundsFVA(:,i) ./ Vmax);
    dGrFcTMFA = -fluxDir.*(maxBoundsFcTMFA(:,i) ./ Vmax);
    % Optimization
    tempModel = tModel;
    tempModel.lb = boundsFVA(1:nRxns,1,i);
    tempModel.ub = boundsFVA(1:nRxns,2,i);
    LPproblem = newbuildTherModels(model,tempModel,0);
    dGrMatrix = LPproblem.dGrMatrix;
    dGrVector = LPproblem.dGrVector;
    LPproblem = rmfield(LPproblem,{'dGrMatrix' 'dGrVector'});
    count = 0;
    for j=find((~exchange)')
        count = count + 1;
        backUpLB = LPproblem.lb(j); backUpUB = LPproblem.ub(j);
        if boundsFcTMFA(j,2,i)<=1e-8 % Negative reaction, dGr >0
            fixValue = boundsFcTMFA(j,1,i);
            optsense = -1;
            LPproblem.lb(j) = fixValue;
            LPproblem.ub(j) = fixValue;
        else
            fixValue = boundsFcTMFA(j,2,i);
            optsense = 1;
            LPproblem.lb(j) = fixValue;
            LPproblem.ub(j) = fixValue;
        end
        LPproblem.f = optsense*dGrMatrix(count,:);
        [x, fval, flag] = linprog(LPproblem);
        if flag > 0
            maxdGr(count,i) = optsense*fval + dGrVector(count);
            X(:,count,i) = x;
        else
            disp("Problem with dGr optimization")
            break
        end
        LPproblem.lb(j) = backUpLB; LPproblem.ub(j) = backUpUB;
    end
    % Clasification
    clasification(abs(dGrFcTMFA-maxdGr(:,i))<=1e-6 & ~NA,i) = 1; % Thermodynamic bottlenecks
    clasification(abs(maxdGr(:,i))-abs(dGrFcTMFA)>1e-6 & abs(dGrFVA) > abs(maxdGr(:,i)) & ~NA,i) = 2; % Thermodynamic restricted
    clasification(abs(dGrFVA) <= abs(maxdGr(:,i)) & ~NA,i) = 3; % Thermodynamic affected
end

% Thermodynamic sensitivity
Tsenv2 = false(ndGr,nTop);
for i=1:nTop
    Tsenv2(:,i) = min(abs(dGrboundsFcTMFA(:,:,i)),[],2) < 4.3478;
end

load("files\ec3\DATAv1.mat","dGrboundsFcTMFA")
Tsenv1 = false(ndGr,nTop);
for i=1:nTop
    Tsenv1(:,i) = min(abs(dGrboundsFcTMFA(:,:,i)),[],2) < 4.3478;
end

% Save on excel
filename = strcat('Results\AdditionalInformation_E_coli.xlsx');
sheet = 'Bottleneck analysis';
writecell(dGrRxns,filename,'Sheet',sheet,'Range','B4')
writematrix(1:nTop,filename,'Sheet',sheet,'Range','C3')
writematrix(clasification,filename,'Sheet',sheet,'Range','C4')

sheet = 'Thermodynamic Sensitivity';
writecell(dGrRxns,filename,'Sheet',sheet,'Range','B4')
writematrix(1:nTop,filename,'Sheet',sheet,'Range','C3')
writematrix(Tsenv1,filename,'Sheet',sheet,'Range','C4')

writecell(dGrRxns,filename,'Sheet',sheet,'Range','J4')
writematrix(1:nTop,filename,'Sheet',sheet,'Range','K3')
writematrix(Tsenv2,filename,'Sheet',sheet,'Range','K4')


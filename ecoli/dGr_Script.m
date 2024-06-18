% dGr script

% Performs the 3 ways to calculate dGr.
% 1: By reaction. Only accounting for concentration bounds. Two versions:
% with and without variation of dGf.
% 2: Using Top-TMFA. This fixes rxn direction and adds the metabolic net 
% constraints.
% 3: top-TMFA + maxFlux. When maximizing |dGr| the maximum flux is fixed.

clear; clc;
% Load data
load("files\ec3\DATAv2.mat","bigModel","boundsFcTMFA","feasibleFcTMFA", ...
    "dGrboundsFcTMFA","boundsFVA","fluxSolutionsFcTMFA")

% Define variables
model    = bigModel.mass;
tModel   = bigModel.thermo;
exchange = findExRxns(tModel);
tModel   = removeProtons(tModel); % Removing protons
dGrRxns  = tModel.rxns(~exchange);
ndGr     = numel(dGrRxns);
St       = transpose(tModel.S(:,~exchange));
Tt       = transpose(tModel.T(:,~exchange));
nTop     = sum(feasibleFcTMFA);
nRxns    = size(tModel.rxns,1);
% Parameters
R        = 8.3144626 * 1e-3;                                               % (KJ/mol K)
T        = 310.15;                                                         % (K) | 25 °C
F        = 96485;                                                          % Faraday constant C/mol
dpH = 0.7;
potential = 1e-3 * (33.33 * dpH - 143.33);                                 % V  = J/C
potential = 1e-3 * potential;                                              % kV = kJ/C
vpH = zeros(size(tModel.mets));
vpH(contains(tModel.mets,'_c')) = R*T*log(10^(-7.0));
vpH(contains(tModel.mets,'_e')) = R*T*log(10^(-6.3));
% Outputs
case_1_withVar    = zeros(ndGr,2);
case_1_withoutVar = zeros(ndGr,2);
case_2            = zeros(ndGr,2);
case_3            = zeros(ndGr,2);
case_3_complete   = zeros(ndGr,2,nTop);

% dGr matrix
dGrMatrix = [zeros(size(St,1),size(model.S,2)) R*T*St St];
dGrVector = Tt*(tModel.dGf0_H + (tModel.NH .* vpH)) + F*potential*tModel.charVec(~exchange); % Constant vector

%% Case 1 (with variation of dGf)

% Optimization
LPproblem.lb = [bigModel.thermo.lb; bigModel.thermo.lbc; bigModel.thermo.dGftLB];
LPproblem.ub = [bigModel.thermo.ub; bigModel.thermo.ubc; bigModel.thermo.dGftUB];
LPproblem.A  = []; LPproblem.Aeq  = [];
LPproblem.b  = []; LPproblem.beq  = [];
LPproblem.options = optimoptions('linprog','Display','none');
LPproblem.solver  = 'linprog';

for i=1:ndGr
    % Minimization
    LPproblem.f = dGrMatrix(i,:);
    [x, fval, flag] = linprog(LPproblem);
    if flag > 0
        case_1_withVar(i,1) = fval + dGrVector(i);
    else
        disp("Problem with dGr optimization")
        break
    end
    % Maximization
    LPproblem.f = -dGrMatrix(i,:);
    [x, fval, flag] = linprog(LPproblem);
    case_1_withVar(i,2) = -fval + dGrVector(i);
end

%% Case 1 (without variation of dGf)

% Optimization
LPproblem.lb = [bigModel.thermo.lb; bigModel.thermo.lbc; bigModel.thermo.dGft_mu];
LPproblem.ub = [bigModel.thermo.ub; bigModel.thermo.ubc; bigModel.thermo.dGft_mu];
LPproblem.A  = []; LPproblem.Aeq  = [];
LPproblem.b  = []; LPproblem.beq  = [];
LPproblem.options = optimoptions('linprog','Display','none');
LPproblem.solver  = 'linprog';

for i=1:ndGr
    % Minimization
    LPproblem.f = dGrMatrix(i,:);
    [x, fval, flag] = linprog(LPproblem);
    if flag > 0
        case_1_withoutVar(i,1) = fval + dGrVector(i);
    else
        disp("Problem with dGr optimization")
        break
    end
    % Maximization
    LPproblem.f = -dGrMatrix(i,:);
    [x, fval, flag] = linprog(LPproblem);
    case_1_withoutVar(i,2) = -fval + dGrVector(i);
end

%% Case 2 (min/max of all topologies)

myLB = reshape(dGrboundsFcTMFA(:,1,:),[ndGr nTop]);
myUB = reshape(dGrboundsFcTMFA(:,2,:),[ndGr nTop]);
case_2(:,1) = min(myLB,[],2);
case_2(:,2) = max(myUB,[],2);

%% Case 3 complete

for i=find(feasibleFcTMFA')
    tempModel    = bigModel.thermo;
    tempModel.lb = boundsFcTMFA(1:nRxns,1,i);
    tempModel.ub = boundsFcTMFA(1:nRxns,2,i);
    LPproblem = newbuildTherModels(bigModel.mass,tempModel,0); % We obtain matrices from top-TMFA
    % Optimization
    count = 0;
    for j=find((~exchange)')
        count = count + 1;
        backUpLB = LPproblem.lb(j); backUpUB = LPproblem.ub(j);
        if boundsFcTMFA(j,2,i)<=1e-8 % Negative reaction, dGr >0
            fixValue = boundsFcTMFA(j,1,i);
            optsense = -1;
            LPproblem.lb(j) = fixValue;
            LPproblem.ub(j) = fixValue;
            col = 2;
            case_3_complete(count,1,i) = dGrboundsFcTMFA(count,1,i);
        else
            fixValue = boundsFcTMFA(j,2,i);
            optsense = 1;
            LPproblem.lb(j) = fixValue;
            LPproblem.ub(j) = fixValue;
            col = 1;
            case_3_complete(count,2,i) = dGrboundsFcTMFA(count,2,i);
        end
        LPproblem.f = optsense*dGrMatrix(count,:);
        [x, fval, flag] = linprog(LPproblem);
        if flag > 0
            case_3_complete(count,col,i) = optsense*fval + dGrVector(count);
        else
            disp("Problem with dGr optimization")
            disp(i)
            disp(count)
            disp(j)
            break
        end
        LPproblem.lb(j) = backUpLB; LPproblem.ub(j) = backUpUB;
    end
end

%% Case 3 (min/max of all topologies)

myLB = reshape(case_3_complete(:,1,:),[ndGr nTop]);
myUB = reshape(case_3_complete(:,2,:),[ndGr nTop]);
case_3(:,1) = min(myLB,[],2);
case_3(:,2) = max(myUB,[],2);

%% Save on Excel

filename = strcat('Results\Excels\final_ec3.xlsx');
sheet = 'dGr_comparison';
writecell(dGrRxns,filename,'Sheet',sheet,'Range','B4')
% writematrix(1:nTop,filename,'Sheet',sheet,'Range','C3')
writematrix([case_1_withoutVar case_1_withVar case_2 case_3],filename,'Sheet',sheet,'Range','C4')

sheet = 'dGr_by_Topology_case2';
writecell(dGrRxns,filename,'Sheet',sheet,'Range','B4')
topVec = zeros(1,2*nTop);
for i=1:nTop; topVec(2*i-1:2*i) = i; end
writematrix(topVec,filename,'Sheet',sheet,'Range','C3')
writematrix(dGrboundsFcTMFA,filename,'Sheet',sheet,'Range','C4')

sheet = 'dGr_by_Topology_case3';
writecell(dGrRxns,filename,'Sheet',sheet,'Range','B4')
writematrix(topVec,filename,'Sheet',sheet,'Range','C3')
writematrix(case_3_complete,filename,'Sheet',sheet,'Range','C4')

%% Save on Matlab
save("DGR_cases.mat","case_1_withoutVar","case_1_withVar","case_2", ...
    "case_3_complete","case_3","dGrRxns")

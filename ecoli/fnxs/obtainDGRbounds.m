function dGrbounds = obtainDGRbounds(path_to_data)
% Loads data from path_to_data and calculates dGr of models reactions.
% Saves and returns dGrboudns.

%% Load data
load(path_to_data,"bigModel","boundsTMFA","boundsFcTMFA","feasibleTMFA","feasibleFcTMFA");

%% Clean bounds for TMFA and tc-TMFA
feasibleTopologies = 1:size(feasibleFcTMFA,1);
if isequal(feasibleFcTMFA,feasibleTMFA)
    boundsFcTMFA(:,:,~feasibleFcTMFA) = [];
    boundsTMFA(:,:,~feasibleTMFA) = [];
    feasibleTopologies(~feasibleFcTMFA) = [];
else
    newFeasible = logical(feasibleFcTMFA .* feasibleTMFA);
    boundsFcTMFA(:,:,~newFeasible) = [];
    boundsTMFA(:,:,~newFeasible) = [];
    feasibleTopologies(~newFeasible) = [];
end

%% Build dGr Matrix
tModel = bigModel.thermo;
exchange = findExRxns(tModel);
tModel = removeProtons(tModel); % Removing protons
St = transpose(tModel.S(:,~exchange));
Tt = transpose(tModel.T(:,~exchange));

% Parameters
R         = 8.3144626 * 1e-3;                                              % (KJ/mol K)
T         = 298.15;                                                        % (K) | 25 °C
F         = 96485;                                                         % Faraday constant C/mol
vpH = zeros(size(tModel.mets));
vpH(contains(tModel.mets,'_c')) = R*T*log(1e-7);
vpH(contains(tModel.mets,'_e')) = R*T*log(10^(-6.3));
dpH = 0.7;
potential = 1e-3 * (33.33 * dpH - 143.33);                                 % V  = J/C
potential = 1e-3 * potential;                                              % kV = kJ/C

% dGr matrix
dGrMatrix = [R*T*St St Tt*(tModel.dGf0_H + (tModel.NH .* vpH)) + F*potential*tModel.charVec(~exchange)]; % Last column is constant
dGrRxns = tModel.rxns(~exchange);

%% Calculate dGr bounds
% Pre-allocate memory
nVars = size(dGrMatrix,2);
nTop  = size(feasibleTopologies,2);
dGrboundsTMFA   = zeros(size(dGrRxns,1),2,nTop);
dGrboundsFcTMFA = dGrboundsTMFA;
flagLBTMFA = dGrboundsTMFA; flagUBTMFA = dGrboundsTMFA;
XTMFA = zeros(nVars,2*nVars,nTop);
flagLBFcTMFA = dGrboundsTMFA; flagUBFcTMFA = dGrboundsTMFA;
XFcTMFA = zeros(nVars,2*nVars,nTop);

% dGr calculation
nRxns = size(tModel.rxns,1);
LPproblem.Aeq     = [];
LPproblem.beq     = [];
LPproblem.Aineq   = [];
LPproblem.bineq   = [];
LPproblem.options = optimoptions('linprog','Display','none');
LPproblem.solver  = 'linprog';
for i=1:size(feasibleTopologies,2)
    fluxDir = ones(nRxns,1);
    fluxDir(boundsTMFA(1:nRxns,2,i)<=1e-8) = -1;
    fluxDir(exchange) = [];
    tempDGRmatrix = fluxDir .* dGrMatrix;
    % Optimization problem
    % TMFA
    LPproblem.lb = [boundsTMFA(nRxns+1:end,1,i); 1];
    LPproblem.ub = [boundsTMFA(nRxns+1:end,2,i); 1];
    for j=1:size(tempDGRmatrix,1)
        LPproblem.f = tempDGRmatrix(j,:);
        [x, fval, flag] = linprog(LPproblem);
        flagLBTMFA(j,1,i) = flag;
        if flag > 0
            dGrboundsTMFA(j,1,i) = fval;
            XTMFA(:,2*j-1,i) = x;
        else
            disp("no feasible solution found")
            break
        end
        % Maximization
        LPproblem.f = -1*LPproblem.f;
        [x, ~, flag] = linprog(LPproblem);
        flagUBTMFA(j,2,i) = flag;
        if flag > 0
            dGrboundsTMFA(j,2,i) = fval;
            XTMFA(:,2*j,i) = x;
        else
            disp("no feasible solution found")
            break
        end
    end
    % tc-TMFA
    LPproblem.lb = [boundsFcTMFA(nRxns+1:end,1,i); 1];
    LPproblem.ub = [boundsFcTMFA(nRxns+1:end,2,i); 1];
    for j=1:size(tempDGRmatrix,1)
        LPproblem.f = tempDGRmatrix(j,:);
        [x, fval, flag] = linprog(LPproblem);
        flagLBTMFA(j,1,i) = flag;
        if flag > 0
            dGrboundsFcTMFA(j,1,i) = fval;
            XTMFA(:,2*j-1,i) = x;
        else
            disp("no feasible solution found")
            break
        end
        % Maximization
        LPproblem.f = -1*LPproblem.f;
        [x, ~, flag] = linprog(LPproblem);
        flagUBTMFA(j,2,i) = flag;
        if flag > 0
            dGrboundsFcTMFA(j,2,i) = fval;
            XTMFA(:,2*j,i) = x;
        else
            disp("no feasible solution found")
            break
        end
    end

end






%% Save bounds


end
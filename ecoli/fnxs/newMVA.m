function newModel = newMVA(model,tModel,relax)
% Performs Metabolite Variability Analysis. Updates model bounds and adds
% thermodinamic restrictions. Returns the updated model and 2*(n°
% variables) feasible points.

% Build Thermodynamic Model
LPproblem = newbuildTherModels(model,tModel,relax);

nVar         = size(LPproblem.lb,1);
X            = zeros(nVar,2*nVar);
flagLB       = zeros(nVar,1);
flagUB       = zeros(nVar,1);
newLB        = LPproblem.lb;
newUB        = LPproblem.ub;
LPproblem.f  = zeros(1,nVar);
dGrMatrix = LPproblem.dGrMatrix;
dGrVector = LPproblem.dGrVector;
LPproblem = rmfield(LPproblem,{'dGrMatrix' 'dGrVector'});
dGrBounds = zeros(size(dGrVector,1),2);
fact = 1;

for i=1:1:nVar
    LPproblem.f(i) = 1;
    [x, ~, flag] = linprog(LPproblem);
    flagLB(i) = flag;
    if flag > 0
        newLB(i) = x(i);
        X(:,2*i-1) = x;
    else
        fact = 0;
        break
    end
    LPproblem.f(i) = -1;
    [x, ~, flag] = linprog(LPproblem);
    flagUB(i) = flag;
    if flag > 0
        newUB(i) = x(i);
        X(:,2*i) = x;
    else
        fact = 0;
        break
    end
    LPproblem.f(i) = 0;
end

%% dGr Optimization

if fact
    for i=1:size(dGrMatrix,1)
        % Minimization
        LPproblem.f = dGrMatrix(i,:);
        [~, fval, flag] = linprog(LPproblem);
        if flag > 0
            dGrBounds(i,1) = fval;
        else
            disp("Problem with dGr optimization")
            break
        end
        % Maximization
        LPproblem.f = -1*dGrMatrix(i,:);
        [~, fval, flag] = linprog(LPproblem);
        if flag > 0
            dGrBounds(i,2) = -fval;
        else
            disp("Problem with dGr optimization")
            break
        end
    end
end

disp("MVA Complete!")

% Update model
newModel = tModel;
newModel.tlb    = newLB;
newModel.tub    = newUB;
newModel.lb     = newLB(1:size(newModel.S,2));
newModel.ub     = newUB(1:size(newModel.S,2));
newModel.lbc    = newLB(size(newModel.S,2)+1:end);
newModel.ubc    = newUB(size(newModel.S,2)+1:end);
newModel.flagLB = flagLB;
newModel.flagUB = flagUB;
newModel.dGrBounds = dGrBounds + dGrVector;
newModel.X      = X;

end
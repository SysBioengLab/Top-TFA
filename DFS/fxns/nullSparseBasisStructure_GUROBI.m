function modelSNP = nullSparseBasisStructure_GUROBI(modelSNP,x)
% Defines model structure for Fast-SNP
% Inputs:     model structure
%
% Outputs:    model structure for Fast-SNP (gurobi)
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2016 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter initialization
[m,n] = size(modelSNP.S);
epsilon = 1e-3;

% Definition of model structure
modelSNP.A = sparse([full(modelSNP.S),zeros(m,n);...
                    -eye(n),eye(n);...
                    eye(n),eye(n);...
                    x,zeros(1,n)]);

% Redefine LP problem:
% (1) RHS
modelSNP.rhs = zeros(m+2*n+1,1);
modelSNP.rhs(end-1) = epsilon;

% (2) Sense
modelSNP.sense = blanks(m+2*n+1);
for ix = 1:m+2*n+1
    if ix<=m
        modelSNP.sense(ix) = '=';
    else
        modelSNP.sense(ix) = '>';
    end
end

% (3) Bounds
modelSNP.lb = [modelSNP.lb;-1e3*ones(n,1)];
modelSNP.ub = [modelSNP.ub;1e3*ones(n,1)];

% (4) Model sense
modelSNP.modelsense = 'min';

% (5) Variables
modelSNP.vtype = 'C';

% (6) Obj fxn
modelSNP.obj = [zeros(n,1);ones(n,1)];
function findUnbalancedMets(model)
params.outputflag = 0;
fluxRef = 100;
flowCheck = fluxCheck(model);
fluxIdxs = 1:size(model.S,2);
flowCheck.lb(fluxIdxs) = -fluxRef*((sign(model.lb(fluxIdxs))<0));%+1e-1*(sign(model.lb)>=0))';
flowCheck.ub(fluxIdxs) = fluxRef*((sign(model.ub(fluxIdxs))>0));%-1e-1*(sign(model.ub)<=0))';

% Solve QP problem and update counter
sol = gurobi(flowCheck,params)
% delta = sol.x(size(model.S,2):end);
% Remove unfeasible flow patterns
if abs(sol.objval)>1e-9
    unbalancedMets = find(abs(sol.x(size(model.S,2)+1:end))>1e-5)
    [~,idxRxns] = find(model.S(unbalancedMets,:)~=0);
    model.rxns(idxRxns);
    sol.x(idxRxns);
end

function feasibilityProblem = fluxCheck(model)
% Extract relevant information
[m,n] = size(model.S);
K = 1e4;

% Problem definition
feasibilityProblem.modelsense = 'min';
feasibilityProblem.obj   = zeros(n+m,1);
feasibilityProblem.A     = sparse([model.S,eye(m)]); % Loop condition and slack condition
feasibilityProblem.Q     = sparse([zeros(n,m+n);zeros(m,n),eye(m)]);
feasibilityProblem.sense = '=';
feasibilityProblem.vtype = 'C';
feasibilityProblem.rhs   = zeros(m,1);
feasibilityProblem.lb    = -K*ones(m+n,1);
feasibilityProblem.ub    = K*ones(m+n,1);

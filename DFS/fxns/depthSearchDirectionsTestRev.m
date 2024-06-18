function depthSearchDirectionsTestRev(model,directionSpace,filename)
% Performs a depth-first search through the feasible space of directions
% using a branch and bound algorithm
% Inputs:  model structure
% Outputs: tfm (topological flux modes)
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2015 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of global variables
optimCounter = 0;
tfm          = [];
tol          = 1e-6;
K            = 1e4;

% Optimization and general parameters
params.outputflag = 0;

% Search parameters
reverRxns = find(model.rev==1);
dirTemp = zeros(1,model.numRxns);
depth   = length(reverRxns);
disp(['Tree depth: ',num2str(depth)])

% Define directionalities for each
for ix = 1:model.numRxns
    % reversible reactions can be either + or -
    if ~any(ix==reverRxns)
        if model.lb(ix)<0
            dirTemp(ix) = -1;
        elseif model.ub(ix)>0            
            dirTemp(ix) = 1;
        end
    end
end

% Initialize directionality problem
looplessProblem = looplessStructureLP(model);

% Perform depth-first search: recursive call
disp('Initializing tfm search...');
depthSearchTest(directionSpace,dirTemp,0);
disp(['Topological flux modes found = ',num2str(size(tfm,1))]);
clearvars -except tfm optimCounter model unfeasiblePattern filename
save([filename,'.mat']);

    function depthSearchTest(directionSpace,dirTemp,level)
        
        % Solve LP in the bottom
        if depth == level
            
            if ~mod(100*optimCounter,2^depth)
                fprintf('Progress status... %d\n',100*optimCounter/2^depth)
            end
            
            % Fix the boundaries of the zero, positive and negative reactions
            looplessProblem.lb(1:model.numRxns) = model.lb(1:model.numRxns).*(dirTemp==-1)'+tol*(dirTemp==1)';
            looplessProblem.ub(1:model.numRxns) = model.ub(1:model.numRxns).*(dirTemp==1)'-tol*(dirTemp==-1)';
            
            % Fix the boundaries of the loop variables (only for internal rxns)
            looplessProblem.lb(model.numRxns+model.internal) = -K*(1+sign(dirTemp(model.internal)))-sign(dirTemp(model.internal));
            looplessProblem.ub(model.numRxns+model.internal) = K*(1-sign(dirTemp(model.internal)))-sign(dirTemp(model.internal));
            
            % Solve and update optimization counter
            sol = gurobi(looplessProblem,params);
            optimCounter = optimCounter+1;
            
            % Check feasibility and save if tfm feasible
            if strcmp(sol.status,'OPTIMAL')
                tfm = [tfm;dirTemp];
            end
            
            % Recursive call
        else
            level = level+1;
            for j = 1:length(directionSpace{reverRxns(level)})
                dirTemp(reverRxns(level)) = directionSpace{reverRxns(level)}(j);
                depthSearchTest(directionSpace,dirTemp,level);
            end
        end
    end
end
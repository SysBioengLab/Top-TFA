function depthSearchUnfeasibleBalance(model,directionSpace,constIdx)
% Performs a depth-first search through the feasible space of directions
% using a branch and bound algorithm
% Inputs:  model structure
% Outputs: tfm (topological flux modes)
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2015 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of global variables
global unfeasiblePattern counter
counter = 0;
unfeasiblePattern = [];

% Optimization and general parameters
lbRef = model.lb(1:model.numRxns);
ubRef = model.ub(1:model.numRxns);
params.outputflag = 0;

% Initialize model structure
[~,activeFluxes] = find(model.S(constIdx,:)~=0);
tempModel.A   = sparse(model.S(constIdx,activeFluxes));
tempModel.rhs = 0;
tempModel.sense = '=';
tempModel.obj = zeros(length(activeFluxes),1);
tempModel.vtype = 'C';

% Search parameters
depth    = length(directionSpace(activeFluxes));
dirTemp  = 2*ones(1,depth);

% Perform depth-first search: recursive call
depthSearchUnbalance(directionSpace(activeFluxes),dirTemp,0)
disp(['Unfeasible modes found in balance ',num2str(constIdx),' = ',num2str(size(unfeasiblePattern,1))])
clearvars -except counter model unfeasiblePattern constIdx directionSpace
save(['unfeasibleBalance_',num2str(constIdx),'.mat'])

    function depthSearchUnbalance(directionSpace,dirTemp,level)
        % Solve LP problem at the bottom node
        if depth == level                   
                        
            % Fix the boundaries of the zero, positive and negative reactions
            tempModel.lb = lbRef(activeFluxes).*(dirTemp==-1)'+1e-1*(dirTemp==1)';
            tempModel.ub = ubRef(activeFluxes).*(dirTemp==1)'-1e-1*(dirTemp==-1)';           
             
            % Solve and update optimization counter
            sol = gurobi(tempModel,params);
            counter = counter + 1;
            
            % Check feasibility and save if tfm feasible
            if ~strcmp(sol.status,'OPTIMAL')
                tempPattern = 2*ones(1,model.numRxns);
                tempPattern(activeFluxes) = dirTemp;
                unfeasiblePattern = [unfeasiblePattern;tempPattern];
            end
    
            % Recursive call
        else
            level = level+1;
            for j = 1:length(directionSpace{level})
                dirTemp(level) = directionSpace{level}(j);
                depthSearchUnbalance(directionSpace,dirTemp,level);
            end
        end
    end
end
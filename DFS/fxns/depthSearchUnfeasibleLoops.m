function depthSearchUnfeasibleLoops(model,directionSpace,constIdx)
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2015 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of global variables
global unfeasiblePattern counter
counter = 0;
unfeasiblePattern = [];

% Optimization and general parameters
params.outputflag = 0;

% Initialize model structure
[~,activeFluxes] = find(model.Nint(constIdx,:)~=0);
tempModel.A   = sparse(model.Nint(constIdx,activeFluxes));
tempModel.rhs = 0;
tempModel.sense = '=';
tempModel.obj = zeros(length(activeFluxes),1);
tempModel.vtype = 'C';

% Search parameters
depth    = length(directionSpace(activeFluxes));
dirTemp  = 2*ones(1,depth);

% Perform depth-first search: recursive call
depthSearchLoops(directionSpace(activeFluxes),dirTemp,0)
disp(['Unfeasible modes found in loop ',num2str(constIdx),' = ',num2str(size(unfeasiblePattern,1))])
clearvars -except counter model unfeasiblePattern constIdx directionSpace
save(['unfeasibleLoop_',num2str(constIdx),'.mat'])

    function depthSearchLoops(directionSpace,dirTemp,level)
        % Solve LP problem at the bottom node
        if depth == level                   
                        
            % Fix the boundaries of the loop variables            
            tempModel.lb = -1e4*(1+sign(dirTemp'))-sign(dirTemp');
            tempModel.ub = 1e4*(1-sign(dirTemp'))-sign(dirTemp');      
             
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
                depthSearchLoops(directionSpace,dirTemp,level);
            end
        end
    end
end
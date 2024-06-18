function depthSearchUnfeasibleBalanceRev(model,directionSpace,constIdx,Stoic,folderName)
% Finds unfeasible flux patterns based on mass balance constraints
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2015 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of global variables
unfeasiblePatterns    = [];
counter = 0;

% Check if we have both active and rev reactions for this constraint
[~,activeFluxes] = find(Stoic~=0);
reverRxns        = find(model.rev==1);
compoundVector   = [activeFluxes(:);reverRxns(:)];
if length(unique(compoundVector))==length(compoundVector)
    return;
else
    % Determine the reversible rxns belonging to the active set
    activeRxns = [];
    dirTemp  = zeros(1,model.numRxns);
    for ix = activeFluxes
        if any(ix==reverRxns)
            activeRxns = [activeRxns,ix];
        else
            dirTemp(ix) = directionSpace{ix};
        end
    end
end

% Search parameters
depth = length(activeRxns);

% Perform depth-first search: recursive call
depthSearchUnbalance(directionSpace,dirTemp,0);

if ~isempty(unfeasiblePatterns)
    disp(['Unfeasible modes found in balance ',num2str(constIdx),' = ',num2str(size(unfeasiblePatterns,1))]);
    clearvars -except counter model constIdx directionSpace unfeasiblePatterns folderName
    save([folderName,'/unbalancedTopology_',num2str(constIdx),'.mat']);
end
    function depthSearchUnbalance(directionSpace,dirTemp,level)
        if depth == level
            counter = counter+1;
            
            % Check mass balances topologicaly
            dirProd = Stoic(activeFluxes).*dirTemp(activeFluxes);
            if all(dirProd>0) || all(dirProd<0)
                tempPattern = zeros(1,model.numRxns);
                tempPattern(activeRxns) = dirTemp(activeRxns);
                unfeasiblePatterns = [unfeasiblePatterns;tempPattern];
            end
            
            % Recursive call
        else
            level = level+1;
            for j = 1:length(directionSpace{activeRxns(level)})
                dirTemp(activeRxns(level)) = directionSpace{activeRxns(level)}(j);
                depthSearchUnbalance(directionSpace,dirTemp,level);
            end
        end
    end
end
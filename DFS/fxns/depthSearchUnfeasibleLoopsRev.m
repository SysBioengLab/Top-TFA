function depthSearchUnfeasibleLoopsRev(model,directionSpace,constIdx,Nint,folderName)
% Finds unfeasible flux patterns based on loopless constraint
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2015 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of global variables
counter = 0;
unfeasiblePatterns = [];

% Check if we have both active and rev reactions for this constraint
[~,activeFluxes] = find(Nint~=0);
reverRxns        = find(model.rev==1);
compoundVector   = [activeFluxes(:);reverRxns(:)];
if length(unique(compoundVector))==length(compoundVector);
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
depthSearchLoops(directionSpace,dirTemp,0);

if ~isempty(unfeasiblePatterns)
    disp(['Unfeasible modes found in loop ',num2str(constIdx),' = ',num2str(size(unfeasiblePatterns,1))]);
    clearvars -except counter model constIdx directionSpace unfeasiblePatterns folderName
    save([folderName,'/loopyTopology_',num2str(constIdx),'.mat']);
end

    function depthSearchLoops(directionSpace,dirTemp,level)
        if depth == level
            counter = counter+1;
            
            % Check loop condition topologicaly
            dirProd = sign(Nint(activeFluxes).*dirTemp(activeFluxes));
            if all(dirProd==1) || all(dirProd==-1)
                tempPattern = zeros(1,model.numRxns);
                tempPattern(activeRxns) = dirTemp(activeRxns);
                unfeasiblePatterns = [unfeasiblePatterns;tempPattern];
            end
            
            % Recursive call
        else
            level = level+1;
            for j = 1:length(directionSpace{activeRxns(level)})
                dirTemp(activeRxns(level)) = directionSpace{activeRxns(level)}(j);
                depthSearchLoops(directionSpace,dirTemp,level);
            end
        end
    end
end
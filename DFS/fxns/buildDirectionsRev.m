function [model,directionSpace] = buildDirectionsRev(model)
% Extract important information
model.metsIdxs  = model.numRxns+length(model.internal)+(1:size(model.S,1));
model.loopsIdxs = model.numRxns+length(model.internal)+size(model.S,1)+(1:size(model.Nint,1));

% Find reaction nature
reverRxns = find(model.rev==1);

% Define directionalities for each
for ix = 1:model.numRxns
    % reversible reactions can be either + or -
    if sum(reverRxns==ix)
        directionSpace{ix} = [-1,1];
        
        % irreversible reactions have one possibility
    else
        if model.lb(ix)<0
            directionSpace{ix} = -1;
        elseif model.ub(ix)>0
            directionSpace{ix} = 1;
        end
    end
end
disp('Directional space built!')
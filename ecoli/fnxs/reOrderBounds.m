function bounds = reOrderBounds(prevBounds,feasible)
% Change the format of the bounds. Only includes feasible topologies.

if size(feasible,1) == 1
    if feasible
        bounds = prevBounds;
    else
        bounds = [];
    end
else
    prevBounds(:,:,~feasible) = [];
    bounds = zeros(size(prevBounds,1),2*size(prevBounds,3));
    for i=1:size(prevBounds,3)
        bounds(:,[2*i-1 2*i]) = prevBounds(:,:,i);
    end
end
end
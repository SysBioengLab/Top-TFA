function rDF = rangeDifference(newBounds,prevBounds,nBf,pBf,model,tol)
% Calculates the range difference between two bounds. If second range is
% zero, zero is returned. newBounds must be at least equally restrictive
% than prevBounds.

% prevBounds selection
if size(nBf,1) > 1
    if ~isequal(nBf,pBf)
        nBf(~pBf) = [];
        toErase = find(~nBf);
        nToErase = zeros(1,2*size(toErase,1));
        i=1;
        for n=toErase'
            nToErase(i) = 2*n-1;
            i = i+1;
            nToErase(i) = 2*n;
            i = i+1;
        end
        prevBounds(:,nToErase) = [];
    end
end

% Concentration transformation
nRxns = size(model.lb,1);
nMets = size(model.lbc,1);
prevBounds(nRxns+1:nRxns+nMets,:) = exp(prevBounds(nRxns+1:nRxns+nMets,:));
newBounds(nRxns+1:nRxns+nMets,:)  = exp(newBounds(nRxns+1:nRxns+nMets,:));

% Calculus
rDF = zeros(size(newBounds,1),size(newBounds,2)/2);
for i=1:size(rDF,2)
    rDF(:,i) = 1 - abs((newBounds(:,2*i-1)-newBounds(:,2*i))./(prevBounds(:,2*i-1)-prevBounds(:,2*i)));
end
rDF(isinf(rDF)|isnan(rDF)) = 0;
rDF(abs(rDF)<tol) = 0;

end
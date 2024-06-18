function model = parseInternalRxns(model,numRxns)
internalRxns = [];
exchangeRxns = [];
for ix = 1:numRxns
    if sum(full(model.S(:,ix)~=0)) > 1
        internalRxns = [internalRxns,ix];
    else
        exchangeRxns = [exchangeRxns,ix];
    end
end

% Reorganize reaction order
model.S        = model.S(:,[internalRxns,exchangeRxns]);
model.rxns     = model.rxns([internalRxns,exchangeRxns]);
model.rxnNames = model.rxnNames([internalRxns,exchangeRxns]);
model.lb       = model.lb([internalRxns,exchangeRxns]);
model.ub       = model.ub([internalRxns,exchangeRxns]);
model.internal = 1:length(internalRxns);
model.exchange = length(internalRxns)+1:length(internalRxns)+length(exchangeRxns);

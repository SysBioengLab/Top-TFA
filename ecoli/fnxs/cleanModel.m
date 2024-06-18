function model = cleanModel(model)
% Performs a Flux-Variability Analysis to remove all blocked rxns and
% metabolites. It also corrects the direction of irreversible rxns.

% model = parseInternalRxns(model);

% FVA
% [minFlux,maxFlux] = FVA(model);          % 10/03: Usaremos la funcion de
% cobra al 5% del valor máximo del objetivo
[minFlux,maxFlux] = fluxVariability(model,5);

% Update model
if ~isfield(model,'vTol')
    model.vTol = 1e-8;
end
minFlux(abs(minFlux)<model.vTol) = 0;
maxFlux(abs(maxFlux)<model.vTol) = 0;
model.lb = minFlux;
model.ub = maxFlux;

% Compress model
model = compressModel(model);
% model = parseInternalRxns(model);  % Omitir por ahora

% MVA (TMFA)
model = MVA(model);
model = compressModel(model);

% Update thermo, transport and exchange rxns
model.exchange_rxns = findExcRxns(model);
model.thermo_rxns = ~model.exchange_rxns;
ext = logical(sum(abs(model.S) .* model.comp,1));
int = logical(sum(abs(model.S) .* ~model.comp,1));
model.transport_rxns = transpose(logical(int & ext));                      % Rxns that involve intra and extra celular metabolites

end
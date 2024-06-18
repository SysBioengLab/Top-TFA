function sample = addFieldsToSample(sample,model,zeroRxns)
% Traspases important fields from model to sample

% Common fields
sample.uniform        = model.uniform;
sample.discreteFactor = model.discreteFactor;

% BuildProbParams fields
% Find indexes to delete from structures
indxTodel = [];
conc_indxTodel = [];
for i=zeroRxns
    if model.normalVars(i) == 1
        indxTodel(end+1) = sum(model.normalVars(1:i));
    end
    if ~isequal(model.concNormal,0)
        if model.concNormal(i) == 1
            conc_indxTodel(end+1) = sum(model.concNormal(1:i));
        end
    end
end
% Delete indexes
model.normalVars(zeroRxns) = [];
model.bigMu(indxTodel) = [];
model.bigCov(indxTodel,:) = [];
model.bigCov(:,indxTodel) = [];

sample.normalVars = model.normalVars;
sample.bigMu      = model.bigMu;
sample.bigCov     = model.bigCov;

if model.concNormal == 0
    sample.concNormal = 0;
else
    model.concNormal(zeroRxns)   = [];
    model.concMu(conc_indxTodel) = [];
    model.concStd(conc_indxTodel) = [];
    sample.concNormal = model.concNormal;
    sample.concMu     = model.concMu;
    sample.concStd    = model.concStd;
end

% Fxns
sample.isFeasible = @(points) closa_constrainRevision(points,sample.S,sample.lb,sample.ub,sample.feasTol);
sample.densityFxn = @(points) customDensityFxn(points,sample);             % Define fxn that evaluates de pdf
% sample.densityFxn = @(points) testDensityFxn(points);             % test

end
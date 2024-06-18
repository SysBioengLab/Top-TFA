function sample = buildProbParams(model,sample)
% Multi-Variate Normal (dGf)
sample.bigMu      = model.dGft_mu;
sample.bigCov     = model.dGft_cov;
sample.normalVars = contains(sample.vars,'dGft');
sample.normalVars = logical(sample.normalVars - contains(sample.vars,'dGft_c242'));

% exp(log(concentrations)) ~ Gaussian
if ~isempty(model.concNormalMets)
    sample.concMu     = model.concMu;
    sample.concStd    = model.concStd;
    sample.concNormal = ismember(sample.vars,model.concNormalMets);
else
    sample.concNormal = 0;
end
end
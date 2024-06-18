function sample = buildProbParams(tModel,sample)

% Multi-Variate Normal (dGf)
sample.bigMu      = tModel.dGft_mu;
sample.bigCov     = tModel.dGft_cov;
sample.normalVars = contains(sample.vars,'dGft');
sample.normalVars = logical(sample.normalVars - contains(sample.vars,'dGft_h_'));

% exp(log(concentrations)) ~ Gaussian
if ~isempty(tModel.concNormalMets)
    metInModel = ismember(tModel.concNormalMets,tModel.mets); % Check mets in model
    sample.concMu     = tModel.concMu(metInModel);
    sample.concStd    = tModel.concStd(metInModel);
    sample.concNormal = ismember(sample.vars,tModel.concNormalMets(metInModel));
else
    sample.concNormal = 0;
end
end
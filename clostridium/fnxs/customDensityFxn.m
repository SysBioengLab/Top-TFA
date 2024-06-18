function normF = customDensityFxn(points,sample)
% Returns the probability vector of points. Importantly, the points are
% adjusted so the top probability equals one for prob1.

% MVN
prob1 = transformed_mvnpdf(points(sample.normalVars,:)',sample.bigMu',sample.bigCov);

% Meditions
if sum(sample.concNormal) > 0
    prob2 = mvnpdf(exp(points(sample.concNormal,:))',sample.concMu', sample.concStd');
    prob2 = prob2/sum(prob2);
else
    prob2 = 1;
end

% Normalization
prob3 = prob1 .* prob2;
normF = prob3/sum(prob3);

end
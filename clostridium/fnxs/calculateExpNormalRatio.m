function ratioExp = calculateExpNormalRatio(sample,prevPoint,nextPoint)
% Calculates the ratio of probabilities between 2 points. Returns
% log(probability).

if sum(sample.concNormal) > 0
    x1    = exp(prevPoint(sample.concNormal));
    x2    = exp(nextPoint(sample.concNormal));
    mu    = sample.concMu;
    sigma = sample.concStd;
    allratiosExp = 0.5 * ((x1 - mu)./sigma).^2 -0.5 * ((x2 - mu)./sigma).^2;
    ratioExp = sum(allratiosExp);
else
    ratioExp = 0; % log(1) = 0
end
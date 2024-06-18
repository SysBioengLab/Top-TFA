function [probability] = modelMVN(point,bigMu,bigCov,normalVars)
% Aplies mvnpdf function to the normal distributed variables.
probability = mvnpdf(point(normalVars),bigMu,bigCov);
end
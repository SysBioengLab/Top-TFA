function normF = testDensityFxn(points)
% Returns the probability vector of points. Importantly, the points are
% adjusted so the top probability equals one for prob1.

% MVN
prob1 = transformed_mvnpdf(points([83 84],:)',[-2351 -1478],[2.2 1.7; 1.7 1.5]);
% prob1 = normpdf(points(7,:)',9,0.1);
normF = prob1/sum(prob1);

end
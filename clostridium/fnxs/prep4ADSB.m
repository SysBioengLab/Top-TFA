function sample = prep4ADSB(sample)
% Preprocess of sample to use ADSB

sample.loopless = 0;
sample.numRxns = sample.numVars;
prevPoints = sample.warmUpPoints;

% Fold number of particles in the population relative to dim(Omega) (only for ADSB, default 3)
if ~isfield(sample,'populationScale'); sample.populationScale = 3; end

% Determine size of the feasible space
omegaSize = rank(sample.warmUpPoints);
if (sample.populationScale == 1); nDim = sample.populationScale*omegaSize+1;
else nDim = sample.populationScale*omegaSize; end
% Define number of chains and points per chain (at least three points are required)
minPoints             = 3;
sample.pointsPerChain = max([nDim,minPoints]);
sample.numChains      = ceil(sample.numSamples/sample.pointsPerChain);

sample.points = zeros(sample.numRxns,nDim,sample.numChains);
for ix = 1:sample.numChains
    sample.points(:,:,ix) = prevPoints(:,sampleSet(size(prevPoints,2),nDim,'false'));
end

% Define fxn to ensure samples are kept within the bounds
sample.keepWithinBounds = @(points) bringToBoundary(points,sample.lb,sample.ub);

% Define fxn that evaluates de pdf
sample.densityFxn = @(points) customDensityFxn(points,sample);

end
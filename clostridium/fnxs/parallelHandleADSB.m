function sample = parallelHandleADSB(sample)
% Parallel run of ADSB

fprintf('Sampling in progress (multiple cores)...\n');

% Delete active workers
if ~isempty(gcp('nocreate')); delete(gcp('nocreate')); end        

% Try to use the number of workers requested. If not posible, use the number of workers allocated to the 'local' profile
if ~isempty(sample.numCores)
    try
        parpool(sample.numCores);
    catch
        parWorkers = parpool('local');
        sample.numCores = parWorkers.NumWorkers;
    end
else
    parWorkers = parpool('local');
    sample.numCores = parWorkers.NumWorkers;
end

% Replicate sample structure for parallel sampling
chainsPerWorker = fix(sample.numChains/sample.numCores);
chainsPerWorker = chainsPerWorker*ones(1,sample.numCores-1);
chainsPerWorker = [chainsPerWorker,sample.numChains-sum(chainsPerWorker)];

% Assign appropriate number of structures to each core
idx1 = 0; idx2 = 0;
for ix = 1:sample.numCores            
    idx1 = idx1+1;                                 % Update indices
    idx2 = idx2+chainsPerWorker(ix);
    samples{ix}           = sample;
    samples{ix}.numChains = chainsPerWorker(ix);
    samples{ix}.points    = sample.points(:,:,idx1:idx2);
    idx1 = idx2;                           % Update starting points
end

% Perform parallel sampling using ADSB
parfor workerIdx = 1:sample.numCores
    rng('shuffle');
    [samples{workerIdx}.points,samples{workerIdx}.samplingTime] = nonUniformADSB(samples{workerIdx},workerIdx);
end
delete(gcp('nocreate'));

% Build definitive structure
sample.points       = [];
sample.samplingTime = 0;
for ix = 1:sample.numCores
    sample.points = [sample.points,reshape(samples{ix}.points,...
                            sample.numRxns,sample.pointsPerChain*chainsPerWorker(ix))];
    sample.samplingTime = max([sample.samplingTime,samples{ix}.samplingTime]);
end

% Fill with the exact number of samples
sample.points = sample.points(:,1:sample.numSamples);

end
% Compute MCMC diagnostics
% -------------------- Copyright (C) 2023 Pedro A. Saa --------------------
clearvars,clc

disp("Running diagnostics from Top-TFA\clostridium...")
filename='results';
fileList = {'samples\bigUni.mat','samples\bigNouni.mat'};

for ix = 1:2
    fileTemp = fileList{ix};  
    load(fileTemp)
    if ix == 1
        sample = bigUni;
    else
        sample = bigNouni;
    end
    sample.numSamples = size(sample.points,2);

    % Running diagnostics: Split samples resulting from the chain in 10 segments for
    % convergence analysis
    numChainsForDiagnostics      = 10;
    pointsPerChainForDiagnostics = fix(sample.numSamples/numChainsForDiagnostics);
    sample.points = reshape(sample.points,sample.numRxns,pointsPerChainForDiagnostics,numChainsForDiagnostics);
    
    % Calculate potential scale reduction statistics and recover original chain
    sample.points = permute(sample.points,[2,1,3]);
    [sample.R,sample.Rint,sample.Neff,sample.tau,sample.thin] = psrf(sample.points);
    sample.points = permute(sample.points,[2,3,1]);
    sample.points = reshape(sample.points,sample.numRxns,pointsPerChainForDiagnostics*numChainsForDiagnostics);
    sample.points = sample.points(:,1:sample.numSamples);
    
    % Calculate chains statistics
    sample.mu     = mean(sample.points,2);
    sample.sigma  = std(sample.points')';
    sample.median = median(sample.points,2);
    sample.IQR    = iqr(sample.points,2);
    
    % Update file
    save(fileTemp,'sample');

end
function tcTMFAsampler(bigModel,bounds,volProp,name,options)
% Performs random sampling of metabolic models
% Based on LooplessFluxSampler (Saa 2019).
% 
% USAGE:
%              sample = runSampling(model, options)
%
% INPUTS:
%              model (structure):    (the following fields are required - others can be supplied)
%                                    * S  - `m x 1` Stoichiometric matrix
%                                    * lb - `n x 1` Lower bounds
%                                    * ub - `n x 1` Upper bounds
%                                    * rxns - `n x 1` rxn identifiers (cell array)
%              options (structure):  (the following fields are required - others can be supplied)
%                                    * numSamples - number of points (double)
%
% OPTIONAL INPUTS:
%              options (structure):  (the following fields are optional)
%                                    * numDiscarded - Burn-in (double) (only used in ll-ACHRB) (default 0)
%                                    * stepsPerPoint - Thinning or number of steps per effective point (double). In ADSB, it refers to the expected number of times a point is moved 
%                                    * algorithm - 'CHRR' (default), metropolisHRB, kernelHRB, generalHR_kernel
%                                    * vTol - Numerical flux tolerance (default 1e-8)
%                                    * diagnostics - MCMC diagnostics flag, true (default) or false
%
% OUTPUT:
%              sample (structure):   sampling structure containing #numSamples random points
%
% -------------------- Copyright (C) 2023 Pedro A. Saa --------------------

% Check inputs
if (nargin<5)
    fprintf('Not enough input arguments.');
    return;
else
    % Check if COBRA toolbox is in the path
    if ~exist('initCobraToolbox','file')
        fprintf('COBRA Toolbox is not in the path. Please refer to https://github.com/opencobra/cobratoolbox for more information on how to install this toolbox.');
        return;
    end

    % Number of samples
    if ~isfield(options,'numSamples')
        fprintf('Field numSamples has not been defined.'); return; end

    % Samples discarded or burn-in (default 0)
    if isfield(options,'numDiscarded'); sample.numDiscarded = options.numDiscarded;
    else sample.numDiscarded = 0; end

    % Steps per point
    if isfield(options,'stepsPerPoint'); sample.stepsPerPoint = options.stepsPerPoint;
    else sample.stepsPerPoint = 1e2; end

    % Sampling algorithm 'HRB' (default) or 'll_ACHRB'
    if isfield(options,'algorithm'); sample.algorithm = options.algorithm;
    else sample.algorithm = 'CHRR'; end

    % Numerical flux tolerance  (default 1e-8)
    if isfield(options,'vTol'); sample.vTol = options.vTol;
    else sample.vTol = 1e-8; end
    
    % Previous centroid (default empty)
    if isfield(options,'prevCentroid'); sample.prevCentroid = options.prevCentroid;
    else sample.prevCentroid = []; end

    % tcTMFA or TMFA (default tcTMFA)
    if isfield(options,'TMFA'); TMFA = options.TMFA; else TMFA = 0; end

    % Run MCMC diagnostics
    if isfield(options,'diagnostics'); sample.diagnostics = options.diagnostics;
    else sample.diagnostics = 1; end
end

% Define solver parameters
changeCobraSolverParams('LP','optTol',1e-9);
changeCobraSolverParams('LP','feasTol',1e-9);

% Initialization
sample.loopless = 0;                                 % EDITAR

if ~(size(bounds,3) == size(volProp,1))
    fprintf('Error. Size of bounds and volProp are incompatible')
    return
else
    nT = size(bounds,3);                             % Number of topologies
end

% Define fields to clear on each iteration
if sample.diagnostics
    otherFields   = {'R','Rint','Neff','tau','thin','mu','sigma',...
            'median','IQR'}';
else
    otherFields = {};
end
fields = {'vars','numVars','S','lb','ub','centroid',...
        'points','samplingTime','rejectionRate','numSamples'}';
fields = [fields; otherFields];

% Define numSamples per topology
nSamples = ceil(options.numSamples*volProp);

% Define cumulative rejection rate
cumRR = 0;

% Main loop
for i=1:nT
    %% I. Preparation phase
    % Time model pre-processing and sampling preparation
%     t0 = cputime;
    fprintf('Working with topology %d of %d.\n',i,nT);
    % Define number of samples
    sample.numSamples = nSamples(i);
    % Define sampling matrix and other fields
    matrix_struct    = createSamplingMatrix(bigModel,bounds(:,:,i),TMFA);
    sample.dGrMatrix = matrix_struct.dGrMatrix;
    sample.dGrRxns   = matrix_struct.dGrRxns;
    % Homogenize system
    matrix_struct = homogenizer(matrix_struct);
    % Run FVA to ensure bounds are correctly determined
    matrix_struct = updateBounds(matrix_struct);

    % Define fields related to the model in sample structure
    sample.vars     = matrix_struct.vars;
    sample.numVars  = numel(matrix_struct.lb);
    sample.S        = matrix_struct.S;
    sample.lb       = matrix_struct.lb;
    sample.ub       = matrix_struct.ub;
    sample.b        = matrix_struct.b;
    sample.uTol     = 1e-10;                                               % Direction tolerance
    sample.bTol     = 1e-8;                                                % Minimum allowed distance to closest constraint
    sample = buildProbParams(bigModel.thermo,sample);

    % Define function to check feasibility of proposed points
    sample.isFeasible = @(points) constrainRevision(points,sample.S,sample.lb,sample.ub,sample.bTol*10);

     % Calculate centroid
    if isempty(sample.prevCentroid)
        % Generate initial warmup points
        [sample,sample.warmUpPoints] = warmupFcTMFA(sample);
        centroid = mean(sample.warmUpPoints,2);
        centroid(abs(centroid)<sample.vTol) = 0;
        sample.centroid = centroid;
        fprintf('%d warmup points spanning the feasible space have been defined.\n',size(sample.warmUpPoints,2));
    else
        sample.centroid = sample.prevCentroid(:,i);
    end
    
    clearvars matrix_struct centroid
    
    %% II. Sampling
    % Run appropriate sampler in either single core or parallel mode
    if strcmp(sample.algorithm,'CHRR')
        fprintf('Sampling in progress (single core)...\n--------------------------------------------\n');
        sample.c = zeros(size(sample.lb)); sample.rxns = sample.vars;
        startCHRR = tic();
        [sample.points, sample.roundedPolytope] = chrrSampler(sample,sample.stepsPerPoint,sample.numSamples,1);
        sample.samplingTime = toc(startCHRR);
        sample.rejectionRate = sum(~(sample.isFeasible(sample.points)))/sample.numSamples;

    elseif strcmp(sample.algorithm,'metropolisHRB')
        fprintf('Sampling in progress (single core)...\n--------------------------------------------\n');
        sample.generalHR = 0;
        [sample.points,sample.samplingTime,sample.rejectionRate] = metropolisHitAndRun(sample,sample.numSamples,sample.numDiscarded,sample.stepsPerPoint);

    elseif strcmp(sample.algorithm,'generalHR')
        fprintf('Sampling in progress (single core)...\n--------------------------------------------\n');
        sample.generalHR = 1;
        [sample.points,sample.samplingTime,sample.rejectionRate] = metropolisHitAndRun(sample,sample.numSamples,sample.numDiscarded,sample.stepsPerPoint);

    elseif strcmp(sample.algorithm,'kernelHRB')
        fprintf('Sampling in progress (single core)...\n--------------------------------------------\n');
        [sample.points,sample.samplingTime,sample.rejectionRate] = kernelHitAndRun(sample,sample.numSamples,sample.numDiscarded,sample.stepsPerPoint);
    
    elseif strcmp(sample.algorithm,'generalHR_kernel')
        fprintf('Sampling in progress (single core)...\n--------------------------------------------\n');
        [sample.points,sample.samplingTime,sample.rejectionRate] = generalHR_kernel(sample);
    
    elseif strcmp(sample.algorithm,'HRB')
        fprintf('Sampling in progress (single core)...\n--------------------------------------------\n');
        [sample.points,sample.samplingTime,sample.rejectionRate] = HRB(sample,sample.numSamples,sample.numDiscarded,sample.stepsPerPoint);
    end
    
    %%  IV. Run MCMC diagnostics
    if sample.diagnostics
        fprintf('Starting MCMC diagnostics...\n');
        % Split samples resulting from the chain in 10 segments for
        % convergence analysis of ll-ACHRB (this is not necessary for ADSB)
        sample.numChains = 10;
        sample.pointsPerChain = fix(sample.numSamples/sample.numChains);
        sample.points = reshape(sample.points,sample.numVars,sample.pointsPerChain,sample.numChains);
        
        % Calculate potential scale reduction statistics and recover original chain
        sample.points = permute(sample.points,[2,1,3]);
        [sample.R,sample.Rint,sample.Neff,sample.tau,sample.thin] = psrf(sample.points);
        sample.points = permute(sample.points,[2,3,1]);
        sample.points = reshape(sample.points,sample.numVars,sample.numChains*sample.pointsPerChain);
        sample.points = sample.points(:,1:sample.numSamples);
        
        % Calculate chains statistics
        sample.mu     = mean(sample.points,2);
        sample.sigma  = std(sample.points')';
        sample.median = median(sample.points,2);
        sample.IQR    = iqr(sample.points,2);

    else
        sample.points = sample.points(:,1:sample.numSamples);
    end
    % Save and reset sample struct for next iteration WARNING WITH PATH
    cumRR = cumRR + sample.rejectionRate;
    save(strcat("samples\sample_",name,"_top_",num2str(i),".mat"),"sample")
    sample = rmfield(sample,fields);

end

fprintf('\n Sampling Finalized! \n Total Rejection Rate: %d \n',cumRR/nT)

end
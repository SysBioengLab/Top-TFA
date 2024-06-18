function [points,samplingTime,rejectionRate] = alpha_sampling(model,alpha_dim,alpha_options)

%% I. Check inputs
% [COMMENT] QUITÉ LA OPCIÓN LOOPLESS
% Number of samples
if isfield(alpha_options,'numSamples'); sample.numSamples = alpha_options.numSamples;
else fprintf('Field numSamples has not been defined.'); return; end

% Samples discarded or burn-in (default 0)
if isfield(alpha_options,'numDiscarded'); sample.numDiscarded = alpha_options.numDiscarded;
else sample.numDiscarded = 0; end

% Steps per point (default 1e2)
if isfield(alpha_options,'stepsPerPoint'); sample.stepsPerPoint = alpha_options.stepsPerPoint;
else sample.stepsPerPoint = 1e2; end

% Sampling algorithm 'HRB' [COMMENT] (¿USAMOS HRB?)
if isfield(alpha_options,'algorithm'); sample.algorithm = alpha_options.algorithm;
else sample.algorithm = 'ADSB'; end

% Numerical flux tolerance  (default 1e-8)
if isfield(alpha_options,'vTol'); sample.vTol = alpha_options.vTol;
else sample.vTol = 1e-8; end

% Parallel sampling option (default false or 0)
if isfield(alpha_options,'parallelFlag'); sample.parallelFlag = alpha_options.parallelFlag;
else sample.parallelFlag = 0; end

% Number of cores for parallel sampling (default empty)
if isfield(alpha_options,'numCores'); sample.numCores = alpha_options.numCores;
else sample.numCores = []; end

% Run MCMC diagnostics [COMMENT] NO SÉ SI SE PUEDEN REALIZAR ESTOS DIAGNÓSTICOS
if isfield(alpha_options,'diagnostics'); sample.diagnostics = alpha_options.diagnostics;
else sample.diagnostics = 1; end

% Fold number of particles in the population relative to dim(Omega) (only
% for ADSB, default 1) [COMMENT] ¿QUITAR?
if strcmp(sample.algorithm,'ADSB') && isfield(alpha_options,'populationScale'); sample.populationScale = alpha_options.populationScale;
else sample.populationScale = 3; end

% Restart sampler from previous results [COMMENT] ¿QUITAR?
if ~isfield(model,'points'); sample.points = [];
else sample.points = model.points; end

% Define solver parameters [COMMENT] ¿CORREPONDE ESTA PARTE?
changeCobraSolverParams('MILP','intTol',1e-9);
changeCobraSolverParams('MILP','relMipGapTol',1e-12);
changeCobraSolverParams('LP','optTol',1e-9);
changeCobraSolverParams('LP','feasTol',1e-9);

%% II. Preparation phase
% Define fields related to the model in sample structure
sample.rxns     = model.rxns;
sample.numRxns  = numel(model.lb);
sample.S        = model.S;
sample.lb       = model.lb;
sample.ub       = model.ub;
sample.intRxns  = model.intRxns;
sample.exchRxns = model.exchRxns;
sample.uTol     = 1e-10;                                                   % Direction tolerance
sample.bTol     = 1e-10;                                                   % Minimum allowed distance to closest constraint
% [COMMENT] QUITÉ LA PARTE LOOPLESS

% Define fxn to ensure samples are kept within the bounds
sample.keepWithinBounds = @(points) bringToBoundary(points,sample.lb,sample.ub);

% Generate initial warmup points
sample.warmUpPoints = warmupLooplessACHRB(sample);

% Calculate centroid
centroid = mean(sample.warmUpPoints,2);
centroid(abs(centroid)<sample.vTol) = 0;

% % Check feasibility of the initial point. If centroid infeasible, find closest feasible point
% if sample.loopless && ~sample.isFeasible(centroid) % [COMMENT] ¿PUEDO QUITAR LA CONDICIÓN DE QUE SEA LOOPLESS?
%     sample.centroid = findNearestFeasiblePoint(sample,centroid,Nint',sample.vTol);
% else
%     sample.centroid = centroid;
% end

sample.centroid = centroid;

fprintf('%d warmup points spanning the feasible space have been defined.\n',size(sample.warmUpPoints,2));
clearvars -except sample



% copiar HRB




end



% %% I. Check inputs
% % [COMMENT] QUITÉ LA OPCIÓN LOOPLESS
% if (nargin<2)
%     fprintf('Not enough input arguments.');
%     return;
% else
%     % Check if COBRA toolbox is in the path
%     if ~exist('initCobraToolbox','file')
%         fprintf('COBRA Toolbox is not in the path. Please refer to https://github.com/opencobra/cobratoolbox for more information on how to install this toolbox.');
%         return;
%     end
% 
%     % Number of samples
%     if isfield(alpha_options,'numSamples'); sample.numSamples = alpha_options.numSamples;
%     else fprintf('Field numSamples has not been defined.'); return; end
% 
%     % Samples discarded or burn-in (default 0)
%     if isfield(alpha_options,'numDiscarded'); sample.numDiscarded = alpha_options.numDiscarded;
%     else sample.numDiscarded = 0; end
% 
%     % Steps per point (default 1e2)
%     if isfield(alpha_options,'stepsPerPoint'); sample.stepsPerPoint = alpha_options.stepsPerPoint;
%     else sample.stepsPerPoint = 1e2; end
% 
%     % Sampling algorithm 'HRB' [COMMENT] (¿USAMOS HRB?)
%     if isfield(alpha_options,'algorithm'); sample.algorithm = alpha_options.algorithm;
%     else sample.algorithm = 'ADSB'; end
% 
%     % Numerical flux tolerance  (default 1e-8)
%     if isfield(alpha_options,'vTol'); sample.vTol = alpha_options.vTol;
%     else sample.vTol = 1e-8; end
%     
%     % Parallel sampling option (default false or 0)
%     if isfield(alpha_options,'parallelFlag'); sample.parallelFlag = alpha_options.parallelFlag;
%     else sample.parallelFlag = 0; end
%     
%     % Number of cores for parallel sampling (default empty)
%     if isfield(alpha_options,'numCores'); sample.numCores = alpha_options.numCores;
%     else sample.numCores = []; end
%     
%     % Run MCMC diagnostics [COMMENT] NO SÉ SI SE PUEDEN REALIZAR ESTOS DIAGNÓSTICOS
%     if isfield(alpha_options,'diagnostics'); sample.diagnostics = alpha_options.diagnostics;
%     else sample.diagnostics = 1; end
%     
%     % Fold number of particles in the population relative to dim(Omega) (only for ADSB, default 1) 
%     if strcmp(sample.algorithm,'ADSB') && isfield(alpha_options,'populationScale'); sample.populationScale = alpha_options.populationScale;
%     else sample.populationScale = 3; end
%     
% 	% Restart sampler from previous results
% 	if ~isfield(model,'points'); sample.points = [];
%     else sample.points = model.points; end
% end
% 
% % Define solver parameters [COMMENT] ¿CORREPONDE ESTA PARTE?
% changeCobraSolverParams('MILP','intTol',1e-9);
% changeCobraSolverParams('MILP','relMipGapTol',1e-12);
% changeCobraSolverParams('LP','optTol',1e-9);
% changeCobraSolverParams('LP','feasTol',1e-9);
% 
% %% II. Preparation phase
% % Define fields related to the model in sample structure
% sample.rxns     = model.rxns;
% sample.numRxns  = numel(model.lb);
% sample.S        = model.S;
% sample.lb       = model.lb;
% sample.ub       = model.ub;
% sample.intRxns  = model.intRxns;
% sample.exchRxns = model.exchRxns;
% sample.uTol     = 1e-10;                                                   % Direction tolerance
% sample.bTol     = 1e-10;                                                   % Minimum allowed distance to closest constraint
% % [COMMENT] QUITÉ LA PARTE LOOPLESS
% 
% % Define fxn to ensure samples are kept within the bounds
% sample.keepWithinBounds = @(points) bringToBoundary(points,sample.lb,sample.ub);
% 
% % Generate initial warmup points
% sample.warmUpPoints = warmupLooplessACHRB(sample);
% 
% % Calculate centroid
% centroid = mean(sample.warmUpPoints,2);
% centroid(abs(centroid)<sample.vTol) = 0;
% 
% % % Check feasibility of the initial point. If centroid infeasible, find closest feasible point
% % if sample.loopless && ~sample.isFeasible(centroid) % [COMMENT] ¿PUEDO QUITAR LA CONDICIÓN DE QUE SEA LOOPLESS?
% %     sample.centroid = findNearestFeasiblePoint(sample,centroid,Nint',sample.vTol);
% % else
% %     sample.centroid = centroid;
% % end
% 
% sample.centroid = centroid;
% 
% fprintf('%d warmup points spanning the feasible space have been defined.\n',size(sample.warmUpPoints,2));
% clearvars -except sample
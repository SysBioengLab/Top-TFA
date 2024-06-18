function [points,samplingTime,rejectionRate,kernel] = metropolisHitAndRun(sample,numSamples,numDiscarted,stepsPerPoint,verbose,workerIdx)
% Hit and Run on a Box
%
% Uses HRB to generate a random sample from the loopless flux solution
% space
%
% USAGE:
%              points = HRB(sample,numSamples,numDiscarted,stepsPerPoint)
%
% INPUTS:
%              sample (structure):    (the following fields are required - others can be supplied)
%                                    * centroid - Initial centroid estimation
%                                    * warmupPoints:  Set of (loopless) flux solutions (2 x n)
%                                    * isFeasible - Fxn handle for looplessCheck.m fxn
%                                    * loopless - Loop removal option, true (default) or false
%              numSamples:           Number of points to sample
%              numDiscarded:         Burn-in (double)
%              stepsPerPoint:        Thining (double)
%
% OPTIONAL INPUTS:
%              verbose:    Display progress of sampling run, true (default) or false
%              workerIdx:  ID of the working processor, 1 (default)
%              cWeight:    Cumulative centroid weight from previous
%                          iterations
%
% OUTPUT:
%              points:   Matrix with n x numSamples (loopless) flux solutions
%
% OPTIONAL OUTPUT:
%              samplingTime:   Runtime of ll-ACHRB
%              centroid:       Estimated centroid
%              cWeight:        Cumulative centroid weight
%              rejectionRate:  Rejection rate estimated by the number of
%                              times the `hyperbox` is shrunk
%
% -------------------- Copyright (C) 2019 Pedro A. Saa --------------------
disp("Running MHR from clostridium")

% Check inputs
if nargin<4; disp('Not enough input parameters'); return; end

% Hit-And-Run parameters set-up
t0 = cputime;
if nargin<5
    verbose   = true;
    workerIdx = 1;
elseif nargin<6
    workerIdx = 1;
end

% Set-up HRB parameters
prevPoint = sample.centroid;

% Determine the dimension of the feasible space
Aeq = full(sample.S);
ix_fixed = (abs(sample.ub-sample.lb)<= 1e-6);
if any(ix_fixed)
    Aeq(:,ix_fixed) = 0;
end

% Allocating memory for the next points
nVars  = size(sample.S,2);
points = zeros(nVars,numSamples);
kernel = null(Aeq,'rational');
kernel(ix_fixed,:) = 0;
kernel(:,sum(abs(kernel),1)<=1e-6) = [];
nK     = size(kernel,2);             % Number of kernel vectors


% Define sampler parameters
uTol      = sample.uTol;               % Direction tolerance
maxMinTol = sample.bTol;               % Min distance to closest bound

% Initialize counters
counter       = 0;
sampleCount   = 0;
totalCount    = 0;
rejectedCount = 0;
rejectionRate = 0;

% Main loop
if verbose && (workerIdx==1); fprintf('%%Prog \t Time \t Time left \t Rejection rate\n--------------------------------------------\n'); end
while sampleCount<numSamples

    % Update count
    totalCount    = totalCount+1;
    rejectionRate = rejectedCount/totalCount;

    % Print step information
    if verbose && (workerIdx==1)
        if totalCount>1
            timeElapsed = (cputime-t0)/60;
            timePerStep = timeElapsed/sampleCount;
            if ~mod(totalCount,500*stepsPerPoint) && counter>0
                fprintf('%d\t%8.2f\t%8.2f\t%8.2f\n',round(1e2*sampleCount/numSamples),timeElapsed,(numSamples-sampleCount)*timePerStep,rejectionRate);
            elseif ~mod(totalCount,500*stepsPerPoint) && counter==0
                fprintf('%d\t%8.2f\t%8.2f\t%8.2f\n',round(1e2*sampleCount/numSamples),timeElapsed,(numSamples-sampleCount)*timePerStep,1);
            end
        end
    end

    % Return if the dynamics is frozen
    if rejectionRate>.9999 && totalCount>1e5
        disp('Hit-And-Run dynamics is frozen!');
        % save("debug.mat","rejectionRate","totalCount","prevPoint","nextPoint")
        points = []; break;
    end

    % Sample random direction from kernel
    u = kernel(:,randi([1 nK],1));

    % Figure out the distances to upper and lower bounds
    distUb = (sample.ub-prevPoint);
    distLb = (prevPoint-sample.lb);

    % Figure out positive and negative directions
    posDirn = (u>uTol);
    negDirn = (u<-uTol);

    % Figure out all the possible maximum and minimum step sizes
    maxStepTemp = distUb./u;
    minStepTemp = -distLb./u;
    maxStepVec  = [maxStepTemp(posDirn);minStepTemp(negDirn)];
    minStepVec  = [minStepTemp(posDirn);maxStepTemp(negDirn)];

    % Figure out the true max & min step sizes
    maxStep = min(maxStepVec);
    minStep = max(minStepVec);

    % Find new direction if we're getting too close to a constraint
    if (abs(minStep)<maxMinTol && abs(maxStep)<maxMinTol) || (minStep>maxStep) || (maxStep<0) || (minStep>0)
        rejectedCount = rejectedCount+1;
        continue;
    end

    % Pick a random step distance
    stepDist = rand(1)*(maxStep-minStep)+minStep;

    % Make move and check whether it is feasible
    nextPoint = prevPoint + stepDist*u;

    % Calculate probability of next point and ratio
    if sample.generalHR
        realRatio = 0; % ln(0) = 1;
    else
        realRatio = myMVNPDF(sample,prevPoint,nextPoint) + calculateExpNormalRatio(sample,prevPoint,nextPoint);
    end

    % Save current point if ratio>1 or greater than test
    if realRatio >= 0 || realRatio >= log(unifrnd(0,1))
        counter   = counter+1;
        prevPoint = nextPoint;
        
        if counter>numDiscarted && ~mod(counter,stepsPerPoint)
            sampleCount = sampleCount+1;
            points(:,sampleCount) = nextPoint;
        end
    else
        % Update rejected count
        rejectedCount = rejectedCount+1;
    end
end
samplingTime = (cputime-t0)/60;
if verbose && (workerIdx==1); fprintf('--------------------------------------------\n'); end

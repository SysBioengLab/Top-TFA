function findTopologiesNewClusters3(model,clusterID)
% Workflow to enumerate all possible directionalities

%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2015 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edited by Sebastián Zapararte
clc
folderName = ['clusterFiles\',clusterID];
try
    rmdir(folderName,'s');
    mkdir(folderName);
catch
    mkdir(folderName);
end
%% Step 0. Extract clusters
% model = ecModel1;

% Check if the matrix is OK (i.e. stoic. coeffs. similar)
% Scluster = Clusters_all{clusterID,1};
% prerevRxns  = Clusters_all{clusterID,2};
% external =  Clusters_all{clusterID,5};

% Balance the stoichiometric matrix if necesary
% while true
%     [maxCoeff,ix] = max(abs(Scluster(:)));
%     % check arbitrary cut-off
%     if maxCoeff>10        
%         % stoic. coeff. replaced by 1
%         Scluster(ix) = sign(Scluster(ix)); 
%     else
%         break;
%     end
% end

% Build model structure
% model.S    = sparse(Scluster);
% model.c    = ones(size(model.S,2),1);
% model.b    = zeros(size(model.S,1),1);
% model.lb   = -1e2*(prerevRxns==1);
% model.ub   = 1e2*ones(size(model.S,2),1);
% model.rxns = Clusters_all{clusterID,3};
% model.mets = Clusters_all{clusterID,4};
% model.rxnNames = model.rxns;
% model.metNames = model.mets;
% model.description = ['Cluster_',num2str(clusterID)];
% premodel=model;

%pre-processing of model, model compression
% [model,revRxns]=RemoveDuplicateRxns(premodel,prerevRxns,external);

%% Step 1. Build loopless structure
model = looplessStructureMILP(model,0);
disp(['Loopless structure succesfully loaded: rxns = ',num2str(length(model.rxns)),', mets = ',num2str(length(model.mets))]);

%% Step 2. Run ll-FVA
model = looplessFVA(model);
disp(['ll-FVA succesfully run: rxns = ',num2str(length(model.rxns)),', mets = ',num2str(length(model.mets))]);

%% Step 3. Find plausible directionalities
% Build directionality map
[model,directionSpace] = buildDirectionsRev(model);

%% 3.A Brut force method (this can take very long!) (comment otherwise)
tic
depthSearchDirectionsTestRev(model,directionSpace,['clusterFiles\tfm_cluster_full_exploration',clusterID])
toc

%% 3.B Topological tree search strategy
disp('############ Topological tree search ###########################')
tic
for constIdx  = 1:size(model.S,1)
    depthSearchUnfeasibleBalanceRev(model,directionSpace,constIdx,full(model.S(constIdx,:)),folderName);
    depthSearchUnfeasibleBalanceRev2(model,directionSpace,constIdx,folderName);
%     depthSearchUnfeasibleBalanceRev3(model,directionSpace,constIdx,folderName);
end

% Find unfeasible loops
for constIdx  = 1:size(model.Nint,1)
    depthSearchUnfeasibleLoopsRev(model,directionSpace,constIdx,full(model.Nint(constIdx,:)),folderName);
end

% Compile information
reorderPatterns = 0;    % (optional, default = 0) can be changed to 1 (smaller first), order of using patterns 
[model,unfeasiblePattern,directionSpace] = compilePatternsRev(model,directionSpace,folderName,reorderPatterns);
unfeasiblePattern = unique(unfeasiblePattern, 'rows');
% Perform global search
depthSearchDirectionsRev(model,unfeasiblePattern,directionSpace,['clusterFiles\tfm_cluster',clusterID]);
toc
tic
if ~isempty(unfeasiblePattern)
    FullunfeasiblePattern = unfeasiblePattern(:, logical(model.rev));
    FullunfeasiblePattern = unique(FullunfeasiblePattern, 'rows');
    % remove unusefull patterns, prior of reorganizing
    unfeasiblePattern=RemoveNoUsefullPatterns(FullunfeasiblePattern); 
    %Convert the patterns to 16 bit integers to be used in the c program
    unfeasiblePattern = int16(unfeasiblePattern);
    %Correct the format of the direction space
    %Find the position of the last non-zero element in each pattern
    lastNonZero = [];

    for i = 1:size(unfeasiblePattern,1)
        lastNonZero=[lastNonZero;find(unfeasiblePattern(i,:), 1, 'last')];
    end
    lastNonZero=int16(lastNonZero);
   [validConfigs, numCfgs] = juicy3_true_binary(unfeasiblePattern, lastNonZero, clusterID);
    %[validConfigs, numCfgs] = juicy3_true_binary_no_save(unfeasiblePattern, lastNonZero, clusterID);
end
toc

% Read and translate the final result from c program
% preFileName=num2str(clusterID);
% for i=1:3-length(preFileName)
%     preFileName=strcat('0',preFileName);
% end
% TotalConfig=0;
% Numconfig=zeros(692,1);
% % for i=1:692
% i=1;
%     preFileName2=num2str(i);
%     for j=1:4-length(preFileName2)
%         preFileName2=strcat('0',preFileName2);
%     end
% binaryFile=strcat('TFM_',preFileName,'_',preFileName2,'.bin');
%      BinaryOutput = readBinaryCfgs(binaryFile, sum(revRxns));
% %     
% %     Numconfig(i)=size(BinaryOutput,1);
% %     TotalConfig=TotalConfig+Numconfig(i);
% %     fprintf(1,'File #: %d Number of configs: %d.\n',i,TotalConfig);
% % end
% 
% Output_tfms = TranslateBinaryCfgs(BinaryOutput, sum(revRxns));
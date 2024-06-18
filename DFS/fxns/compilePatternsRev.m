function [model,unfeasiblePattern,directionSpace] = compilePatternsRev(model,directionSpace,folderName,optimizePatternOrder)
unfeasiblePattern = [];
files = dir(folderName);
% Compile unfeasible flux patterns
for j = 3:length(files)
    try load([folderName,'/',files(j).name]);
        disp([folderName,'/',files(j).name,' loaded']);
        unfeasiblePattern = [unfeasiblePattern;unfeasiblePatterns];        
    catch 
        break;
    end
end

% Remove extra entries of the model
model.lb(model.numRxns+1:end) = [];
model.ub(model.numRxns+1:end) = [];

% Order unfeasible list from the smallest to the tallest patterns
if optimizePatternOrder && ~isempty(unfeasiblePattern)  
    tempList = unfeasiblePattern;
    
    % Find active constraints and inactivate wild-cards (=0)
    tempList(tempList~=0)=1;
    tempList(tempList==0)=0;
    
    % Find optimal order
    [~,ix] = sort(sum(tempList),'ascend');
    
    % Reset parameters and model fields
    unfeasiblePattern = unfeasiblePattern(:,ix);
    directionSpace    = directionSpace(ix);
    model.lb          = model.lb(ix);
    model.ub          = model.ub(ix);
    model.rev         = model.rev(ix);
    model.rxns        = model.rxns(ix);
    model.rxnNames    = model.rxnNames(ix);
    model.S           = model.S(:,ix);
end

% Remove unnecessary fields
fields = {'internal','exchange','Nint','obj','modelsense','A','rhs','sense',...
        'vtype','metsIdxs','loopsIdxs'};
model  = rmfield(model,fields);
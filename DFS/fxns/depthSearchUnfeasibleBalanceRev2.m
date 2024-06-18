function depthSearchUnfeasibleBalanceRev2(model,directionSpace,constIdx,folderName)
% Finds unfeasible flux patterns based on mass balance constraints with two
% metabolits conected
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2015 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of global variables
unfeasiblePatterns    = [];
counter = 0;

Stoic=full(model.S(constIdx,:)); %%%added Vero
% Check if we have both active and rev reactions for this constraint
[~,preactiveFluxes] = find(Stoic~=0);
Stoic2=full(model.S(:,preactiveFluxes)); %%%added Vero
Stoic2(constIdx,:)=zeros(1,size(Stoic2,2)); %%%added Vero, remove from the matrix the studied metabolite
Stoic2(logical(Stoic2))=1;    %%%added Vero
check=find(sum(Stoic2,2)>1);    %%%added Vero, find metabolites that share more than one reaction with the studied metabolite
reverRxns      = find(model.rev==1);

for i=1:length(check)       %%%added Vero
 if check(i)>constIdx       %%%added Vero to do not repeat the previous rules again
     Stoic=sum(full(model.S([constIdx,check(i)],:)));   %%%added Vero external reactions to that two metabolites
     [~,activeFluxes] = find(Stoic~=0);             %%%added Vero
     compoundVector   = [activeFluxes(:);reverRxns(:)]; 
     if length(unique(compoundVector))==length(compoundVector)
         return;
     else
    % Determine the reversible rxns belonging to the active set
         activeRxns = [];
         dirTemp  = zeros(1,model.numRxns);
         for ix = activeFluxes
             if any(ix==reverRxns)
                 activeRxns = [activeRxns,ix];
             else
                 dirTemp(ix) = directionSpace{ix};
             end
         end
     end
     
     % Search parameters
     depth = length(activeRxns);
     
     % Perform depth-first search: recursive call
     depthSearchUnbalance(directionSpace,dirTemp,0);
     
     if ~isempty(unfeasiblePatterns)
         disp(['Unfeasible modes found in double balance ',num2str(constIdx),'  ',num2str(check(i)),'  = ',num2str(size(unfeasiblePatterns,1))]);
         save([folderName,'/unbalancedTopology2_',num2str(constIdx),'_',num2str(check(i)),'.mat']);
     end
 end
end
 
function depthSearchUnbalance(directionSpace,dirTemp,level)
        if depth == level
            counter = counter+1;
            
            % Check mass balances topologicaly
            dirProd = Stoic(activeFluxes).*dirTemp(activeFluxes);
            if all(dirProd>0) || all(dirProd<0)
                tempPattern = zeros(1,model.numRxns);
                tempPattern(activeRxns) = dirTemp(activeRxns);
                unfeasiblePatterns = [unfeasiblePatterns;tempPattern];
            end
            
            % Recursive call
        else
            
            level = level+1;
            for j = 1:length(directionSpace{activeRxns(level)})
                dirTemp(activeRxns(level)) = directionSpace{activeRxns(level)}(j);
                depthSearchUnbalance(directionSpace,dirTemp,level);
            end
        end
    end
end


  
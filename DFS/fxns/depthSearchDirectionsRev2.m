function depthSearchDirectionsRev2(model,unfeasiblePattern,directionSpace,filename)
% Performs a depth-first search through the feasible space of directions
% using a table of unfeasible patterns
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2015 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified by Veronica Martinez
% Definition of global variables
counter  = 0;
progress = 0;
tfm      = zeros(floor(3E9/(8*length(find(model.rev)))),length(find(model.rev))); 
counter2 = 0;

% Search parameters
reverRxns = find(model.rev==1);
dirTemp = zeros(1,model.numRxns);
depth   = length(reverRxns);
disp(['Tree depth: ',num2str(depth)]);

% Define directionalities for each
for ix = 1:model.numRxns
    % reversible reactions can be either + or -
    if all(reverRxns~=ix)
        if model.lb(ix)<0
            dirTemp(ix) = -1;
        elseif model.ub(ix)>0
            dirTemp(ix) = 1;
        end
    end
end

% Perform depth-first search: recursive call
disp('Initializing tfm search...');
depthSearchRev(directionSpace,dirTemp,0);
disp(['Topological flux modes found = ',num2str(size(tfm,1))]);
clearvars -except tfm counter model unfeasiblePattern filename
filename2=strcat(filename,'.mat');
save(filename2,'-v7.3')

    function depthSearchRev(directionSpace,dirTemp,level)
        % Add solution to the feasible set
        if depth == level
            counter = counter+1;
            counter2=counter2+1;
            temp = dirTemp(find(model.rev));
            tfm(counter2,:) = temp;     
            if ~mod(counter,size(tfm,1))
                filename2=strcat(filename,'_',num2str(counter/counter2),'.mat');
                save(filename2,'-v7.3','tfm')
                tfm  = zeros(floor(3E9/(8*length(find(model.rev)))),length(find(model.rev)));
                counter2=0;
            end
            
            if ~mod(counter,1e4)
                progress = progress+counter;
                disp(['Current TFMs: ',num2str(counter),', Progress so far: ',num2str(fix(100*progress/2^depth)),'%']);
            end
            % Recursive call
        else
            level = level+1;
            for j = 1:length(directionSpace{reverRxns(level)})
                dirTemp(reverRxns(level)) = directionSpace{reverRxns(level)}(j);
                if isempty(unfeasiblePattern)
                    depthSearchRev(directionSpace,dirTemp,level);
                    % Check if it is worthwhile to further deepen
                else
                    checkPattern = dirTemp(ones(1,size(unfeasiblePattern,1)),:);
                    checkPattern(unfeasiblePattern==0) = 0;
                    if all(sum(checkPattern==unfeasiblePattern,2)~=model.numRxns)
                        depthSearchRev(directionSpace,dirTemp,level);
                    else
                        if depth ~= level
                            progress = progress+2^(depth-level);
                        end
                    end
                end
            end
        end
    end
end
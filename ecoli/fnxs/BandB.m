function [sol] = BandB(IM,RM,MM,MMf)

disp('Starting Branch and Bound')
t0 = cputime;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization
N          = size(RM,2);                              % Number of reversible rxns (variables)
% sol        = sparse(N,sumComb(N));                    % Pre-allocate memory (pretty big)
sol        = zeros(N,1);
options    = zeros(N,1);
x          = zeros(N,1);
options(1) = 1;                                       % Starting point
x(1)       = 1;                                       % Starting point
sol_count  = 1;
node       = 1;
lastnodes  = [];
subnodes   = [];
SNM        = zeros(N,N);                              % Sub-Node Matrix
cyclesMade = 0;

% Loop
% while sum(abs(options)) > 0
while 1
    cyclesMade = cyclesMade +1;
    if rem(cyclesMade,1000) == 0
        cyclesMade
        size(sol,2)
    end
%     options
%     x
    lastnodes(end+1) = node;
%     lastnodes
    % Acople
    mergeNodes = [node subnodes];
    mergeNodes = mergeNodes(ismember(mergeNodes,MM));
%     mergeNodes
    groupsUsed = [];
    for n=mergeNodes
        [r,c] = find(MM==n);
        if ~ismember(r,groupsUsed)
            groupsUsed(end+1) = r;
            group = MM(r,:);
            group(group==0) = [];
            mergeFactor = MMf(r,c);
            group(c) = [];
            fgroup = MMf(r,:);
            fgroup(fgroup==0) = [];
            fgroup(c) = [];
            if mergeFactor == x(n)
                x(group) =  fgroup;
            else
                x(group) = -fgroup;
            end
        end
    end
%     disp('x dsps del acople:')
%     x'
    subnodes = [];
    if sum(abs(x)) == N                               % Check if x is complete
%         disp('X completo')
        if FEASIBLE(IM,RM,x)
%             disp('X completo factible')
            sol(:,sol_count) = x;
            sol_count = sol_count +1;
        end
        if sum(options) == 0                          % All options have been checked
%             disp('Todas las opciones revisadas')
            break
        else
%             disp('Llego al final y ahora hace comeback')
            [node,lastnodes,options,x,SNM] = comeback(lastnodes,options,x,SNM,MM);
        end
    else
        if UNFEASIBLE(IM,RM,x)
%             disp('Es infactible:')
%             x
            if sum(options) == 0                          % All options have been checked
                break
            else
%                 disp('Como es infactible hace comeback')
                [node,lastnodes,options,x,SNM] = comeback(lastnodes,options,x,SNM,MM);
            end
        else
%             disp('No es infactible')
            [X,f] = case1(IM,RM,x);
            if isempty(X)
%                 disp('No hay caso 1')
                for i=node+1:N
                    if x(i) == 0
                        nextNode = i;
                        x(i) = 1;
                        options(i) = 1;
                        break
                    end
                end
                node = nextNode;
            else
%                 disp('Hay caso 1')
                z_count = 0;
                for i=1:N
                    if x(i) == 0
                        z_count = z_count + 1;
                        pos = find(X==z_count);
                        if ~isempty(pos)
                            x(i) = f(pos);
                            subnodes(end+1) = i;
                        end
                    end
                end
%                 subnodes
                SNM(node,:) = [subnodes zeros(1,N-size(subnodes,2))];
                lastnodes(end) = [];
                % Repeat loop with case1 subnodes in it
            end
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prepTime = (cputime-t0)/60;
disp(['Branch and Bound ended after ',num2str(prepTime),' seconds'])
disp([num2str(size(sol,2)),' feasible topologies found'])
end
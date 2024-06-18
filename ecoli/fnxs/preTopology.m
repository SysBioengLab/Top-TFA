function [IM,RM,MM,MMf] = preTopology(model)
% Returns the matrices needed to do the topological analysis of the network
lb = model.lb;
ub = model.ub;
S = model.S;
S = sign(S);              % Only flux direction matters
rev = zeros(size(lb));
rev(lb < 0 & ub > 0) = 1;
rev = logical(rev);
RM  = S(:,rev);           % Reversible matrix
IM  = S(:,~rev);          % Irreversible matrix

Nm   = size(RM,1);        % Number of mets
Nr   = size(RM,2);        % Number of rev rxns
MM  = zeros(Nm,Nr);       % Merge matrix
MMf = zeros(Nm,Nr);       % Factors of merge matrix

merge_candidates = zeros(Nm,1);
merge_candidates(sum(abs(IM),2) == 0 & sum(abs(RM),2) == 2) = 1;
mc_positions = find(merge_candidates);
merge_vector = zeros(sum(merge_candidates),2);
for i=1:sum(merge_candidates)
    merge_vector(i,:) = find(RM(mc_positions(i),:));
end

c = 0;
level = 1;

while ~isempty(merge_vector)
    if c == 0
        m1 = merge_vector(1,:);
    end
    n = size(merge_vector,1);
    c = 0;
    for i=2:n
        m2 = merge_vector(i-c,:);
        if any(ismember(m1,m2))
            m1 = union(m1,m2);
            merge_vector(i-c,:) = [];
            c = c+1;
        end
    end
    if c == 0
        MM(level,:) = [m1 zeros(1,Nr-size(m1,2))];
        level = level + 1;
        merge_vector(1,:) = [];
    end
end
% MMf

for i=1:Nm
    if MM(i,:) == 0
        break
    end
    line = MM(i,:);
    line(line==0) = [];
    prior  = line(1);
    fs = zeros(size(line));
    fs(1) = 1;
    checked = [];
    while sum(logical(fs)) < size(line,2)
        check = find(fs==0);
        for j=check
            rxnToCheck = line(j);
            inter = intersect(find(RM(:,prior)),find(RM(:,rxnToCheck)));
            if ~isempty(inter)
                for k2=1:size(inter,1)
                    k = inter(k2);
                    if (isempty(IM) || all(IM(k,:) == 0)) && sum(abs(RM(k,:)))==2
                        if sign(RM(k,prior)) ~= sign(RM(k,rxnToCheck))
                            fs(j) =  1;
                        else
                            fs(j) = -1;
                        end
                        break
                    end
                end
            end
        end
        if sum(logical(fs)) < size(line,2)
            checked(end+1) = prior;
            nextOptions = line(logical(fs));
            for p=1:size(nextOptions,2)
                if ~ismember(nextOptions(p),checked)
                    prior = nextOptions(p);
                    break
                end
            end
            if prior == checked(end)
                prior = line(find(fs==0,1));
            end
        end
    end
    MMf(i,:) = [fs zeros(1,Nr-size(fs,2))];
end

end
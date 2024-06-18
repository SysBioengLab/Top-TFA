function F = FEASIBLE(IM,RM,x)
% Checks if vector x is feasible
nRM = (x .* RM.').';
M = [IM nRM];
F = 1;
for i=1:size(M,1)
    if all(M(i,:) <= 0) || all(M(i,:) >= 0)
        F = 0;
        break
    end
end
end

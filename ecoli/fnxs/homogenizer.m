function [model] = homogenizer(model)
% Transforms non-homogeneus problem: Ax = b into homogeneus problem Ax=0.

nB       = sum(~(model.b==0));
rhMatrix = -eye(size(model.b,1));
rhMatrix(:,model.b==0) = [];
model.S  = [model.S rhMatrix];
bBound   = model.b;
bBound(model.b==0) = [];
model.lb = [model.lb; bBound];
model.ub = [model.ub; bBound];
model.b  = zeros(size(model.b));

% Add new "variables" names
bname1 = cell(nB,1);
bname1(:) = {'rhs_'};
bname2 = cellstr(string(1:nB)');
bname3 = strcat(bname1,bname2);
model.vars = [model.vars; bname3];

end
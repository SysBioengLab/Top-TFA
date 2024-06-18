function [cons, b] = closa_addReductionConstraints(tModel)
% Returns reduction constraints
tempModel = tModel;

% Find metabolites
r1 = 0; r2 = 0;

if any(ismember(tempModel.metNames,'NAD'))
    r1 = 1;
    c1 = double(ismember(tempModel.metNames,'NADH') - ismember(tempModel.metNames,'NAD'))';
    c2 = -c1;
    b12 = [log(10); 4*log(10)];
end

if any(ismember(tempModel.metNames,'NADP'))
    r2 = 1;
    c3 = double(ismember(tempModel.metNames,'NADPH') - ismember(tempModel.metNames,'NADP'))';
    c4 = -c3;
    b34 = [log(10); 4*log(10)];
end

if r1==1 && r2==1
    cons = [c1;c2;c3;c4];
    b = [b12;b34];
elseif r1==1 && r2==0
    cons = [c1;c2];
    b = b12;
elseif r1==0 && r2==1
    cons = [c3;c4];
    b = b34;
else
    cons=[];
    b=[];
end

end
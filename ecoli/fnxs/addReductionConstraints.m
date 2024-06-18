function [cons, b] = addReductionConstraints(tModel)
% Returns reduction constraints
model = tModel;

% Find metabolites
r1 = 0; r2 = 0;

if any(ismember(model.mets,'nad_c'))
    r1 = 1;
    c1 = double(ismember(model.mets,'nadh_c') - ismember(model.mets,'nad_c'))';
    c2 = -c1;
    b12 = [log(10); 4*log(10)];
end

if any(ismember(model.mets,'nadp_c'))
    r2 = 1;
    c3 = double(ismember(model.mets,'nadph_c') - ismember(model.mets,'nadp_c'))';
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
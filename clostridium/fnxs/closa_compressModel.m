function model = closa_compressModel(model)
% Removes blocked reactions, metabolites and fix the sense of irreversible
% backward reactions to the positive direction.
% Fields S, lb and ub are required.

% Remove blocked rxns (by bounds)
bR = zeros(size(model.S,2),1);
bR(model.lb==0 & model.ub==0) = 1;
bR = logical(bR);
if sum(bR)>0
    disp(["compressModel (",model.description,"): Removing the following reactions:"])
    disp(model.rxns(bR))
end
model.grRules(bR)      = [];
model.rules(bR)        = [];
model.rxns(bR)         = [];
model.rxnNames(bR)     = [];
model.rxnRheaID(bR)    = [];
model.rxnisabbreviated_nameID(bR) = [];
model.subSystems(bR)   = [];
model.S(:,bR)          = [];
model.T(:,bR)          = [];
model.charVec(bR)      = [];
model.lb(bR)           = [];
model.ub(bR)           = [];
model.c(bR)            = [];
model.rev(bR)          = [];
model.vmax(bR)         = [];
model.vmin(bR)         = [];

% Check if there are blockedMets
ibM = sum(abs(sign(model.S)),2) == 0;
% Raise warning if there are unbalanced mets
ubM = sum(abs(sign(model.S)),2) == 1;
if sum(ubM) > 0
    disp('Warning:')
    disp(['In the model ',model.description,...
        ' the following metabolites are unbalanced:'])
    disp(model.mets(ubM)')
end

if sum(bR) > 0 || sum(ibM) > 0
    nbR = 1;
    while sum(nbR) > 0
        % Remove blocked rxns and metabolites
        bM = zeros(size(model.S,1),1);
        notZ = model.S ~= 0;
        bM(sum(notZ,2)==0) = 1;
        bM = logical(bM);
        if sum(bM)>0
            disp(["compressModel (",model.description,"): Removing the following metabolites:"])
            disp(model.metNames(bM))
        end
        model.S(bM,:)    = [];
        model.T(bM,:)    = [];
        model.csense(bM) = [];
        model.mets(bM)        = [];
        model.metNames(bM)    = [];
        model.metFormulas(bM) = [];
        model.metCharges(bM)  = [];
        model.metChEBIID(bM)  = [];
        model.metisabbreviated_nameID = [];
        model.b(bM)          = [];
        model.lbc(bM)        = [];
        model.ubc(bM)        = [];
        model.dGft_mu(bM)    = [];
        model.dGft_std(bM)   = [];
        model.dGft_cov(bM,:) = [];
        model.dGft_cov(:,bM) = [];
        model.extMets(bM)    = [];
        model.NH(bM)         = [];
        model.conc_mu(bM)    = [];
        model.conc_std(bM)   = [];

        if isfield(model,"unbalancedMets")
            model.unbalancedMets(bM) = [];
        end
%         model.comp(bM)= [];           % 10.03: Borrado (temporalmente?)
        nbR = zeros(size(model.S,2),1);
        nbR(sum(abs(model.S),1)==0) = 1;
        nbR = logical(nbR);
        if sum(nbR)>0
            disp(["compressModel (",model.description,"): Removing the following reactions:"])
            disp(model.rxnNames(nbR))
        end
        model.grRules(nbR)      = [];
        model.rules(nbR)        = [];
        model.rxns(nbR)         = [];
        model.rxnNames(nbR)     = [];
        model.rxnRheaID(nbR)    = [];
        model.rxnisabbreviated_nameID(nbR) = [];
        model.subSystems(nbR)   = [];
        model.S(:,nbR)          = [];
        model.T(:,nbR)          = [];
        model.charVec(nbR)      = [];
        model.lb(nbR)           = [];
        model.ub(nbR)           = [];
        model.c(nbR)            = [];
        model.rev(nbR)          = [];
        model.vmax(nbR)         = [];
        model.vmin(nbR)         = [];
    end
end

end
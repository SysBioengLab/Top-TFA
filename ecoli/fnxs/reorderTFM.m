function new_tfm = reorderTFM(tfm,model,ecModel)
% Reorders tfm values based on the rxns IDs.

[~,idx] = ismember(ecModel.rxns(logical(ecModel.rev)),model.rxns(logical(model.rev)));
new_tfm = tfm(:,idx);

end
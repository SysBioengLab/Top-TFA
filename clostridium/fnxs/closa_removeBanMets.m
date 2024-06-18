function newModel = closa_removeProtons(oldModel)
% Replaces with zero the rows with banMets. Intended for tModel, NOT
% massModel. It also makes zero any rxn that involves a banMet.
% It compresses the tModel.

newModel = oldModel;
banRxns = any(newModel.S(oldModel.banMets,:));
newModel.S(:,banRxns) = zeros(size(newModel.S,1),sum(banRxns));
newModel.S(oldModel.banMets,:) = zeros(sum(oldModel.banMets),size(newModel.rxns,1));

newModel = closa_compressModel(newModel);

end
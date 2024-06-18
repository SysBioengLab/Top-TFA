function model = removeProtons(tModel)

model = tModel;
protons = logical(ismember(model.mets,'h_c') + ismember(model.mets,'h_e'));
model.S(protons,:) = model.T(protons,:);

end
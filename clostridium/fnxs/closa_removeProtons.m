function newModel = closa_removeProtons(oldModel)
% Intended for tModel, NOT massModel.

newModel = oldModel;
protons = ismember(newModel.metNames,'H+');
newModel.S(protons,:) = newModel.T(protons,:); % Removes H+ except transported H+

end
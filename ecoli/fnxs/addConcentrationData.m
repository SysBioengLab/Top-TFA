function model = addConcentrationData(model)
% Adds concentration data obtained from Bennett 2009 [Manually].

model.lbc(2)  = log(3.5e-3); model.ubc(2) = log(4.15e-3); % L-glutamine (citosol)
model.lbc(4)  = log(9.24e-2); model.ubc(4) = log(9.98e-2); % L-Glutamate (citosol)
model.lbc(14) = log(1.66e-3); model.ubc(14) = log(1.7e-3); % L-Malate
model.lbc(16) = log(2.32e-3); model.ubc(16) = log(2.8e-3); % NAD+
model.lbc(17) = log(5.45e-5); model.ubc(17) = log(1.27e-4); % NADH
model.lbc(18) = log(1.4e-7); model.ubc(18) = log(3.11e-5); % NADP+
model.lbc(19) = log(1.1e-4); model.ubc(19) = log(1.34e-4); % NADPH
model.lbc(28) = log(1.46e-4); model.ubc(28) = log(2.31e-4); % Phosphoenolpyruvate
model.lbc(29) = log(3.69e-3); model.ubc(29) = log(3.85e-3); % 6-Phospho-D-gluconate
model.lbc(38) = log(1e-4); model.ubc(38) = log(1.77e-3); % Alpha-D-Ribose 5-phosphate (lb divided by ~9)
model.lbc(39) = log(1e-4); model.ubc(39) = log(1.77e-3); % D-Ribulose 5-phosphate (lb divided by ~9)
model.lbc(44) = log(5.29e-4); model.ubc(44) = log(6.94e-4); % Acetyl-CoA
model.lbc(45) = log(3.41e-4); model.ubc(45) = log(9.49e-4); % Succinate
model.lbc(47) = log(1.42e-4); model.ubc(47) = log(3.83e-4); % Succinyl-CoA
model.lbc(48) = log(1.38e-5); model.ubc(48) = log(1.88e-5); % Cis-Aconitate
model.lbc(49) = log(1e-4); model.ubc(49) = log(1.77e-3); % D-Xylulose 5-phosphate (lb divided by ~9)
model.lbc(50) = log(1.02e-3); model.ubc(50) = log(1.13e-3); % Acetyl phosphate
model.lbc(51) = log(4.37e-4); model.ubc(51) = log(7.04e-4); % ADP
model.lbc(54) = log(2.32e-4); model.ubc(54) = log(3.41e-4); % AMP
model.lbc(55) = log(8.13e-3); model.ubc(55) = log(1.14e-2); % ATP
model.lbc(56) = log(1.1e-3); model.ubc(56) = log(3.48e-3); % Citrate
model.lbc(59) = log(8.83e-5); model.ubc(59) = log(2.12e-2); % Coenzyme A
model.lbc(60) = log(3.44e-4); model.ubc(60) = log(4.05e-4); % Dihydroxyacetone phosphate
model.lbc(64) = log(1e-3); model.ubc(64) = log(9.08e-3); % D-Fructose 6-phosphate (lb reduced)
model.lbc(65) = log(1.4e-2); model.ubc(65) = log(1.64e-2); % D-Fructose 1,6-bisphosphate
model.lbc(69) = log(3e-6); model.ubc(69) = log(4.42e-3); % Fumarate
model.lbc(72) = log(1e-3); model.ubc(72) = log(9.08e-3); % D-Glucose 6-phosphate (lb reduced)

% Mean and std of mets (In order)
meditionIndexes = [2 4 14 16 17 18 19 28 29 38 39 44 45 47 48 49 50 51 54 55 56 59 60 64 65 69 72];
model.concMu     = (model.lbc(meditionIndexes)+model.ubc(meditionIndexes))/2;
model.concStd    = (model.lbc(meditionIndexes)+model.ubc(meditionIndexes))/4;
model.concNormalMets = model.mets(meditionIndexes); % Leave empty if no meditions available
end
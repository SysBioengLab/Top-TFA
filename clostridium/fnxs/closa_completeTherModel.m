function model = closa_completeTherModel(model)
% Adds transport matrix and charge vector.
% This is manual and hardcoded

model.extMets = contains(model.mets,'[e]');
model.T       = zeros(size(model.S));
model.charVec = zeros(size(model.rxns));

% There are 8 inter-compartment reactions
% Indexes of reactions and protons
ihc      = find(ismember(model.mets,'c242[c]'));
ihe      = find(ismember(model.mets,'c242[e]'));
iCO2t    = find(ismember(model.rxns,'r1567'));                             % CO2 transport via diffusion
iH2Ot    = find(ismember(model.rxns,'r2030'));                             % H2O transport
iETHt    = find(ismember(model.rxns,'r2499'));                             % Ethanol transport
iH2t     = find(ismember(model.rxns,'r2522'));                             % Hydrogen transport
iACt     = find(ismember(model.rxns,'r2574'));                             % Acetate transport
iBDOt    = find(ismember(model.rxns,'r3066'));                             % Butane-2,3-diol transport
iRNFc    = find(ismember(model.rxns,'r3107'));                             % Rnf complex
iFIFOATP = find(ismember(model.rxns,'r3130'));                             % F1FO ATP synthase (3.66 protons for one ATP)

%% T and charVec construction (se lee de izq a der)
% CO2t: Simple difussion, no charge
model.T(find(ismember(model.mets,'c1466[c]')),iCO2t) =  1;
model.T(find(ismember(model.mets,'c1466[e]')),iCO2t) = -1;
% H2Ot: Simple difussion, no charge
model.T(find(ismember(model.mets,'c1392[c]')),iH2Ot) =  1;
model.T(find(ismember(model.mets,'c1392[e]')),iH2Ot) = -1;
% ETHt: Simple difussion, no charge
model.T(find(ismember(model.mets,'c1298[e]')),iETHt) =  1;
model.T(find(ismember(model.mets,'c1298[c]')),iETHt) = -1;
% H2t: Simple difussion, no charge
model.T(find(ismember(model.mets,'c1311[c]')),iH2t) =  1;
model.T(find(ismember(model.mets,'c1311[e]')),iH2t) = -1;
% ACt: 1 proton 1 acetate, net charge = 0
model.T(ihe,iACt) =  1;
model.T(ihc,iACt) = -1;
model.T(find(ismember(model.mets,'c241[e]')),iACt) =  1;
model.T(find(ismember(model.mets,'c241[c]')),iACt) = -1;
% BDOt: Simple difussion, no charge
model.T(find(ismember(model.mets,'c1305[e]')),iBDOt) =  1;
model.T(find(ismember(model.mets,'c1305[c]')),iBDOt) = -1;
% RNF complex: 2 protons transported
model.T(ihc,iRNFc) = -2;
model.T(ihe,iRNFc) =  2;
model.charVec(iRNFc) = -2;
% iFIFOATP: 3.66 protons transported (probably 11 for 3 ATP)
model.T(ihc,iFIFOATP) =  3.66;
model.T(ihe,iFIFOATP) = -3.66;
model.charVec(iFIFOATP) = 3.66;

end
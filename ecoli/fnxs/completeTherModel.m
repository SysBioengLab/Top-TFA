function model = completeTherModel(model)
% Adds transport matrix and charge vector.
% This is manual and hardcoded

model.extMets = contains(model.mets,'_e');
model.T       = zeros(size(model.S));
model.charVec = zeros(size(model.rxns));
% V 1.0: Only the following reactions are going to be adressed:
% PYRt2, CO2t, GLCpts, H2Ot, O2t, PIt2r, ATPS4r, CYTBD, THD2, NADH16
% Equations and definitions are based on e_coli_core paper (Orth, 2009).

% Indexes of reactions and protons
ihc     = find(ismember(model.mets,'h_c'));
ihe     = find(ismember(model.mets,'h_e'));
iPYR    = find(ismember(model.rxns,'PYRt2'));
iCO2t   = find(ismember(model.rxns,'CO2t'));
iGLCpts = find(ismember(model.rxns,'GLCpts'));
iH2Ot   = find(ismember(model.rxns,'H2Ot'));
iO2t    = find(ismember(model.rxns,'O2t'));
iPIt2r  = find(ismember(model.rxns,'PIt2r'));
iATPS4r = find(ismember(model.rxns,'ATPS4r'));
iCYTBD  = find(ismember(model.rxns,'CYTBD'));
iTHD2   = find(ismember(model.rxns,'THD2'));
iNADH16 = find(ismember(model.rxns,'NADH16'));

%% T and charVec construction
% NADH16: 3 protons are transported
model.T(ihc,iNADH16) = -3;
model.T(ihe,iNADH16) = 3;
model.charVec(iNADH16) = -3;
% CYTBD: 2 protons are transported
model.T(ihc,iCYTBD) = -2;
model.T(ihe,iCYTBD) = 2;
model.charVec(iCYTBD) = -2;
% O2t: 1 molecule of o2 is transported [REVERSIBLE]
model.T(find(ismember(model.mets,'o2_c')),iO2t) = -1;
model.T(find(ismember(model.mets,'o2_e')),iO2t) = 1;
% ATPS4r: 4 protons are transported [REVERSIBLE]
model.T(ihc,iATPS4r) =  4;
model.T(ihe,iATPS4r) = -4;
model.charVec(iATPS4r) = 4;
% THD2: 2 protons
model.T(ihc,iTHD2) =  2;
model.T(ihe,iTHD2) = -2;
model.charVec(iTHD2) = 2;
% PIt2r: 1 proton
model.T(ihc,iPIt2r) =  1;
model.T(ihe,iPIt2r) = -1;
model.T(find(ismember(model.mets,'pi_c')),iPIt2r) =  1;
model.T(find(ismember(model.mets,'pi_e')),iPIt2r) = -1;
model.charVec(iPIt2r) = -1;
% GLCpts: glucose to g6p
model.T(find(ismember(model.mets,'glc__D_e')),iGLCpts) = -1;
model.T(find(ismember(model.mets,'g6p_c')),iGLCpts) =  1;
% PYRt2: 1 proton and 1 pyruvate
model.T(ihc,iPYR) =  1;
model.T(ihe,iPYR) = -1;
model.T(find(ismember(model.mets,'pyr_c')),iPYR) =  1;
model.T(find(ismember(model.mets,'pyr_e')),iPYR) = -1;
% Sum of charges = 0
% CO2t: No charge. Simple diffusion
model.T(find(ismember(model.mets,'co2_c')),iCO2t) =  1;
model.T(find(ismember(model.mets,'co2_e')),iCO2t) = -1;
% H2Ot: No charge. Simple diffusion
model.T(find(ismember(model.mets,'h2o_c')),iH2Ot) =  1;
model.T(find(ismember(model.mets,'h2o_e')),iH2Ot) = -1;

end
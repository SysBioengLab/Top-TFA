function charge = chargeMatrix(model)
% Build charge matrix of the model

charge = zeros(size(model.S));

% Find indexes of rxns and mets
ACALDt = contains(model.rxns,'ACALDt');
acald_c = contains(model.mets,'acald_c');
ACt2r = contains(model.rxns,'ACt2r');
ac_c = contains(model.mets,'ac_c');
h_c = contains(model.mets,'h_c');
AKGt2r = contains(model.rxns,'AKGt2r');
akg_c = contains(model.mets,'akg_c');
ATPS4r_r = contains(model.rxns,'ATPS4r_r');
CO2t_r = contains(model.rxns,'CO2t_r');
co2_c = contains(model.mets,'co2_c');
CYTBD = contains(model.rxns,'CYTBD');
D_LACt2 = contains(model.rxns,'D_LACt2');
lac__D_c = contains(model.mets,'lac__D_c');
ETOHt2r = contains(model.rxns,'ETOHt2r');
etoh_c = contains(model.mets,'etoh_c');
FORt = contains(model.rxns,'FORt');
for_c = contains(model.mets,'for_c');
GLUt2r = contains(model.rxns,'GLUt2r');
glu__L_c = contains(model.mets,'glu__L_c');
H2Ot_r = contains(model.rxns,'H2Ot_r');
h2o_c = contains(model.mets,'h2o_c');
NADH16 = contains(model.rxns,'NADH16');
PYRt2 = contains(model.rxns,'PYRt2');
pyr_c = contains(model.mets,'pyr_c');
SUCCt3 = contains(model.rxns,'PYRt2');
succ_c = contains(model.mets,'succ_c');

% Manually add corresponding charges
charge(acald_c,ACALDt) = model.charge(acald_c);
charge(ac_c,ACt2r) = model.charge(ac_c);
charge(h_c,ACt2r)  = model.charge(h_c);
charge(akg_c,AKGt2r) = model.charge(akg_c);
charge(h_c,AKGt2r) = model.charge(h_c);
charge(h_c,ATPS4r_r) = model.charge(h_c);
charge(co2_c,CO2t_r) = model.charge(co2_c);
charge(h_c,CYTBD) = model.charge(h_c);
charge(h_c,D_LACt2) = model.charge(h_c);
charge(lac__D_c,D_LACt2) = model.charge(lac__D_c);
charge(etoh_c,ETOHt2r) = model.charge(etoh_c);
charge(h_c,ETOHt2r) = model.charge(h_c);
charge(for_c,FORt) = model.charge(for_c);
charge(glu__L_c,GLUt2r) = model.charge(glu__L_c);
charge(h_c,GLUt2r) = model.charge(h_c);
charge(h2o_c,H2Ot_r) = model.charge(h2o_c);
charge(h_c,NADH16) = model.charge(h_c);
charge(pyr_c,PYRt2) = model.charge(pyr_c);
charge(h_c,PYRt2) = model.charge(h_c);
charge(succ_c,SUCCt3) = model.charge(succ_c);

end
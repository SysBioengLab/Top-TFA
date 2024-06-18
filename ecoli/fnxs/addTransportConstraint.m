function [charge,h_trans]= addTransportConstraint(tModel)
% Search manually selected rxns. Returns charge and h_trans vector.

charge  = zeros(size(tModel.rxns));
h_trans = charge;
if any(ismember(tModel.rxns,'PYRt2'))
    nPYR = ismember(tModel.rxns,'PYRt2');
    charge(nPYR)  = tModel.metCharge(ismember(tModel.mets,'pyr_e'));
    h_trans(nPYR) = 1;
end
if any(ismember(tModel.rxns,'ATPS4r'))
    nATP = ismember(tModel.rxns,'ATPS4r');
    charge(nATP)  = 0;
    h_trans(nATP) = 4;
end
if any(ismember(tModel.rxns,'CYTBD'))
    nCYT = ismember(tModel.rxns,'CYTBD');
    charge(nCYT)  = 0;
    h_trans(nCYT) = -2;
end
if any(ismember(tModel.rxns,'NADH16'))
    nNAD = ismember(tModel.rxns,'NADH16');
    charge(nNAD)  = 0;
    h_trans(nNAD) = -4;
end

end
function dGrMatrix = createdGrMatrix(tModel)
% Generates matrix that multiplies variables in order to obtain dGr.
% No se si incluir el fluxDirection o no. Por ahora no.

exchange = findExRxns(tModel);
tModel = removeProtons(tModel); % Removing protons
St = transpose(tModel.S(:,~exchange));
Tt = transpose(tModel.T(:,~exchange));

% Parameters
R         = 8.3144626 * 1e-3;                                              % (KJ/mol K)
T         = 298.15;                                                        % (K) | 25 °C
F         = 96485;                                                         % Faraday constant C/mol
vpH = zeros(size(tModel.mets));
vpH(contains(tModel.mets,'_c')) = R*T*log(1e-7);
vpH(contains(tModel.mets,'_e')) = R*T*log(10^(-6.3));
dpH = 0.7;
potential = 1e-3 * (33.33 * dpH - 143.33);                                 % V  = J/C
potential = 1e-3 * potential;                                              % kV = kJ/C

nRxns    = size(tModel.S,2);

% dGr matrix
dGrMatrix = [zeros(size(St,1),nRxns) R*T*St St Tt*(tModel.dGf0_H + (tModel.NH .* vpH)) + F*potential*tModel.charVec(~exchange)];

end
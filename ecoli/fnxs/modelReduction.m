function modelReduction(N,vX)
%% Model Reduction
disp(strcat('*************************',num2str(N),'*',num2str(vX),'***************************************'))

load(strcat('models/COREv',num2str(vX),'.mat'),"model");

%% Add g6p exchange and make EX_pyr reversible [DELETED: No piruvate on medium]
% if N == 1
%     model = addExchangeRxn(model,'g6p_c',0,1000);
%     model.subSystems(end) = {'Extracellular exchange'};  % Correct this field
% end
% model = changeRxnBounds(model,'EX_pyr_e',-1000,'l');

%% Addition (manually) of rxns of interest and unbalanced metabolites
% Glycolysis/gluconeogenesis
glicoRxns = {'EX_glc__D_e' 'GLCpts' 'PGI' 'PFK' 'FBP' 'FBA'...
    'TPI' 'GAPD' 'PGK' 'PGM' 'ENO' 'PYK' 'PPS' 'PYRt2' 'EX_pyr_e'}';

% if N == 1 [DELETED: No virtual fructose exchange]
%     glicoRxns = {'EX_glc__D_e' 'GLCpts' 'PGI' 'PFK' 'FBP' 'FBA' 'EX_g6p_c'...
%         'TPI' 'GAPD' 'PGK' 'PGM' 'ENO' 'PYK' 'PPS' 'PYRt2' 'EX_pyr_e'}';
% else
%     glicoRxns = {'EX_glc__D_e' 'GLCpts' 'PGI' 'PFK' 'FBP' 'FBA'...
%         'TPI' 'GAPD' 'PGK' 'PGM' 'ENO' 'PYK' 'PPS' 'PYRt2' 'EX_pyr_e'}';
% end
pppRxns = model.rxns(ismember(model.subSystems,{'Pentose Phosphate Pathway'}));
m2Rxns = vertcat(glicoRxns,pppRxns);

% Rxns included in ecModelv21
ssv21 = model.rxns(ismember(model.subSystems,{'Citric Acid Cycle','Anaplerotic reactions'}));
otherv21 = {'PDH' 'SUCDi' 'FRD7'}';
rxnsv21 = vertcat(m2Rxns,ssv21,otherv21);

% Rxns included in ecModel3
newSS = model.rxns(ismember(model.subSystems,{'Citric Acid Cycle', ...
    'Oxidative Phosphorylation','Anaplerotic reactions'})); % subSystems that are now included
otherRxns = {'O2t' 'EX_o2_e' 'EX_h_e' 'CO2t' 'EX_co2_e' 'H2Ot'...
    'EX_h2o_e' 'PIt2r' 'EX_pi_e' 'PDH'}'; % 'ATPM' erased
m3Rxns = vertcat(m2Rxns,newSS,otherRxns);

% Unbalanced mets
glicoUBM = {'h2o_c' 'h_c' 'h_e' 'nad_c' 'nadh_c' 'atp_c' 'adp_c' 'amp_c' 'pi_c'};
pppUBM   = {'co2_c' 'nadp_c' 'nadph_c'};
m2UBM    = horzcat(glicoUBM,pppUBM);
v21UBM   = horzcat(glicoUBM,pppUBM,{'q8_c' 'q8h2_c'});
v3UBM    = {'q8_c' 'q8h2_c'};

%% Function call

if N == 1
    [ecModel1,ecTherModel1] = generatesubModel(model,glicoRxns,glicoUBM); % Read "easy model"
    ecModel1.description = strcat("ecModel1v",num2str(vX));
    ecTherModel1.description = strcat("ecTherModel1v",num2str(vX));
    save(strcat("models/subModels/ec1v",num2str(vX),".mat"),"ecModel1","ecTherModel1")
elseif N == 2
    [ecModel2,ecTherModel2] = generatesubModel(model,m2Rxns,m2UBM); % Read "easy model"
    ecModel2.description = strcat("ecModel2v",num2str(vX));
    ecTherModel2.description = strcat("ecTherModel2v",num2str(vX));
    save(strcat("models/subModels/ec2v",num2str(vX),".mat"),"ecModel2","ecTherModel2")
elseif N == 21 % [DISCONTINUED]
    [ecModelv21,ecTherModelv21] = generatesubModel(model,rxnsv21,v21UBM); % Read "easy model"
    ecModelv21.description = "ecModelv21";
    ecTherModelv21.description = "ecTherModelv21";
    save("models/ecModelv21.mat","ecModelv21","ecTherModelv21")
elseif N == 3
    [ecModel3,ecTherModel3] = generatesubModel(model,m3Rxns,{}); % Read "easy model"
    ecModel3.description = strcat("ecModel3v",num2str(vX));
    ecTherModel3.description = strcat("ecTherModel3v",num2str(vX));
    save(strcat("models/subModels/ec3v",num2str(vX),".mat"),"ecModel3","ecTherModel3")
else
    disp("N is invalid argument")
end


end


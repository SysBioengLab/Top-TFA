% Discrimitation graph

% Load data
load("..\files\ec3\DATAv2.mat","bigModel","boundsFcTMFA","boundsTMFA", ...
    "feasibleFcTMFA","feasibleTMFA")
% Manually changing labels
myLabel = bigModel.thermo.rxns;
myLabel{29} = 'EX\_CO2'; myLabel{30} = 'EX\_D-GLC'; myLabel{31} = 'EX\_H+';
myLabel{32} = 'EX\_H2O'; myLabel{33} = 'EX\_O2'; myLabel{34} = 'EX\_PYR';

nRxns = size(bigModel.thermo.lb,1);
affectedTop = [3 4 9 10];
nonAffectedTop = [1 2 5 6 7 8 11 12];
%% Upper bounds
f1 = figure(1);

hold on
afs  = scatter(1:nRxns,reshape(boundsFcTMFA(1:nRxns,2,affectedTop),nRxns,size(affectedTop,2))',100,[0 0.4470 0.7410],'filled','MarkerFaceAlpha',0.4);
nafs = scatter(1:nRxns,reshape(boundsFcTMFA(1:nRxns,2,nonAffectedTop),nRxns,size(nonAffectedTop,2))',[],[0.8500 0.3250 0.0980],'filled','MarkerFaceAlpha',0.4);
hold off
obs = findobj(gca); % 1st is "Axes"
legend([obs(2),obs(end)],"Non affected topologies","Affected topologies",Location="best")
ylim([min(boundsFcTMFA(1:nRxns,2,:),[],"all")-10 max(boundsFcTMFA(1:nRxns,2,:),[],"all")+10])
xticks(1:nRxns)
xticklabels(myLabel)
title("Upper bounds")

%% Lower bounds
f2 = figure(2);

hold on
afs  = scatter(1:nRxns,reshape(boundsFcTMFA(1:nRxns,1,affectedTop),nRxns,size(affectedTop,2))',100,[0 0.4470 0.7410],'filled','MarkerFaceAlpha',0.4);
nafs = scatter(1:nRxns,reshape(boundsFcTMFA(1:nRxns,1,nonAffectedTop),nRxns,size(nonAffectedTop,2))',[],[0.8500 0.3250 0.0980],'filled','MarkerFaceAlpha',0.4);
hold off
obs = findobj(gca); % 1st is "Axes"
legend([obs(2),obs(end)],"Non affected topologies","Affected topologies",Location="best")
ylim([min(boundsFcTMFA(1:nRxns,1,:),[],"all")-10 max(boundsFcTMFA(1:nRxns,1,:),[],"all")+10])
xticks(1:nRxns)
xticklabels(myLabel)
title("Lower bounds")
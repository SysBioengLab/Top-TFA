function plotMedianDGR(pathSample1,pathSample2,nTop)
% Plots median dGr for all topologies

pathModel = "results\DATA_v2.mat"; % The version does not affect
load(pathModel,"thermoRCA")
nMets = size(thermoRCA.S,1);
nRxns = size(thermoRCA.S,2);

%% Topologies
figure(1)
hold on
for i=1:nTop
    load(strcat(pathSample1,num2str(i),".mat"),"sample")
    dGr = sample.dGrMatrix * [sample.points(nRxns+1:nRxns+2*nMets,:); ones(1,sample.numSamples)];
    dGrMed = median(dGr,2)';
    histogram(dGrMed,41)
    xlim([-60,10])
end
hold off


%% Metropolis VS CHRR
if pathSample2 == ' '
    return
end
figure(2)
t = tiledlayout(2,2);
for i=1:4
    load(strcat(pathSample1,num2str(i),".mat"),"sample") % CHRR
    dGrCHRR = sample.dGrMatrix * [sample.points(nRxns+1:nRxns+2*nMets,:); ones(1,sample.numSamples)];
    dGrMedCHRR = median(dGrCHRR,2)'; clear sample;
    load(strcat(pathSample2,num2str(i),".mat"),"sample") % Metropolis
    dGrMET = sample.dGrMatrix * [sample.points(nRxns+1:nRxns+2*nMets,:); ones(1,sample.numSamples)];
    dGrMedMET = median(dGrMET,2)'; clear sample;
    % Plot
    nexttile
    hold on
    histogram(dGrMedCHRR,41)
    histogram(dGrMedMET,41)
    hold off
    caption = strcat("Topology ",num2str(i));
    title(caption,'FontSize',10)
end
t.Padding = 'compact';
t.TileSpacing = 'compact';
title(t,'median dGr distribution')


%% Metropolis VS CHRR (v2)
figure(3)
t = tiledlayout(2,2);
for i=1:4
    load(strcat(pathSample1,num2str(i),".mat"),"sample") % CHRR
    dGrCHRR = sample.dGrMatrix * [sample.points(nRxns+1:nRxns+2*nMets,:); ones(1,sample.numSamples)];
    dGrMedCHRR = median(dGrCHRR,2); clear sample;
    load(strcat(pathSample2,num2str(i),".mat"),"sample") % Metropolis
    dGrMET = sample.dGrMatrix * [sample.points(nRxns+1:nRxns+2*nMets,:); ones(1,sample.numSamples)];
    dGrMedMET = median(dGrMET,2); clear sample;
    % Plot
    nexttile
    hold on
    bar([dGrMedCHRR dGrMedMET])
    hold off
    caption = strcat("Topology ",num2str(i));
    title(caption,'FontSize',10)
end
t.Padding = 'compact';
t.TileSpacing = 'compact';
title(t,'median dGr distribution')
lg = legend('US','MHR');
lg.Layout.Tile = 'North';

end
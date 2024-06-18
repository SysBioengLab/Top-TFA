% Flux graph

% Load data for selected topologies
load("..\files\ec3\DATAv2.mat","bigModel","boundsFcTMFA","boundsTMFA", ...
    "feasibleFcTMFA","feasibleTMFA") % We use the one with more information

% Manually changing labels
myLabel = bigModel.thermo.rxns;
myLabel{29} = 'EX\_CO2'; myLabel{30} = 'EX\_D-GLC'; myLabel{31} = 'EX\_H+';
myLabel{32} = 'EX\_H2O'; myLabel{33} = 'EX\_O2'; myLabel{34} = 'EX\_PYR';

%% Grafico de flujos (Version 1)
selectedTop = [3 4 9 10];
f = figure();
t = tiledlayout(2,2);
nRxns = size(bigModel.thermo.lb,1);
for i=selectedTop
    nexttile
    hold on
    boxchart([boundsTMFA(1:nRxns,1,i)'; boundsTMFA(1:nRxns,2,i)'],"BoxMedianLineColor",'none','BoxWidth',1);
    boxchart([boundsFcTMFA(1:nRxns,1,i)'; boundsFcTMFA(1:nRxns,2,i)'],"BoxMedianLineColor",'none','BoxWidth',1);
    hold off
    ax = gca;
    ax.FontSize = 6; 
    xticklabels(myLabel)
    xtickangle(45)
    caption = strcat("Topology ",num2str(i));
    title(caption,'FontSize',8)
end
t.Padding = 'compact';
t.TileSpacing = 'compact';
leg = legend("TFA","tc-TMFA",'Orientation', 'Horizontal','Fontsize',8);
leg.Layout.Tile = 'north';
title(t,'tc-TMFA VS TFA')

%% Grafico de flujos (Version 2)

% f = figure();
% nRxns = size(bigModel.thermo.lb,1);
% hold on
% for i=selectedTop
%     boxchart([boundsTMFA(1:nRxns,1,i)'; boundsTMFA(1:nRxns,2,i)'],"BoxMedianLineColor",'none','BoxWidth',1);
%     boxchart([boundsFcTMFA(1:nRxns,1,i)'; boundsFcTMFA(1:nRxns,2,i)'],"BoxMedianLineColor",'none','BoxWidth',1);
% end
% hold off
% ax = gca;
% ax.FontSize = 6; 
% xticklabels(myLabel)
% xtickangle(45)
% leg = legend("TFA-Top3","tc-TMFA-Top3","TFA-Top4","tc-TMFA-Top4","TFA-Top9","tc-TMFA-Top9", ...
%     "TFA-Top10","tc-TMFA-Top10",'Fontsize',8,'Location','Best');
% title('tc-TMFA VS TFA')

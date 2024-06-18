% dGr graph

% Load data for selected topologies
load("..\files\ec3\DATAv2.mat","bigModel","dGrboundsTMFA","dGrboundsFcTMFA", ...
    "feasibleFcTMFA","feasibleTMFA")

feasibleTopologies = 1:size(feasibleFcTMFA,1);
if isequal(feasibleFcTMFA,feasibleTMFA)
    dGrboundsFcTMFA(:,:,~feasibleFcTMFA) = [];
    dGrboundsTMFA(:,:,~feasibleTMFA) = [];
    feasibleTopologies(~feasibleFcTMFA) = [];
else
    newFeasible = logical(feasibleFcTMFA .* feasibleTMFA);
    dGrboundsFcTMFA(:,:,~newFeasible) = [];
    dGrboundsTMFA(:,:,~newFeasible) = [];
    feasibleTopologies(~newFeasible) = [];
end

% Manually changing labels
exchange = findExRxns(bigModel.thermo);
myLabel = bigModel.thermo.rxns(~exchange);

%% Plot 1
f = figure(1);
t = tiledlayout(3,4);

for i=feasibleTopologies
    nexttile
    hold on
    boxchart(dGrboundsTMFA(:,:,i)',"BoxMedianLineColor",'none','BoxWidth',1);
    boxchart(dGrboundsFcTMFA(:,:,i)',"BoxMedianLineColor",'none','BoxWidth',1);
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
title(t,'dGr (kJ/mol) tc-TMFA VS TFA')

%% Plot 2
f = figure(2);
t = tiledlayout(1,2);

% TFA
nexttile
hold on
for i=feasibleTopologies
    boxchart(dGrboundsTMFA(:,:,i)',"BoxMedianLineColor",'none','BoxWidth',1);
end
hold off
ax = gca;
ax.FontSize = 6; 
xticklabels(myLabel)
xtickangle(45)
title("TFA",'FontSize',8)

% tcTMFA
nexttile
hold on
for i=feasibleTopologies
    boxchart(dGrboundsFcTMFA(:,:,i)',"BoxMedianLineColor",'none','BoxWidth',1);
end
hold off
ax = gca;
ax.FontSize = 6; 
xticklabels(myLabel)
xtickangle(45)
title("tcTMFA",'FontSize',8)
% General
t.Padding = 'compact';
t.TileSpacing = 'compact';
leg = legend("Topology 1","Topology 2","Topology 3","Topology 4", ...
    "Topology 5","Topology 6","Topology 7","Topology 8","Topology 9", ...
    "Topology 10","Topology 11","Topology 12",'Orientation', 'Vertical','Fontsize',8);
leg.Layout.Tile = 'east';
title(t,'dGr (kJ/mol)')

%% Plot 3 (segregation | only tc-TMFA)
f = figure(3);

affectedTop = [3 4 9 10];
nonAffectedTop = [1 2 5 6 7 8 11 12];
hold on
boxchart(reshape(dGrboundsFcTMFA(:,:,affectedTop),size(dGrboundsFcTMFA,1),2*size(affectedTop,2))',"BoxMedianLineColor",'none','BoxWidth',1);
boxchart(reshape(dGrboundsFcTMFA(:,:,nonAffectedTop),size(dGrboundsFcTMFA,1),2*size(nonAffectedTop,2))',"BoxMedianLineColor",'none','BoxWidth',1);
hold off
ax = gca;
ax.FontSize = 6; 
xticklabels(myLabel)
xtickangle(45)
title("dGr (kJ/mol) by group",'FontSize',8)
legend("Affected","Non affected","Location","best")

%% Plot 4 (segregation | only tc-TMFA)
f = figure(4);

affectedTop = [3 4 9 10];
nonAffectedTop = [1 2 5 6 7 8 11 12];
hold on
scatter(1:size(dGrboundsFcTMFA,1),reshape(dGrboundsFcTMFA(:,:,affectedTop),size(dGrboundsFcTMFA,1),2*size(affectedTop,2))',[],[0 0.4470 0.7410],'filled','MarkerFaceAlpha',0.4);
scatter(1:size(dGrboundsFcTMFA,1),reshape(dGrboundsFcTMFA(:,:,nonAffectedTop),size(dGrboundsFcTMFA,1),2*size(nonAffectedTop,2))',[],[0.8500 0.3250 0.0980],'filled','MarkerFaceAlpha',0.4);
hold off
ax = gca;
ax.FontSize = 6; 
xticklabels(myLabel)
xtickangle(45)
title("dGr (kJ/mol) by group",'FontSize',8)
legend("Affected","Non affected","Location","best")

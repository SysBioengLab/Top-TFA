% Figure 5

%% Load data
% A
load("..\files\ec3\DATAv2.mat","bigModel","boundsFcTMFA","boundsTMFA","feasibleFcTMFA","feasibleTMFA")
% B
load("..\DGR_cases.mat","case_2","dGrRxns")
load("..\files\henry2007_dGr.mat","M")
% C
sTC   = "..\samples\topTFA_ADS_pro.mat";
sTM   = "..\samples\TMFA_ADS_pro.mat";
load(sTC,"sample")
sampleTC = sample; clear sample;
load(sTM,"sample")
sampleTM = sample; clear sample;

%% Parameters definition
% C
model = bigModel.thermo;
nF    = [6 17];


%% General calculations
% Bounds (A)
nRxns = size(bigModel.thermo.S,2);
nMets = size(bigModel.thermo.S,1);
boundsFcTMFA = boundsFcTMFA(1:nRxns+nMets,:,:);
boundsTMFA = boundsTMFA(1:nRxns+nMets,:,:);

Topologies = 1:size(feasibleFcTMFA,1);
if isequal(feasibleFcTMFA,feasibleTMFA)
    boundsFcTMFA(:,:,~feasibleFcTMFA) = [];
    boundsTMFA(:,:,~feasibleTMFA) = [];
    Topologies(~feasibleFcTMFA) = [];
else
    newFeasible = logical(feasibleFcTMFA .* feasibleTMFA);
    boundsFcTMFA(:,:,~newFeasible) = [];
    boundsTMFA(:,:,~newFeasible) = [];
    Topologies(~newFeasible) = [];
end

% Ratio calculation
ratio = zeros(size(boundsFcTMFA,1),size(boundsFcTMFA,3));
for i=1:size(ratio,2)
    ratio(:,i) = abs(boundsFcTMFA(:,1,i)-boundsFcTMFA(:,2,i))./abs(boundsTMFA(:,1,i)-boundsTMFA(:,2,i));
end
ratio(isinf(ratio)|isnan(ratio)) = 1; % This happens when lb=ub, so ratio=1 assumed

% dGr graph (B)
% Remove H2Ot
% wtr_indx = find(contains(dGrRxns,'H2Ot'));
banRxns = find(M(:,1)==0);
indxVec = (1:size(case_2,1))';
case_2(banRxns,:) = [];
M(banRxns,:) = [];
indxVec(banRxns) = [];

% Reorder
orderM = [median(case_2,2) case_2 M indxVec];
orderM = sortrows(orderM);
c2 = orderM(:,2:3)';
henryM = orderM(:,4:5)';
ordererdGrRxns = dGrRxns(orderM(:,6));

% dGr (C)
avoidIndx = [3 4]; % Hardcode, index of water concentration
dGrMatrixTC = sampleTC.dGrMatrix;
dGrMatrixTM = sampleTM.dGrMatrix;
dGrMatrixTC(:,avoidIndx) = [];
dGrMatrixTM(:,avoidIndx) = [];
dGrTC = dGrMatrixTC * [sampleTC.points(nRxns+1:nRxns+2*nMets-2,:); ones(1,size(sampleTC.points,2))];
dGrTM = dGrMatrixTM * [sampleTM.points(nRxns+1:nRxns+2*nMets-2,:); ones(1,size(sampleTM.points,2))];
dGrTC = -dGrTC(nF,:); % -dGr
dGrTM = -dGrTM(nF,:); % -dGr

% Limits calculation
grStar = 4.3478;
exchange = findExRxns(model);
Vmax = model.vmax;
Vmax(exchange) = [];
Vmax = Vmax(nF);
tcRegionX = 0:0.001:4.3478;
tcRegionX12 = [tcRegionX grStar:0.1:max(dGrTM(1,:))];
tcRegionX22 = [tcRegionX grStar:0.1:max(dGrTM(2,:))];
tcRegionY12 = [tcRegionX*Vmax(1)/grStar Vmax(1)*ones(1,size(grStar:0.1:max(dGrTM(1,:)),2))];
tcRegionY22 = [tcRegionX*Vmax(2)/grStar Vmax(2)*ones(1,size(grStar:0.1:max(dGrTM(2,:)),2))];
tmfaRegionX1 = ones(1,size(0:1:max(dGrTM(1,:)),2));
tmfaRegionX2 = ones(1,size(0:1:max(dGrTM(2,:)),2));
tmfaRegionY1 = Vmax(1) .* tmfaRegionX1;
tmfaRegionY2 = Vmax(2) .* tmfaRegionX2;

%% Plots
figure(1)
T = tiledlayout(2,3,"TileIndexing","columnmajor"); % ,"TileSpacing","tight","Padding","compact"
% A
ax1 = nexttile(T);
% boxplot(ratio(1:nRxns,:),'BoxStyle','outline','Colors','k','MedianStyle','line','Jitter',.5,'Symbol','xb')
boxchart(ratio(1:nRxns,:),"JitterOutliers","on","WhiskerLineColor",[0 0.4470 0.7410])
ylim([0.5 1.05])
set(gca,'xticklabel',Topologies)
xlabel('Topology')
ylabel('Ratio of flux ranges')
ax1.FontSize = 14;
title(ax1,'A',"FontSize",18,"FontWeight","bold")
ax1.TitleHorizontalAlignment = 'left';

% B
ax2 = nexttile(T);
boxchart(c2,'BoxWidth',1)
hold on
% boxchart(c3,'BoxWidth',1)
boxchart(henryM,'BoxWidth',1)
hold off
ylim([-230 60]);
xticks([])
xlabel("Reactions")
ylabel("\Delta_{r} G' (kJ/mol)")
ax2.FontSize = 14;
% legend({'Normal' 'Flux-fixed' 'Henry2007'},'Location','best','FontSize',12)
legend({'Top-TFA' 'Henry2007'},'Location','best','FontSize',12)
title(ax2,'B',"FontSize",18,"FontWeight","bold")
ax2.TitleHorizontalAlignment = 'left';

% C
t = tiledlayout(T,2,2);
t.Layout.Tile = 3;
t.Layout.TileSpan = [2 2];
t.TileSpacing = "tight";
t.Padding = "tight";
t.TileIndexing = "columnmajor";

% TMFA - RPE
ax3 = nexttile(t);
dscatter(dGrTM(1,:)',sampleTM.points(nF(1),:)','marker','o','filled',false)
hold on
yline(Vmax(1),"LineWidth",1.5,"Color",[0 0 0])
hold off
ax3.FontSize = 13;
title(strcat("TMFA: ",dGrRxns(nF(1))),"FontSize",12,"FontWeight","normal")

% TMFA - RPE
ax4 = nexttile(t);
dscatter(dGrTM(2,:)',sampleTM.points(nF(2),:)','marker','o','filled',false)
hold on
yline(Vmax(2),"LineWidth",1.5,"Color",[0 0 0])
hold off
ax4.FontSize = 13;
title(strcat("TMFA: ",dGrRxns(nF(2))),"FontSize",12,"FontWeight","normal")

% Top - PPC
ax5 = nexttile(t);
dscatter(dGrTC(1,:)',sampleTC.points(nF(1),:)','marker','o','filled',false)
hold on
plot(tcRegionX12,tcRegionY12,"LineWidth",1.5,"Color",[0 0 0])
hold off
ax5.FontSize = 13;
title(strcat("Top-TMFA: ",dGrRxns(nF(1))),"FontSize",12,"FontWeight","normal")

% Top - RPE
ax6 = nexttile(t);
dscatter(dGrTC(2,:)',sampleTC.points(nF(2),:)','marker','o','filled',false)
hold on
plot(tcRegionX22,tcRegionY22,"LineWidth",1.5,"Color",[0 0 0])
hold off
ax6.FontSize = 13;
title(strcat("Top-TMFA: ",dGrRxns(nF(2))),"FontSize",12,"FontWeight","normal")

cbh = colorbar;
cbh.Layout.Tile = 'east';
linkaxes([ax3 ax5],'xy')
linkaxes([ax4 ax6],'xy')
ax3.YLim = [0 Vmax(1)*1.05];
xlabel(t,"-\DeltaG_{r} (kJ/mol)","FontSize",14)
ylabel(t,"Flux (mmol/gdcw/h)","FontSize",14)

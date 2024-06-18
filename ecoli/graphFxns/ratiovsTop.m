function ratiovsTop()
% Ratio vs Topology graph

% Load data
load("..\files\ec3\DATAv2.mat","bigModel","boundsFcTMFA","boundsTMFA","feasibleFcTMFA","feasibleTMFA")

% Bounds
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

% Plot (double)
% figure(1)
% t = tiledlayout(1,2);
% nexttile
% boxplot(ratio(1:nRxns,:),'BoxStyle','outline','Colors','k','MedianStyle','line')
% set(gca,'xticklabel',Topologies)
% title('Fluxes')
% xlabel('Topology')
% ylabel('Ratio')
% 
% nexttile
% boxplot(ratio(nRxns+1:nRxns+nMets,:),'BoxStyle','outline','Colors','k','MedianStyle','line')
% set(gca,'xticklabel',Topologies)
% title('Concentrations')
% xlabel('Topology')
% ylabel('Ratio')
% title(t,'Ratio between tc-TMFA and TMFA')

%% Real Fluxes
figure(1) % (single)
boxplot(ratio(1:nRxns,:),'BoxStyle','outline','Colors','k','MedianStyle','line','Jitter',.5,'Symbol','xb')
set(gca,'xticklabel',Topologies)
xlabel('Topology')
ylabel('Ratio between fluxes ranges (Top-TMFA/TMFA)')

%% Real concentrations
figure(2)
boxplot(ratio(nRxns+1:nRxns+nMets,:),'BoxStyle','outline','Colors','k','MedianStyle','line','Jitter',.5,'Symbol','xb')
set(gca,'xticklabel',Topologies)
title('Ratio between concentrations ranges of top-TMFA and TMFA')
xlabel('Topology')
ylabel('Ratio')
% %% Violin
% if ~isViolin
%     figure(1) % (single)
%     boxplot(ratio(1:nRxns,:),'BoxStyle','outline','Colors','k','MedianStyle','line','Jitter',.5,'Symbol','xb')
%     set(gca,'xticklabel',Topologies)
%     title('Ratio between fluxes ranges of top-TMFA and TMFA')
%     xlabel('Topology')
%     ylabel('Ratio')
% else
%     figure(1) % (single)
%     [h,L,MX,MED,bw] = violin(ratio(1:nRxns,:));
%     % boxplot(ratio(1:nRxns,:),'BoxStyle','outline','Colors','k','MedianStyle','line')
%     set(gca,'xticklabel',Topologies)
%     title('Ratio between fluxes ranges of top-TMFA and TMFA')
%     xlabel('Topology')
%     ylabel('Ratio')
% end

end
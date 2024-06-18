function plotDGRvsflux(model,sTC,sTM,nF)

load(sTC,"sample")
sampleTC = sample; clear sample;

load(sTM,"sample")
sampleTM = sample; clear sample;

nRxns = size(model.S,2);
nMets = size(model.S,1);
% dGr calculation
dGrTC = sampleTC.dGrMatrix * [sampleTC.points(nRxns+1:nRxns+2*nMets,:); ones(1,sampleTC.numSamples)];
dGrTM = sampleTM.dGrMatrix * [sampleTM.points(nRxns+1:nRxns+2*nMets,:); ones(1,sampleTM.numSamples)];
dGrTC = -dGrTC(nF,:); % -dGr
dGrTM = -dGrTM(nF,:); % -dGr

% Limits calculation
grStar = 4.3478;
exchange = findExRxns(model);
Vmax = model.vmax;
Vmax(exchange) = [];
Vmax = Vmax(nF);
dGrRxns = model.rxns(~exchange);
tcRegionX = 0:0.001:4.3478;
tcRegionX11 = [tcRegionX grStar:0.1:max(dGrTC(1,:))];
tcRegionX12 = [tcRegionX grStar:0.1:max(dGrTM(1,:))];
tcRegionX21 = [tcRegionX grStar:0.1:max(dGrTC(2,:))];
tcRegionX22 = [tcRegionX grStar:0.1:max(dGrTM(2,:))];
tcRegionY11 = [tcRegionX*Vmax(1)/grStar Vmax(1)*ones(1,size(grStar:0.1:max(dGrTC(1,:)),2))];
tcRegionY12 = [tcRegionX*Vmax(1)/grStar Vmax(1)*ones(1,size(grStar:0.1:max(dGrTM(1,:)),2))];
tcRegionY21 = [tcRegionX*Vmax(2)/grStar Vmax(2)*ones(1,size(grStar:0.1:max(dGrTC(2,:)),2))];
tcRegionY22 = [tcRegionX*Vmax(2)/grStar Vmax(2)*ones(1,size(grStar:0.1:max(dGrTM(2,:)),2))];
% tcRegionY2 = tcRegionX*Vmax(2)/grStar;

% Plot
figure(1)
tiledlayout(2,2,"TileSpacing","tight","Padding","tight");
ax1 = nexttile;
dscatter(dGrTM(1,:)',sampleTM.points(nF(1),:)','marker','o','filled',false)
hold on
plot(tcRegionX12,tcRegionY12,"LineWidth",1)
% plot(grStar:0.1:max(dGrTM(1,:)),Vmax(1)*ones(1,size(grStar:0.1:max(dGrTM(1,:)),2)),"LineWidth",1)
hold off
title(strcat("TMFA: ",dGrRxns(nF(1))))
xlabel("-\DeltaG_{r} (kJ/mol)","FontSize",10)
ylabel("Flux (mmol/gdcw/h)","FontSize",10)
ax2 = nexttile;
dscatter(dGrTM(2,:)',sampleTM.points(nF(2),:)','marker','o','filled',false)
hold on
plot(tcRegionX22,tcRegionY22,"LineWidth",1)
% plot(grStar:0.1:max(dGrTM(2,:)),Vmax(2)*ones(1,size(grStar:0.1:max(dGrTM(2,:)),2)),"LineWidth",1)
hold off
title(strcat("TMFA: ",dGrRxns(nF(2))))
xlabel("-\DeltaG_{r} (kJ/mol)","FontSize",10)
ylabel("Flux (mmol/gdcw/h)","FontSize",10)
ax3 = nexttile;
dscatter(dGrTC(1,:)',sampleTC.points(nF(1),:)','marker','o','filled',false)
% hold on
% plot(tcRegionX11,tcRegionY11,"LineWidth",1)
% plot(grStar:0.1:max(dGrTC(1,:)),Vmax(1)*ones(1,size(grStar:0.1:max(dGrTC(1,:)),2)),"LineWidth",1)
% hold off
title(strcat("Top-TMFA: ",dGrRxns(nF(1))))
xlabel("-\DeltaG_{r} (kJ/mol)","FontSize",10)
ylabel("Flux (mmol/gdcw/h)","FontSize",10)
ax4 = nexttile;
dscatter(dGrTC(2,:)',sampleTC.points(nF(2),:)','marker','o','filled',false)
% hold on
% plot(tcRegionX21,tcRegionY21,"LineWidth",1)
% plot(grStar:0.1:max(dGrTC(2,:)),Vmax(2)*ones(1,size(grStar:0.1:max(dGrTC(2,:)),2)),"LineWidth",1)
% hold off
title(strcat("Top-TMFA: ",dGrRxns(nF(2))))
xlabel("-\DeltaG_{r} (kJ/mol)","FontSize",10)
ylabel("Flux (mmol/gdcw/h)","FontSize",10)
cbh = colorbar;
cbh.Layout.Tile = 'east';
linkaxes([ax1 ax3],'xy')
linkaxes([ax2 ax4],'xy')

end
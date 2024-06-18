function plotComparison(model,sampleName1,sampleName2,dir)
% Genera los graficos de flujo de reacciones, concentración de metabolitos
% y valor de dGr. Agrega el campo dGr a sample.

load(strcat("samples\tests\",sampleName1,".mat"),"sample")
sample1 = sample; clear sample;
load(strcat("samples\tests\",sampleName2,".mat"),"sample")
sample2 = sample; clear sample;

% Inicialización de parametros
nRxns = size(model.S,2);
nMets = size(model.S,1) - 2; % Hardcode: water concentration is removed with ADSB
ndGr  = size(sample1.dGrRxns,1);
metNames = model.metNames(~(ismember(model.metNames,'H2O')));
exchange = findExcRxns(model);
bigg_dGr = model.BiGGID(~exchange);

% Grafico de flujos
f1 = figure(1);
t = tiledlayout(5,8);
for i=1:nRxns
    nexttile
    histogram(sample1.points(i,:),40,'DisplayStyle','stairs',"EdgeColor",[0 0.4470 0.7410])
    hold on
    histogram(sample1.points(i,:),40,"EdgeColor","none","FaceColor",[0 0.4470 0.7410])
    histogram(sample2.points(i,:),40,'DisplayStyle','stairs',"EdgeColor",[0.8500 0.3250 0.0980])
    histogram(sample2.points(i,:),40,"EdgeColor","none","FaceColor",[0.8500 0.3250 0.0980])
    hold off
    % caption = createCaption(model.rxnNames{i});
    caption = correct_TEXT(model.BiGGID{i});
    title(caption,'FontSize',10)
end
t.Padding = 'compact';
t.TileSpacing = 'compact';
title(t,'Flux distribution (mmol /(gdcw h))')
lg = legend('','Uniform','','Non Uniform');
lg.Layout.Tile = 'North';

% Grafico de concentraciones
f2 = figure(2);
t = tiledlayout(6,7);
for i=1:nMets
    nexttile
    histogram(1e3 * exp(sample1.points(nRxns+i,:)),40,'DisplayStyle','stairs',"EdgeColor",[0 0.4470 0.7410])                      % mM (mmol)
    hold on
    histogram(1e3 * exp(sample1.points(nRxns+i,:)),40,"EdgeColor","none","FaceColor",[0 0.4470 0.7410])
    histogram(1e3 * exp(sample2.points(nRxns+i,:)),40,'DisplayStyle','stairs',"EdgeColor",[0.8500 0.3250 0.0980])
    histogram(1e3 * exp(sample2.points(nRxns+i,:)),40,"EdgeColor","none","FaceColor",[0.8500 0.3250 0.0980])
    hold off
    caption = createCaption(metNames{i});
    title(caption,'FontSize',8)
end
t.Padding = 'compact';
t.TileSpacing = 'compact';
title(t,'Concentration distributions (mmol)')
lg = legend('','Uniform','','Non Uniform');
lg.Layout.Tile = 'North';

% Grafico de dGr
% CAMBIO
% dGr1 = sample1.dGrMatrix * [sample1.points(nRxns+1:nRxns+2*nMets,:); ones(1,sample1.numSamples)];
% dGr2 = sample2.dGrMatrix * [sample2.points(nRxns+1:nRxns+2*nMets,:); ones(1,sample2.numSamples)];
avoidIndx = [3 34]; % Hardcode, index of water concentration
dGrMatrix1 = sample1.dGrMatrix;
dGrMatrix2 = sample2.dGrMatrix;
dGrMatrix1(:,avoidIndx) = [];
dGrMatrix2(:,avoidIndx) = [];
dGr1 = dGrMatrix1 * [sample1.points(nRxns+1:nRxns+2*nMets+2,:); ones(1,size(sample1.points,2))];
dGr2 = dGrMatrix2 * [sample2.points(nRxns+1:nRxns+2*nMets+2,:); ones(1,size(sample2.points,2))];

h2oIndex = find(contains(sample1.dGrRxns,'r2030'));
f3 = figure(3);
t = tiledlayout(3,10);
for i=1:ndGr
    if ~isequal(i,h2oIndex)
        nexttile
        histogram(dGr1(i,:),40,'DisplayStyle','stairs',"EdgeColor",[0 0.4470 0.7410])
        hold on
        histogram(dGr1(i,:),40,"EdgeColor","none","FaceColor",[0 0.4470 0.7410])
        histogram(dGr2(i,:),40,'DisplayStyle','stairs',"EdgeColor",[0.8500 0.3250 0.0980])
        histogram(dGr2(i,:),40,"EdgeColor","none","FaceColor",[0.8500 0.3250 0.0980])
        hold off
        caption = correct_TEXT(bigg_dGr{i});
        title(caption,'FontSize',8)
    end
end
t.Padding = 'compact';
t.TileSpacing = 'compact';
title(t,'\Delta_r G distributions (kJ/mol)')
lg = legend('','Uniform','','Non Uniform');
lg.Layout.Tile = 'North';

% Guardar archivos
saveas(f1,strcat(dir,sampleName1,"VS",sampleName2,"_fluxes"))
saveas(f2,strcat(dir,sampleName1,"VS",sampleName2,"_concentrations"))
saveas(f3,strcat(dir,sampleName1,"VS",sampleName2,"_dGr"))
%% HARDCODE grafico agrupado
% disp("Generating joint plots")
% nRxns = size(model.S,2);
% nMets = size(model.S,1);
% 
% load("samples\sample_check_top_1.mat","sample")
% sample1 = sample; clear sample;
% load("samples\sample_check_top_2.mat","sample")
% sample2 = sample; clear sample;
% load("samples\sample_check_top_3.mat","sample")
% sample3 = sample; clear sample;
% load("samples\sample_check_top_4.mat","sample")
% sample4 = sample; clear sample;
% 
% ndGr  = size(sample1.dGrRxns,1);
% 
% dGr1 = sample1.dGrMatrix * [sample1.points(nRxns+1:nRxns+2*nMets,:); ones(1,sample1.numSamples)];
% dGr2 = sample2.dGrMatrix * [sample2.points(nRxns+1:nRxns+2*nMets,:); ones(1,sample2.numSamples)];
% dGr3 = sample3.dGrMatrix * [sample3.points(nRxns+1:nRxns+2*nMets,:); ones(1,sample3.numSamples)];
% dGr4 = sample4.dGrMatrix * [sample4.points(nRxns+1:nRxns+2*nMets,:); ones(1,sample4.numSamples)];
% bigDGR1 = [dGr1 dGr2 dGr3 dGr4];
% bigSample1.points = [sample1.points sample2.points sample3.points sample4.points];
% clear sample1 sample2 sample3 sample4 dGr1 dGr2 dGr3 dGr4

% load("results\sample_noMedCHRR_top_1.mat","sample")
% sample1 = sample; clear sample;
% load("results\sample_noMedCHRR_top_2.mat","sample")
% sample2 = sample; clear sample;
% load("results\sample_noMedCHRR_top_3.mat","sample")
% sample3 = sample; clear sample;
% load("results\sample_noMedCHRR_top_4.mat","sample")
% sample4 = sample; clear sample;

% dGr1 = sample1.dGrMatrix * [sample1.points(nRxns+1:nRxns+2*nMets,:); ones(1,sample1.numSamples)];
% dGr2 = sample2.dGrMatrix * [sample2.points(nRxns+1:nRxns+2*nMets,:); ones(1,sample2.numSamples)];
% dGr3 = sample3.dGrMatrix * [sample3.points(nRxns+1:nRxns+2*nMets,:); ones(1,sample3.numSamples)];
% dGr4 = sample4.dGrMatrix * [sample4.points(nRxns+1:nRxns+2*nMets,:); ones(1,sample4.numSamples)];
% bigDGR2 = [dGr1 dGr2 dGr3 dGr4];
% bigSample2.points = [sample1.points sample2.points sample3.points sample4.points];
% clear sample1 sample2 sample3 sample4 dGr1 dGr2 dGr3 dGr4

% Grafico de flujos
% f1 = figure(1);
% t = tiledlayout(5,8);
% for i=1:nRxns
%     nexttile
%     histogram(sample1.points(i,:),20)
%     hold on
%     histogram(sample2.points(i,:),20)
%     histogram(sample3.points(i,:),20)
%     histogram(sample4.points(i,:),20)
%     hold off
%     caption = createCaption(model.rxnNames{i});
%     title(caption,'FontSize',8)
% end
% t.Padding = 'compact';
% t.TileSpacing = 'compact';
% title(t,'Flujos')
% lg = legend('top1','top2','top3','top4');
% lg.Layout.Tile = 'North';
% 
% % Grafico de concentraciones
% f2 = figure(2);
% t = tiledlayout(6,7);
% for i=1:nMets
%     nexttile
%     histogram(1e3 * exp(sample1.points(nRxns+i,:)),20)                      % mM (mmol)
%     hold on
%     histogram(1e3 * exp(sample2.points(nRxns+i,:)),20)                      % mM (mmol)
%     histogram(1e3 * exp(sample3.points(nRxns+i,:)),20)                      % mM (mmol)
%     histogram(1e3 * exp(sample4.points(nRxns+i,:)),20)                      % mM (mmol)
%     hold off
%     caption = createCaption(model.metNames{i});
%     title(caption,'FontSize',8)
% end
% t.Padding = 'compact';
% t.TileSpacing = 'compact';
% title(t,'Concentraciones (mmol)')
% lg = legend('top1','top2','top3','top4');
% lg.Layout.Tile = 'North';
% 
% % Grafico de dGr
% dGrRxns = model.rxnNames(~findExRxns(model));
% 
% f3 = figure(3);
% t = tiledlayout(3,10);
% for i=1:ndGr
%     nexttile
%     histogram(dGr1(i,:),20)
%     hold on
%     histogram(dGr2(i,:),20)
%     histogram(dGr3(i,:),20)
%     histogram(dGr4(i,:),20)
%     hold off
%     caption = createCaption(dGrRxns{i});
%     title(caption,'FontSize',8)
% end
% t.Padding = 'compact';
% t.TileSpacing = 'compact';
% title(t,'dGr')
% lg = legend('top1','top2','top3','top4');
% lg.Layout.Tile = 'North';
% 
% % Guardar archivos (HARDCODE)
% saveas(f1,strcat(dir,"check_flux"))
% saveas(f2,strcat(dir,"check_conc"))
% saveas(f3,strcat(dir,"check_DGR"))

end
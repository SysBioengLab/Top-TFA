function plotFig6(model,sampleName1,sampleName2)
% Generate graph for figure 6

load(strcat("samples\tests\",sampleName1,".mat"),"bigUni")
sample1 = bigUni; clear bigUni;
load(strcat("samples\tests\",sampleName2,".mat"),"bigNouni")
sample2 = bigNouni; clear bigNouni;

% Inicialización de parametros
nRxns = size(model.S,2);
nMets = size(model.S,1) - 2; % Hardcode: water concentration is removed with ADSB
metNames = model.metNames(~(ismember(model.metNames,'H2O')));
exchange = findExcRxns(model);
bigg_dGr = model.BiGGID(~exchange);

% Calculo dGr
avoidIndx = [3 34]; % Hardcode, index of water concentration
dGrMatrix1 = sample1.dGrMatrix;
dGrMatrix2 = sample2.dGrMatrix;
dGrMatrix1(:,avoidIndx) = [];
dGrMatrix2(:,avoidIndx) = [];
dGr1 = dGrMatrix1 * [sample1.points(nRxns+1:nRxns+2*nMets+2,:); ones(1,size(sample1.points,2))];
dGr2 = dGrMatrix2 * [sample2.points(nRxns+1:nRxns+2*nMets+2,:); ones(1,size(sample2.points,2))];

% First-row calculations
mf1 = prctile(sample1.points(1:nRxns,:),50,2);
mf2 = prctile(sample2.points(1:nRxns,:),50,2);
mc1 = prctile(log10(1e3*exp(sample1.points(nRxns+1:nRxns+nMets,:))),50,2);
mc2 = prctile(log10(1e3*exp(sample2.points(nRxns+1:nRxns+nMets,:))),50,2);
md1 = prctile(dGr1,50,2);
md2 = prctile(dGr2,50,2);

% Linear regression
modF = fitlm(mf1,mf2);
modC = fitlm(mc1,mc2);
modD = fitlm(md1,md2);

% Selected variables
nvars = 3;
sel_fluxes = [9 10 19];
sel_conc   = [10 24 23];
sel_dGr    = [3 8 21];

% Matrix of points
s1   = sample1.points(sel_fluxes,:);
s2   = sample2.points(sel_fluxes,:);
s11  = sample1.points(nRxns+sel_conc,:);
s22  = sample2.points(nRxns+sel_conc,:);
dg1 = dGr1(sel_dGr,:);
dg2 = dGr2(sel_dGr,:);

% Graph
figure(1);
t = tiledlayout(4,3,TileIndexing="columnmajor");
% Fluxes
nexttile
c1 = modF.Coefficients.Estimate(1); c2 = modF.Coefficients.Estimate(2); R = modF.Rsquared.Ordinary;
plot(mf1,mf2,'x',"Color",[0.7490, 0.2510, 0.7490])
l1 = 1.1*min([mf1;mf2]);
l2 = 1.1*max([mf1;mf2]);
hold on
plot(l1:l2,l1:l2,"Color",[0 0 0])
plot(l1:l2,c1+(l1:l2)*c2,"Color",[0.7490, 0.2510, 0.7490])
hold off
text(.05,.9,strcat('y = ',num2str(c1),' + ',num2str(c2),' x'),"Units","normalized","FontSize",9);
text(.05,.75,strcat("R^{2} = ",num2str(R)),"Units","normalized","FontSize",9);
xlim([l1 l2])
ylim([l1 l2])
title('Fluxes','FontSize',14)

for i=1:nvars
    nexttile
    h1 = histogram(sample1.points(sel_fluxes(i),:),40,'DisplayStyle','stairs',"EdgeColor",[0 0.4470 0.7410],"Normalization","pdf");
    hold on
    histogram(sample1.points(sel_fluxes(i),:),40,"EdgeColor","none","FaceColor",[0 0.4470 0.7410],"Normalization","pdf")
    h2 = histogram(sample2.points(sel_fluxes(i),:),40,'DisplayStyle','stairs',"EdgeColor",[0.8500 0.3250 0.0980],"Normalization","pdf");
    histogram(sample2.points(sel_fluxes(i),:),40,"EdgeColor","none","FaceColor",[0.8500 0.3250 0.0980],"Normalization","pdf")
    hold off
    ylim([0 1.1*max([h1.Values h2.Values])])
    caption = correct_TEXT(model.BiGGID{sel_fluxes(i)});
    title(caption,'FontSize',10)
    if i ==3
        p = ranksum(s1(i,:),s2(i,:));
        text(.05,.9,strcat('p = ',num2str(p)),"Units","normalized","FontSize",9);
        xlabel("mmol/(gdcw h)","FontSize",12)
    else
        p = ranksum(s1(i,:),s2(i,:));
        text(.6,.9,strcat('p = ',num2str(p)),"Units","normalized","FontSize",9);
    end
end

% Concentrations
nexttile
c1 = modC.Coefficients.Estimate(1); c2 = modC.Coefficients.Estimate(2); R = modC.Rsquared.Ordinary;
plot(mc1,mc2,'x',"Color",[0.7490, 0.2510, 0.7490])
l1 = ceil(min([mc1;mc2]))-1;
l2 = ceil(max([mc1;mc2]));
disp(l1)
disp(l2)
hold on
plot(l1:l2,l1:l2,"Color",[0 0 0])
plot(l1:l2,c1+(l1:l2)*c2,"Color",[0.7490, 0.2510, 0.7490])
hold off
text(.05,.9,strcat('y = ',num2str(c1),' + ',num2str(c2),' x'),"Units","normalized","FontSize",9);
text(.05,.75,strcat("R^{2} = ",num2str(R)),"Units","normalized","FontSize",9);
xlim([l1 l2])
ylim([l1 l2])
title('Concentrations','FontSize',14)

for i=1:nvars
    nexttile
    h1 = histogram(1e3 * exp(sample1.points(nRxns+sel_conc(i),:)),40,'DisplayStyle','stairs',"EdgeColor",[0 0.4470 0.7410],"Normalization","pdf");                      % mM (mmol)
    hold on
    histogram(1e3 * exp(sample1.points(nRxns+sel_conc(i),:)),40,"EdgeColor","none","FaceColor",[0 0.4470 0.7410],"Normalization","pdf")
    h2 = histogram(1e3 * exp(sample2.points(nRxns+sel_conc(i),:)),40,'DisplayStyle','stairs',"EdgeColor",[0.8500 0.3250 0.0980],"Normalization","pdf");
    histogram(1e3 * exp(sample2.points(nRxns+sel_conc(i),:)),40,"EdgeColor","none","FaceColor",[0.8500 0.3250 0.0980],"Normalization","pdf")
    hold off
    ylim([0 1.1*max([h1.Values h2.Values])])
    caption = createCaption(metNames{sel_conc(i)});
    title(caption,'FontSize',10)
    p = ranksum(s11(i,:),s22(i,:));
    text(.8,.9,strcat('p = ',num2str(p)),"Units","normalized");
    if i==3
        xlabel("mmol/L","FontSize",12)
    end
end

% dGr
nexttile
c1 = modD.Coefficients.Estimate(1); c2 = modD.Coefficients.Estimate(2); R = modD.Rsquared.Ordinary;
plot(md1,md2,'x',"Color",[0.7490, 0.2510, 0.7490])
l1 = 1.1*min([md1;md2]);
l2 = 1.1*max([md1;md2]);
hold on
plot(l1:l2,l1:l2,"Color",[0 0 0])
plot(l1:l2,c1+(l1:l2)*c2,"Color",[0.7490, 0.2510, 0.7490])
hold off
text(.05,.9,strcat('y = ',num2str(c1),' + ',num2str(c2),' x'),"Units","normalized","FontSize",9);
text(.05,.75,strcat("R^{2} = ",num2str(R)),"Units","normalized","FontSize",9);
xlim([l1 l2])
ylim([l1 l2])
title("\Delta_r G'",'FontSize',14)

for i=1:nvars
    nexttile
    h1 = histogram(dGr1(sel_dGr(i),:),40,'DisplayStyle','stairs',"EdgeColor",[0 0.4470 0.7410],"Normalization","pdf");
    hold on
    histogram(dGr1(sel_dGr(i),:),40,"EdgeColor","none","FaceColor",[0 0.4470 0.7410],"Normalization","pdf")
    h2 = histogram(dGr2(sel_dGr(i),:),40,'DisplayStyle','stairs',"EdgeColor",[0.8500 0.3250 0.0980],"Normalization","pdf");
    histogram(dGr2(sel_dGr(i),:),40,"EdgeColor","none","FaceColor",[0.8500 0.3250 0.0980],"Normalization","pdf")
    hold off
    ylim([0 1.1*max([h1.Values h2.Values])])
    caption = correct_TEXT(bigg_dGr{sel_dGr(i)});
    title(caption,'FontSize',10)
    p = ranksum(dg1(i,:),dg2(i,:));
    text(.05,.9,strcat('p = ',num2str(p)),"Units","normalized");
    if i==3
        xlabel("kJ/mol","FontSize",12)
    end
end

t.Padding = 'compact';
t.TileSpacing = 'compact';
lg = legend('','Uniform','','Non Uniform','Orientation','horizontal');
lg.Layout.Tile = 'South';

end
%% Figure 7

% Load data
load('samples\tests\bigUni.mat')
uni = bigUni; clear bigUni
load('samples\tests\bigNouni.mat')
nouni = bigNouni; clear bigNouni

% Calculate correlation (A)
corruni = corr(uni.points(1:117,:)'); % 117 real variables (no slack, no rhs)
corrnou = corr(nouni.points(1:117,:)');
corr_diff = abs(corruni-corrnou)>0.3;
corrnou = tril(corrnou,-1);
corruni = triu(corruni,1);
corr_global = corruni+corrnou;

% Choose variables (B)
variables = find(ismember(uni.rxns,{'c741[c]' 'c1156[c]' 'c1095[c]' 'dGft_c13[c]' 'dGft_c16[c]'}));
xNames = {"5-10-MethyleneTHF log_{10}(mM)" "5-10-MethenylTHF log_{10}(mM)" "10-FormylTHF log_{10}(mM)" "\Delta_{f} G' ATP (kJ/mol)" "\Delta_{f} G' ADP (kJ/mol)"};
yNames = {"5-10-MethyleneTHF" "5-10-MethenylTHF" "10-FormylTHF" "\Delta_{f} G' ATP" "\Delta_{f} G' ADP"};

nVars = numel(variables);
x1 = uni.points(variables,:)';
x2 = nouni.points(variables,:)';
x1(:,1:3) = log10(1e3*exp(x1(:,1:3)));
x2(:,1:3) = log10(1e3*exp(x2(:,1:3)));

% Plot A
figure(1)
imagesc(corr_global)
xticks([])
yticks([])
title("Uniform","FontWeight","normal","FontSize",14)
ylabel("Non uniform","FontSize",14)
colormap(myCM)

hold on
[cr,cc] = find(corr_diff);
plot(cr,cc,"Marker","square","Color",[1 0.7490 0],"LineStyle","none")
hold off

% figure(3)
% imagesc(corr_global)
% xticks([])
% yticks([])
% title("Uniform","FontWeight","normal")
% ylabel("Non uniform")
% title("A","FontSize",14,"FontWeight","bold")

%% Plot B
figure(2)
t = tiledlayout(nVars,nVars,"TileSpacing","none","Padding","tight");
row = 1;
col = 1;
for i=1:nVars^2
    % Revisar si es histograma o dispersión
    % Histograma:
    % Copiar histograma de plotComparison, incluir x1 y x2 altiro
    % Dispersión:
    % Si estamos del lado derecho de la diagonal dscatter de x1,
    % si estamos del lado izquierdo, dscatter de x2.
    nexttile
    if row==col
        % Histogram
        histogram(x1(:,col),20,'DisplayStyle','stairs',"EdgeColor",[0 0.4470 0.7410])
        hold on
        histogram(x1(:,col),20,"EdgeColor","none","FaceColor",[0 0.4470 0.7410])
        histogram(x2(:,col),20,'DisplayStyle','stairs',"EdgeColor",[0.8500 0.3250 0.0980])
        histogram(x2(:,col),20,"EdgeColor","none","FaceColor",[0.8500 0.3250 0.0980])
        hold off
    elseif col > row
        % Right (x1)
        % dscatter(x1(:,col),x1(:,row),'marker','o','filled',false,'MSIZE',4)
        dscatter(x1(:,col),x1(:,row),'PLOTTYPE','contour','smoothing',1)
    else
        % Left (x2)
        % dscatter(x2(:,col),x2(:,row),'marker','o','filled',false,'MSIZE',4)
        dscatter(x2(:,col),x2(:,row),'PLOTTYPE','contour','smoothing',5)
    end
    if col==1
        ylabel(yNames(row),"FontSize",12,"Rotation",0,'HorizontalAlignment','right','VerticalAlignment','middle')
    end
    if row==nVars
        xlabel(xNames(col),"FontSize",12)
    end
    xticks([])
    yticks([])
    % Counters
    if col==nVars
        col = 1;
        row = row+1;
    else
        col = col+1;
    end
end

% title(t,"B","FontSize",14,"FontWeight","bold")

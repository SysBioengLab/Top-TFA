% plotMatrix_final
% Own implementation of improved plotMatrix

load('samples\tests\bigUni.mat')
uni = bigUni; clear bigUni
load('samples\tests\bigNouni.mat')
nouni = bigNouni; clear bigNouni

% Acetate kinase
variables = find(ismember(uni.rxns,{'dGft_c13[c]' 'dGft_c16[c]' 'c741[c]' 'c1095[c]' 'c1156[c]'}));
varNames = {"\Delta_{f} G' ATP" "\Delta_{f} G' ADP" "5-10-MethyleneTHF" "10-FormylTHF" "5-10-MethenylTHF"};
% variables = 81:89;
% variables = [4 5 6 23 24 33 34]; % selected fluxes
% variables = find(ismember(uni.rxns,{'c638[c]' 'c1164[c]' 'c1466[c]'...
    % 'c1311][c]' 'c741[c]' 'c1156[c]' 'c1095[c]' 'c1298[c]' 'c1298[e]' 'c1311[e]'})); % selected concentrations
% variables = find(ismember(uni.rxns,{'r1972' 'r1988' 'r2554' 'r2734' 'c741[c]' 'c1156[c]' 'c1095[c]' 'dGft_c741[c]' 'dGft_c1156[c]' 'dGft_c1095[c]'}));
nVars = numel(variables);
x1 = uni.points(variables,:)';
x2 = nouni.points(variables,:)';

% Plot
f1 = figure(1);
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
        dscatter(x1(:,col),x1(:,row),'PLOTTYPE','contour')
    else
        % Left (x2)
        % dscatter(x2(:,col),x2(:,row),'marker','o','filled',false,'MSIZE',4)
        dscatter(x2(:,col),x2(:,row),'PLOTTYPE','contour')
    end
    if col==1
        ylabel(varNames(row),"FontSize",13)
    end
    if row==nVars
        xlabel(varNames(col),"FontSize",14)
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

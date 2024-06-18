function plotAllData(model,sampleName,nTop,join,dir)
% Genera los graficos de flujo de reacciones, concentración de metabolitos
% y valor de dGr. Agrega el campo dGr a sample.

if ~join
    for j=1:nTop
        disp(["Generating plots for topology ",num2str(j)])
        load(strcat("samples\sample_",sampleName,"_top_",num2str(j),".mat"),"sample")
        % Inicialización de parametros
        nRxns = size(model.S,2);
        nMets = size(model.S,1);
        ndGr  = size(sample.dGrRxns,1);
        
        % Grafico de flujos
        f1 = figure(1);
        t = tiledlayout(5,8);
        for i=1:nRxns
            nexttile
            histogram(sample.points(i,:),20)
            caption = createCaption(model.rxnNames{i});
            title(caption,'FontSize',8)
        end
        t.Padding = 'compact';
        t.TileSpacing = 'compact';
        title(t,'Flujos')
        
        % Grafico de concentraciones
        f2 = figure(2);
        t = tiledlayout(6,7);
        for i=1:nMets
            nexttile
            histogram(1e3 * exp(sample.points(nRxns+i,:)),20)              % mM/L (mmol por litro)
            caption = createCaption(model.metNames{i});
            title(caption,'FontSize',8)
        end
        t.Padding = 'compact';
        t.TileSpacing = 'compact';
        title(t,'Concentraciones (mmol)')

        % Grafico de dGr
        dGr = sample.dGrMatrix * [sample.points(nRxns+1:nRxns+2*nMets,:); ones(1,sample.numSamples)];
        
        h2oIndex = find(contains(sample.dGrRxns,'r2030'));
        f3 = figure(3);
        t = tiledlayout(3,10);
        for i=1:ndGr
            if ~isequal(i,h2oIndex)
                nexttile
                histogram(dGr(i,:),20)
                caption = createCaption(sample.dGrRxns{i});
                title(caption,'FontSize',8)
            end
        end
        t.Padding = 'compact';
        t.TileSpacing = 'compact';
        title(t,'dGr')

        % Guardar archivos
        saveas(f1,strcat(dir,sampleName,"\fluxes_top_",num2str(j)))
        saveas(f2,strcat(dir,sampleName,"\concentrations_top_",num2str(j)))
        saveas(f3,strcat(dir,sampleName,"\dGr_top_",num2str(j)))
    end
else
    % HARDCODE
    disp("Generating joint plots")
    load(strcat("samples\sample_",sampleName,"_top_1.mat"),"sample")
    sample1 = sample; clear sample;
    load(strcat("samples\sample_",sampleName,"_top_2.mat"),"sample")
    sample2 = sample; clear sample;
    load(strcat("samples\sample_",sampleName,"_top_3.mat"),"sample")
    sample3 = sample; clear sample;
    load(strcat("samples\sample_",sampleName,"_top_4.mat"),"sample")
    sample4 = sample; clear sample;
    
    % Inicialización de parametros
    nRxns = size(model.S,2);
    nMets = size(model.S,1);
    ndGr  = size(sample1.dGrRxns,1);
    
    % Grafico de flujos
    f1 = figure(1);
    t = tiledlayout(5,8);
    for i=1:nRxns
        nexttile
        histogram(sample1.points(i,:),20)
        hold on
        histogram(sample2.points(i,:),20)
        histogram(sample3.points(i,:),20)
        histogram(sample4.points(i,:),20)
        hold off
        caption = createCaption(model.rxnNames{i});
        title(caption,'FontSize',8)
    end
    t.Padding = 'compact';
    t.TileSpacing = 'compact';
    title(t,'Flujos')

    % Grafico de concentraciones
    f2 = figure(2);
    t = tiledlayout(6,7);
    for i=1:nMets
        nexttile
        histogram(1e3 * exp(sample1.points(nRxns+i,:)),20)                      % mM (mmol)
        hold on
        histogram(1e3 * exp(sample2.points(nRxns+i,:)),20)
        histogram(1e3 * exp(sample3.points(nRxns+i,:)),20)
        histogram(1e3 * exp(sample4.points(nRxns+i,:)),20)
        hold off
        caption = createCaption(model.metNames{i});
        title(caption,'FontSize',8)
    end
    t.Padding = 'compact';
    t.TileSpacing = 'compact';
    title(t,'Concentraciones (mmol)')

    % Grafico de dGr
    dGr1 = sample1.dGrMatrix * [sample1.points(nRxns+1:nRxns+2*nMets,:); ones(1,sample1.numSamples)];
    dGr2 = sample2.dGrMatrix * [sample2.points(nRxns+1:nRxns+2*nMets,:); ones(1,sample2.numSamples)];
    dGr3 = sample3.dGrMatrix * [sample3.points(nRxns+1:nRxns+2*nMets,:); ones(1,sample3.numSamples)];
    dGr4 = sample4.dGrMatrix * [sample4.points(nRxns+1:nRxns+2*nMets,:); ones(1,sample4.numSamples)];
    dGrRxns = model.rxnNames(~findExRxns(model));

    f3 = figure(3);
    t = tiledlayout(3,10);
    for i=1:ndGr
        nexttile
        histogram(dGr1(i,:),20)
        hold on
        histogram(dGr2(i,:),20)
        histogram(dGr3(i,:),20)
        histogram(dGr4(i,:),20)
        hold off
        caption = createCaption(dGrRxns{i});
        title(caption,'FontSize',8)
    end
    t.Padding = 'compact';
    t.TileSpacing = 'compact';
    title(t,'dGr')
    
    % Guardar archivos
    saveas(f1,strcat(dir,sampleName,"\fluxes_joined"))
    saveas(f2,strcat(dir,sampleName,"\concentrations_joined"))
    saveas(f3,strcat(dir,sampleName,"\dGr_joined"))
end
close all
end
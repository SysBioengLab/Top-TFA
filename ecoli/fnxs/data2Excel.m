function data2Excel(ecM,version)
% Saves data and generate tables 
% It is assumed that TMFA and Top-TMFA have the same feasible topologies.

for N=ecM
    filename = strcat('Results\Excels\final_ec',num2str(N),".xlsx");
    for vX=version
        my_file = strcat('files\ec',num2str(N),"\DATAv",num2str(vX),".mat");
        load(my_file)
        sheet = 'Topologies';
        revs = bigModel.thermo.rxns(logical(bigModel.thermo.rev))';
        if vX==1
            writecell(revs,filename,'Sheet',sheet,'Range','B4')
            writematrix(tfm,filename,'Sheet',sheet,'Range','B5')
            writematrix(feasibleTMFA,filename,'Sheet',sheet,'Range',strcat(hCelladd('B',size(revs,2)),'5'))
            writematrix(feasibleFcTMFA,filename,'Sheet',sheet,'Range',strcat(hCelladd('B',size(revs,2)+1),'5'))
            last = hCelladd('B',size(revs,2)+1); % Last letter used;
        else
            ini = hCelladd(last,3);
            writecell(revs,filename,'Sheet',sheet,'Range',strcat(ini,'4'))
            writematrix(tfm,filename,'Sheet',sheet,'Range',strcat(ini,'5'))
            writematrix(feasibleTMFA,filename,'Sheet',sheet,'Range',strcat(hCelladd(ini,size(revs,2)),'5'))
            writematrix(feasibleFcTMFA,filename,'Sheet',sheet,'Range',strcat(hCelladd(ini,size(revs,2)+1),'5'))
            last = hCelladd(ini,size(revs,2)+1); % Last letter used;
        end
        % Headers construction
        hFVA = 1:size(tfm,1);
        hTMFA = hFVA;
        hfcTMFA = hTMFA;
        hTMFA(~feasibleTMFA) = [];
        hfcTMFA(~feasibleFcTMFA) = [];
        rH = zeros(1,2*size(hFVA,2));
        rHT = zeros(1,2*size(hTMFA,2));
        rHf = zeros(1,2*size(hfcTMFA,2));
        for i=1:size(hFVA,2)
            rH([2*i-1 2*i]) = hFVA(i);
        end
        for i=1:size(hTMFA,2)
            rHT([2*i-1 2*i]) = hTMFA(i);
        end
        for i=1:size(hfcTMFA,2)
            rHf([2*i-1 2*i]) = hfcTMFA(i);
        end

        sheet = strcat('Data TR v',num2str(vX));
        nboundsFVA  = reOrderBounds(boundsFVA,true(size(feasibleTMFA)));
        c_bounds    = repmat([bigModel.thermo.lbc bigModel.thermo.ubc],1,sum(true(size(feasibleTMFA))));
        dgf_bounds  = repmat([bigModel.thermo.dGftLB bigModel.thermo.dGftUB],1,sum(true(size(feasibleTMFA))));
        nnboundsFVA = [nboundsFVA;c_bounds;dgf_bounds];
        metsIdx = num2str(5+size(bigModel.thermo.rxns,1));
        dGfIdx  = num2str(5+size(bigModel.thermo.rxns,1)+size(bigModel.thermo.dGftLB,1));
        dGf1 = cell(size(bigModel.thermo.mets,1),1);
        dGf1(:) = {'dGft_'};
        dGfNames = strcat(dGf1,bigModel.thermo.mets);
        for i=1:sum(feasibleTMFA); c_bounds(:,[2*i-1 2*i]) = [bigModel.thermo.lbc bigModel.thermo.ubc]; end
        writematrix(rH,filename,'Sheet',sheet,'Range','C4')
        writecell(bigModel.thermo.rxns,filename,'Sheet',sheet,'Range','B5')
        writecell(bigModel.thermo.mets,filename,'Sheet',sheet,'Range',strcat('B',metsIdx))
        writecell(dGfNames,filename,'Sheet',sheet,'Range',strcat('B',dGfIdx))
        writematrix(nnboundsFVA,filename,'Sheet',sheet,'Range','C5')

        sheet = strcat('TMFA Data v',num2str(vX));
        boundsTMFA = reOrderBounds(boundsTMFA,feasibleTMFA);
        writematrix(rHT,filename,'Sheet',sheet,'Range','C4')
        writecell(bigModel.thermo.rxns,filename,'Sheet',sheet,'Range','B5')
        writecell(bigModel.thermo.mets,filename,'Sheet',sheet,'Range',strcat('B',metsIdx))
        writecell(dGfNames,filename,'Sheet',sheet,'Range',strcat('B',dGfIdx))
        writematrix(boundsTMFA,filename,'Sheet',sheet,'Range','C5')

        sheet = strcat('fcTMFA Data v',num2str(vX));
        boundsFcTMFA = reOrderBounds(boundsFcTMFA,feasibleFcTMFA);
        writematrix(rHf,filename,'Sheet',sheet,'Range','C4')
        writecell(bigModel.thermo.rxns,filename,'Sheet',sheet,'Range','B5')
        writecell(bigModel.thermo.mets,filename,'Sheet',sheet,'Range',strcat('B',metsIdx))
        writecell(dGfNames,filename,'Sheet',sheet,'Range',strcat('B',dGfIdx))
        writematrix(boundsFcTMFA,filename,'Sheet',sheet,'Range','C5')

        % Calculations
        tol = 1e-10;
        sheet = strcat('RangeDiff TMFAvsFVA v',num2str(vX));
        rDF1 = rangeDifference(boundsTMFA,nnboundsFVA,feasibleTMFA,true(size(feasibleTMFA)),bigModel.thermo,tol);
        writematrix(hTMFA,filename,'Sheet',sheet,'Range','C4')
        writecell(bigModel.thermo.rxns,filename,'Sheet',sheet,'Range','B5')
        writecell(bigModel.thermo.mets,filename,'Sheet',sheet,'Range',strcat('B',metsIdx))
        writecell(dGfNames,filename,'Sheet',sheet,'Range',strcat('B',dGfIdx))
        writematrix(rDF1,filename,'Sheet',sheet,'Range','C5')

        sheet = strcat('RangeDiff fcTMFAvsFVA v',num2str(vX));
        rDF2 = rangeDifference(boundsFcTMFA,nnboundsFVA,feasibleFcTMFA,true(size(feasibleTMFA)),bigModel.thermo,tol);
        writematrix(hfcTMFA,filename,'Sheet',sheet,'Range','C4')
        writecell(bigModel.thermo.rxns,filename,'Sheet',sheet,'Range','B5')
        writecell(bigModel.thermo.mets,filename,'Sheet',sheet,'Range',strcat('B',metsIdx))
        writecell(dGfNames,filename,'Sheet',sheet,'Range',strcat('B',dGfIdx))
        writematrix(rDF2,filename,'Sheet',sheet,'Range','C5')

        sheet = strcat('RangeDiff fcTMFAvsTMFA v',num2str(vX));
        rDF3 = rangeDifference(boundsFcTMFA,boundsTMFA,feasibleFcTMFA,feasibleTMFA,bigModel.thermo,tol);
        writematrix(hfcTMFA,filename,'Sheet',sheet,'Range','C4')
        writecell(bigModel.thermo.rxns,filename,'Sheet',sheet,'Range','B5')
        writecell(bigModel.thermo.mets,filename,'Sheet',sheet,'Range',strcat('B',metsIdx))
        writecell(dGfNames,filename,'Sheet',sheet,'Range',strcat('B',dGfIdx))
        writematrix(rDF3,filename,'Sheet',sheet,'Range','C5')

    end
end
end
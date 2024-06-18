function closa_data2Excel()
% Saves data and generate tables 

filename = strcat('results\Excels\final_closa.xlsx');
for vX=[1 2]
    my_file = strcat('results\DATA_v',num2str(vX),'.mat');
    load(my_file,"thermoRCA","tfm","boundsFVA","feasible","boundsTcTMFA")
    sheet = 'Topologies';
    revs = thermoRCA.rxns(logical(thermoRCA.rev))';
    if vX==1
        writecell(revs,filename,'Sheet',sheet,'Range','B4')
        writematrix(tfm,filename,'Sheet',sheet,'Range','B5')
        writematrix(feasible,filename,'Sheet',sheet,'Range',strcat(hCelladd('B',size(revs,2)),'5'))
        last = hCelladd('B',size(revs,2)+1); % Last letter used;
    else
        ini = hCelladd(last,3);
        writecell(revs,filename,'Sheet',sheet,'Range',strcat(ini,'4'))
        writematrix(tfm,filename,'Sheet',sheet,'Range',strcat(ini,'5'))
        writematrix(feasible,filename,'Sheet',sheet,'Range',strcat(hCelladd(ini,size(revs,2)),'5'))
        last = hCelladd(ini,size(revs,2)+1); % Last letter used;
    end
    % Headers construction
    hFVA = 1:size(tfm,1);
    rH = zeros(1,2*size(hFVA,2));
    for i=1:size(hFVA,2)
        rH([2*i-1 2*i]) = hFVA(i);
    end
    sheet = strcat('Data TR v',num2str(vX));
    nboundsFVA  = reshape(boundsFVA,size(boundsFVA,1),2*size(boundsFVA,3));
    c_bounds    = repmat([thermoRCA.lbc thermoRCA.ubc],1,sum(true(size(feasible))));
    dgf_bounds  = repmat([thermoRCA.dGftLB thermoRCA.dGftUB],1,sum(true(size(feasible))));
    nnboundsFVA = [nboundsFVA;c_bounds;dgf_bounds];
    metsIdx = num2str(5+size(thermoRCA.rxns,1));
    dGfIdx  = num2str(5+size(thermoRCA.rxns,1)+size(thermoRCA.dGftLB,1));
    dGf1 = cell(size(thermoRCA.mets,1),1);
    dGf1(:) = {'dGft_'};
    dGfNames = strcat(dGf1,thermoRCA.mets);
    % for i=1:sum(feasibleTMFA); c_bounds(:,[2*i-1 2*i]) = [bigModel.thermo.lbc bigModel.thermo.ubc]; end
    writematrix(rH,filename,'Sheet',sheet,'Range','C4')
    writecell(thermoRCA.rxns,filename,'Sheet',sheet,'Range','B5')
    writecell(thermoRCA.mets,filename,'Sheet',sheet,'Range',strcat('B',metsIdx))
    writecell(dGfNames,filename,'Sheet',sheet,'Range',strcat('B',dGfIdx))
    writematrix(nnboundsFVA,filename,'Sheet',sheet,'Range','C5')

    sheet = strcat('tcTMFA Data v',num2str(vX));
    boundsTcTMFA = reshape(boundsTcTMFA,size(boundsTcTMFA,1),2*size(boundsTcTMFA,3));
    writematrix(rH,filename,'Sheet',sheet,'Range','C4')
    writecell(thermoRCA.rxns,filename,'Sheet',sheet,'Range','B5')
    writecell(thermoRCA.mets,filename,'Sheet',sheet,'Range',strcat('B',metsIdx))
    writecell(dGfNames,filename,'Sheet',sheet,'Range',strcat('B',dGfIdx))
    writematrix(boundsTcTMFA,filename,'Sheet',sheet,'Range','C5')

    % Calculations
    tol = 1e-10;
    sheet = strcat('RangeDiff tcTMFAvsFVA v',num2str(vX));
    rDF = rangeDifference(boundsTcTMFA,nnboundsFVA,feasible,true(size(feasible)),thermoRCA,tol);
    writematrix(hFVA,filename,'Sheet',sheet,'Range','C4')
    writecell(thermoRCA.rxns,filename,'Sheet',sheet,'Range','B5')
    writecell(thermoRCA.mets,filename,'Sheet',sheet,'Range',strcat('B',metsIdx))
    writecell(dGfNames,filename,'Sheet',sheet,'Range',strcat('B',dGfIdx))
    writematrix(rDF,filename,'Sheet',sheet,'Range','C5')

end
end
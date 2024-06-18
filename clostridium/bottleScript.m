% BottleScript

bottleOutputs1 = closa_bottleneckAnalysis("results\DATA_v1.mat");
bottleOutputs2 = closa_bottleneckAnalysis("results\DATA_v2.mat");

% Save on excel
filename = strcat('results\excels\final_closa.xlsx');
sheet = 'Bottleneck analysis';
writecell(bottleOutputs1.dGrRxns,filename,'Sheet',sheet,'Range','B4')
writematrix(1:4,filename,'Sheet',sheet,'Range','C3')
writematrix(bottleOutputs1.clasification,filename,'Sheet',sheet,'Range','C4')

writecell(bottleOutputs2.dGrRxns,filename,'Sheet',sheet,'Range','H4')
writematrix(1:4,filename,'Sheet',sheet,'Range','I3')
writematrix(bottleOutputs2.clasification,filename,'Sheet',sheet,'Range','I4')

sheet = 'Thermodynamic sensitivity';
writecell(bottleOutputs1.dGrRxns,filename,'Sheet',sheet,'Range','B4')
writematrix(1:4,filename,'Sheet',sheet,'Range','C3')
writematrix(bottleOutputs1.Tsen,filename,'Sheet',sheet,'Range','C4')

writecell(bottleOutputs2.dGrRxns,filename,'Sheet',sheet,'Range','H4')
writematrix(1:4,filename,'Sheet',sheet,'Range','I3')
writematrix(bottleOutputs2.Tsen,filename,'Sheet',sheet,'Range','I4')


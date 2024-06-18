% Generates dGr graph

% Load data
load("..\DGR_cases.mat","case_2","case_3","dGrRxns")
load("..\files\henry2007_dGr.mat","M")

% Remove H2Ot
wtr_indx = find(contains(dGrRxns,'H2Ot'));
case_2(wtr_indx,:) = [];
case_3(wtr_indx,:) = [];
M(wtr_indx,:) = [];

% Reorder
orderM = [median(case_2,2) case_2 case_3 M];
orderM = sortrows(orderM);
c2 = orderM(:,2:3)';
c3 = orderM(:,4:5)';
henryM = orderM(:,6:7)';

figure(1)
boxchart(c2,'BoxWidth',1)
hold on
boxchart(c3,'BoxWidth',1)
boxchart(henryM,'BoxWidth',1)
% yl = yline(4.3478,'-.k',{'Sensitive range'},Fontsize=12);
% yl2 = yline(-4.3478,'-.k');
hold off
ylim([-230 60]);
% yl.LabelHorizontalAlignment = 'left';
legend({'Case 2' 'Case 3' 'Henry2007'},'Location','best','FontSize',16) % '\pm g*_{r}'
xticks([])
xlabel("Reactions")
ylabel("\Delta_{r} G' (kJ/mol)")
ax = gca;
ax.FontSize = 16;
% boxchart(case_2','BoxWidth',1)
% hold on
% boxchart(case_3','BoxWidth',1)
% boxchart(M','BoxWidth',1)
% hold off

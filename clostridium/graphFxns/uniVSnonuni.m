% Plot figures.

% Load data
load('../samples/bigUni.mat')
uni = sample; clear sample
load('../samples/bigNouni.mat')
nouni = sample; clear sample

figure()
t = tiledlayout(1,2,"TileSpacing","compact","Padding","compact");

% A: PSRF
x1 = uni.R';
x2 = nouni.R';

group = [ones(numel(x1),1);2*ones(numel(x2),1)];
positions = [1,2];

ax(1) = nexttile(1);
boxplot([x1;x2],group,'positions', positions,'symbol','');
bar(positions,[mean(x1),mean(x2)],.15,'FaceColor',[0,0.45,0.74],'BarWidth',3);
hold on

p1 = prctile(x1,[2.5,97.5]);
p2 = prctile(x2,[2.5,97.5]);

errorbar(positions,[mean(x1),mean(x2)],[p1(1),p2(1)]-[mean(x1),mean(x2)],...
    [p1(2),p2(2)]-[mean(x1),mean(x2)],'.k')

set(gca,'xtick',positions)
set(gca,'xticklabel',{'Uniform','Non-Uniform'})

axis([.5,2.5,0,1.5])
title('A');
ax(1).TitleHorizontalAlignment = 'Left';
ylabel('Potential Scale Redution Factor (psrf)')

% B: Time per Neff
x1 = uni.samplingTime./uni.Neff';
x2 = nouni.samplingTime./nouni.Neff';

group = [ones(numel(x1),1);2*ones(numel(x2),1)];
positions = [1 2];

ax(2) = nexttile(2);
boxplot(([x1;x2]),group,'positions', positions,'symbol','');

set(gca,'YScale','log')
set(gca,'xtick',positions)
set(gca,'xticklabel',{'Uniform','Non-Uniform'})

color = [1,0.41,0.16; 1,0.41,0.16];
% color = [1,0.41,0.16;...
%          0,0.45,0.74;...
%          1,0.41,0.16;...
%          0,0.45,0.74;...
%          1,0.41,0.16;...
%          0,0.45,0.74;...
%          1,0.41,0.16;...
%          0,0.45,0.74];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j,:),'FaceAlpha',.9);
end
c = get(gca,'Children');

axis([.5,2.5,10^-3,10^3])
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');

title('B');
ax(2).TitleHorizontalAlignment = 'Left';
ylabel('Time per effective sample (s)')

title(t,"Comparison between uniform and non-uniform sampling")

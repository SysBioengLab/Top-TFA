% Corrplot script

load('samples\sample_Uniform_ADSB_top_1.mat')
uni = sample; clear sample
load('samples\sample_Non_uniform_ADSB_top_1.mat')
nuni= sample; clear sample
uniConcPoints = uni.points(38:78,:)';
mi_dev = std(uniConcPoints);
uniConcPoints(:,mi_dev<0.001) = [];
nuniConcPoints = nuni.points(38:78,:)';
mi_dev = std(nuniConcPoints);
nuniConcPoints(:,mi_dev<0.001,:) = [];

% f1 = figure(1);
% corrplot(uni.points(1:18,:)')
% f2 = figure(2);
% corrplot(nuni.points(1:18,:)')
% f3 = figure(3);
% corrplot(uni.points(19:37,:)')
% f4 = figure(4);
% corrplot(nuni.points(19:37,:)')
f5 = figure(5);
corrplot(uniConcPoints(:,1:18))
f6 = figure(6);
corrplot(nuniConcPoints(:,1:18))
f7 = figure(7);
corrplot(uniConcPoints(:,19:36))
f8 = figure(8);
corrplot(nuniConcPoints(:,19:36))


% saveas(f1,"figures\corrplots\fluxes_1to18_uni",'png')
% saveas(f2,"figures\corrplots\fluxes_1to18_nuni",'png')
% saveas(f3,"figures\corrplots\fluxes_19to37_uni",'png')
% saveas(f4,"figures\corrplots\fluxes_19to37_nuni",'png')
% saveas(f5,"figures\corrplots\conc_1to18_uni",'png')
% saveas(f6,"figures\corrplots\conc_1to18_nuni",'png')
% saveas(f7,"figures\corrplots\conc_19to37_uni",'png')
% saveas(f8,"figures\corrplots\conc_19to37_nuni",'png')

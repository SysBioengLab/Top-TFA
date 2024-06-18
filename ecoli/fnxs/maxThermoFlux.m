% Calculate flux by thermodynamics

% Load data
load("files\ec3\DATAv2.mat","bigModel","boundsFcTMFA","feasibleFcTMFA", ...
    "dGrboundsFcTMFA","boundsFVA","fluxSolutionsFcTMFA")

% Define variables
tModel = bigModel.thermo;
exchange = findExRxns(tModel);
dGrRxns = tModel.rxns(~exchange);
nTop = sum(feasibleFcTMFA);
nRxns = size(tModel.rxns,1);
ndGrRxns = size(dGrboundsFcTMFA,1);
maxTflux = zeros(ndGrRxns,nTop);
VmaxMatrix = zeros(ndGrRxns,nTop);
dGrValues = zeros(ndGrRxns,nTop);
boundsFVA(exchange,:,:) = [];
FVAMatrix = zeros(ndGrRxns,nTop);
feasibleFcTMFA(~feasibleFcTMFA) = [];
ban = sort([2*find(exchange)' (2*find(exchange)-1)']);
fluxSolutionsFcTMFA(:,ban,:) = [];
fluxSolutionsFcTMFA(:,:,~feasibleFcTMFA) = [];
dGrMatrix = createdGrMatrix(tModel);

% Main loop (per topology)
for i=find(feasibleFcTMFA')
    % Vmax calculation
    fluxDir = ones(size(tModel.rxns));
    fluxDir(boundsFcTMFA(1:nRxns,2,i)<=1e-8) = -1;
    Vmax = tModel.vmax;
    Vmax(fluxDir==-1) = tModel.vmin(fluxDir==-1);
    Vmax(exchange) = [];
    fluxDir(exchange) = [];
    Vmax = Vmax ./ 4.3478;
    % dGr calculation
    dGrOpts = zeros(ndGrRxns,2);
    for j=1:ndGrRxns
        tempdGrOpts = dGrMatrix * [fluxSolutionsFcTMFA(:,[2*j-1 2*j],i); ones(1,2)];
        dGrOpts(j,:) = tempdGrOpts(j,:);
    end
    % Flux obtention
    dGr =  dGrOpts(:,2); % Positive flux
    dGr(fluxDir==-1) = dGrOpts(fluxDir==-1,1); % Negative flux
    dGrValues(:,i) = dGr;
    maxTflux(fluxDir==1,i) = -Vmax(fluxDir==1) .* dGr(fluxDir==1);
    maxTflux(fluxDir==-1,i) = Vmax(fluxDir==-1) .* dGr(fluxDir==-1);
    VmaxMatrix(:,i) = Vmax;
    FVAMatrix(:,i) = boundsFVA(:,2,i);
    FVAMatrix(fluxDir==-1,i) = boundsFVA(fluxDir==-1,1,i);
end

% Flux comparison
thermoLimitedVmax = abs(maxTflux) < abs(VmaxMatrix) - 1e-4;
thermoLimitedFVA  = abs(maxTflux) < abs(FVAMatrix)  - 1e-4;

% Save data
save("files\ec3\DATAv2_BottleNecks.mat","thermoLimitedVmax", ...
    "thermoLimitedFVA","dGrValues","VmaxMatrix")

% Generate excel
filename = strcat('Results\Excels\ec3.xlsx');
sheet = 'Bottleneck analysis 1';
writecell(dGrRxns,filename,'Sheet',sheet,'Range','B4')
writematrix(1:12,filename,'Sheet',sheet,'Range','C3')
writematrix(thermoLimitedFVA,filename,'Sheet',sheet,'Range','C4')
sheet = 'Bottleneck analysis 2';
writecell(dGrRxns,filename,'Sheet',sheet,'Range','B4')
writematrix(1:12,filename,'Sheet',sheet,'Range','C3')
writematrix(dGrValues,filename,'Sheet',sheet,'Range','C4')
writematrix(VmaxMatrix,filename,'Sheet',sheet,'Range',strcat(hCelladd('C',14),'4'))


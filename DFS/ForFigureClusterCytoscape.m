clear all;clc;
addpath('fxns/')

clusterID = 15; % number of the cluster to be studied
preFileName=num2str(clusterID);

% Extract clusters
load('clusterAnalysis.mat');

% Identify matrix, reactions, metabolites and reversibilities
M = Clusters_all{clusterID,1};
RevRxns  = Clusters_all{clusterID,2};
IDrxn = Clusters_all{clusterID,3};
metsName = Clusters_all{clusterID,4};
external =  Clusters_all{clusterID,5};

%create list of relationships for cytoscape
model=ForCytoscape(M,IDrxn,metsName);

clearvars -except model M RevRxns IDrxn metsName external clusterID

eval(strcat('save ToCytoscape',num2str(clusterID)))
clear all;clc;
load Clusters_iJO1366
addpath ./fxns

Clusters_all=cell(max(clusters),5);
Clusters_size=zeros(max(clusters),3);
IDexternalC=[];

 for i=1:max(clusters)
    [M,dirxn2,IDrxn,metabolites,external]=Clusters(dirxn,Mc_b2,namesc_b2,clusters,i,order,IDShort);
    Clusters_size(i,:)=[sum(dirxn2) length(dirxn2)-sum(dirxn2)-length(external) length(external)];
    Clusters_all{i,1}=M;
    Clusters_all{i,2}=dirxn2;
    Clusters_all{i,3}=IDrxn;
    Clusters_all{i,4}=metabolites;
    Clusters_all{i,5}=external;
    IDexternalC=[IDexternalC; IDrxn(external)];
 end
 
 [IDexternalC2,frecuency]=FrecuencyString(IDexternalC);
 
 %save clusterAnalysis
 
function [M,dirxn2,IDrxn,metabolites,external]=Clusters(dirxn,Mc_b2,namesc_b2,clusters,ClusterN,order,IDShort)
% clear all;clc;
% load Clusters_iJO1366
% ClusterN=7;

MRev=Mc_b2(:,dirxn==1);
MRev=MRev(:,order);
IDrev=IDShort(find(dirxn==1));
IDrev=IDrev(order);

MIrrev=Mc_b2(:,dirxn==0);
IDIrrev=IDShort(find(dirxn==0));
 
M=MRev(:,clusters==ClusterN);
IDrxn=IDrev(clusters==ClusterN);

if size(M,2)>1
    delete=(find(sum(abs(M'))==0));
else
    delete=find(M==0);
end
M(delete,:)=[]; MIrrev(delete,:)=[];

if size(MIrrev,1)>1
    delete2=find(sum(abs(MIrrev))==0);
else
    delete2=find(MIrrev==0);
end
MIrrev(:,delete2)=[]; IDIrrev(delete2)=[];

dirxn2=[ones(size(M,2),1);zeros(size(MIrrev,2),1)];
M=[M MIrrev]; IDrxn=[IDrxn;IDIrrev];

metabolites=namesc_b2;
metabolites(delete)=[];

if size(M,1)>1
    external=unique([find(sum(M>0)==0) find(sum(M<0)==0)]);
    external=external(find(external>sum(dirxn2)));
else
    external=sum(dirxn2)+1:1:length(M);
end





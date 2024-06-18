%author Veronica Martinez, The University of Queensland, AIBN
%25-Jun-2015
function [model,revRxns]=RemoveDuplicateRxns(premodel,prerevRxns,external)
%Model compression, to lump repeated reactions wit the same or opposite
%direction
addpath ./fxns 
reversibility=prerevRxns;
stoicMat=full(premodel.S'); ID=premodel.rxnNames; lb=premodel.lb; ub=premodel.ub;
NoToLump=external;
SNTL=stoicMat(NoToLump,:);
IDNTL=ID(NoToLump); 
lbNTL=lb(NoToLump); ubNTL=ub(NoToLump); 
reversibilityNTL=reversibility(NoToLump);
stoicMat(NoToLump,:)=[];ID(NoToLump)=[];lb(NoToLump)=[];ub(NoToLump)=[];
reversibility(NoToLump)=[];stoicMat2=sign(stoicMat);
% % remove repeated reactions in the matrix
[~,index] = unique(stoicMat2,'rows','first');
repeatedIndex = setdiff(1:size(stoicMat2,1),index);
stoicMatSmall=stoicMat;
reversibilitySmall=reversibility;
IDSmall=ID;
lbSmall=lb; ubSmall=ub; 
% fix by repeated rxns
ToRemove=[];
for i=1:length(repeatedIndex)
    [~,indx]=ismember(stoicMat,stoicMat(repeatedIndex(i),:),'rows');
    [~,indx2]=ismember(stoicMat2,stoicMat2(repeatedIndex(i),:),'rows');
    indx=find(indx); indx2=find(indx2); 
    indx3=setdiff(indx2,indx);
    NoIndx=[];
    for j=1:length(indx3)
        check=unique(stoicMat(indx(1),:)./stoicMat(indx3(j),:));
        check(isnan(check))=[];
        if length(check)>1
            NoIndx=[NoIndx; j];
        else
            lbSmall(indx3(j))=lb(indx3(j))*check;
            ubSmall(indx3(j))=ub(indx3(j))*check;
        end
    end
    indx3(NoIndx)=[]; indx2=[indx;indx3];
    ToRemove=[ToRemove; indx2(2:end)];
    preID=ID{indx2(1)};
    for j=1:length(indx2)-1
        preID=strcat(preID,'_',ID{indx2(j+1)});
    end
    IDSmall{indx2(1)}=preID;
    reversibilitySmall(indx2(1))=max(reversibility(indx2));
    ubSmall(indx2(1))=max([ub(indx2);ubSmall(indx3)]);
    lbSmall(indx2(1))=min([lb(indx2);lbSmall(indx3)]);
end
ToRemove=unique(ToRemove);
stoicMatSmall(ToRemove,:)=[];
reversibilitySmall(ToRemove)=[];
IDSmall(ToRemove)=[];
lbSmall(ToRemove)=[];ubSmall(ToRemove)=[];

% remove repeated reactions in the opposite direction
stoicMat2=stoicMatSmall;
stoicMat3=sign(stoicMat2);
repeatedIndex2=[];
for i=1:size(stoicMatSmall,1)-1
    if isempty(find(i==repeatedIndex2, 1))
        stoicMat2(1,:)=-stoicMat2(1,:);
        stoicMat3(1,:)=-stoicMat3(1,:);
        [~,index] = unique(stoicMat3,'rows','first');
        if length(index)<size(stoicMat3,1)
            [~,indx]=ismember(stoicMat2,stoicMat2(1,:),'rows');
            [~,indx2]=ismember(stoicMat3,stoicMat3(1,:),'rows');
            indx=find(indx)+i-1; indx2=find(indx2)+i-1;
            if ~isempty(indx2)
                indx3=setdiff(indx2,indx);
                NoIndx=[];
                for j=1:length(indx3)
                    check=unique(stoicMat2(indx(1)-i+1,:)./stoicMat2(indx3(j)-i+1,:));
                    check(isnan(check))=[];
                    if length(check)>1
                        NoIndx=[NoIndx; j];
                    else
                        lbSmall(indx3(j))=lbSmall(indx3(j))*check;
                        ubSmall(indx3(j))=ubSmall(indx3(j))*check;
                    end
                end
                indx3(NoIndx)=[]; indx2=[indx;indx3];
                indx2(1)=[];
                if ~isempty(indx2)
                    repeatedIndex2=[repeatedIndex2;indx2];
                    preID=IDSmall{i};
                    for j=1:length(indx2)
                        preID=strcat(preID,'_',IDSmall{indx2(j)});
                    end
                    IDSmall{i}=preID;
                    reversibilitySmall(i)=1;
                    lbSmall(i)=-max(ubSmall(indx2),abs(lbSmall(i)));
                    ubSmall(i)=max(abs(lbSmall(indx2)),ubSmall(i));
                end
            end
        end
    end
    stoicMat2(1,:)=[];
    stoicMat3(1,:)=[];
end
stoicMatSmall(repeatedIndex2,:)=[];
reversibilitySmall(repeatedIndex2)=[];
IDSmall(repeatedIndex2)=[];
lbSmall(repeatedIndex2)=[];
ubSmall(repeatedIndex2)=[];
stoicMatSmall=[stoicMatSmall;SNTL];
lbSmall=[lbSmall;lbNTL];
ubSmall=[ubSmall;ubNTL];
IDSmall=[IDSmall;IDNTL];
reversibilitySmall=[reversibilitySmall;reversibilityNTL];
model=premodel;
model.S=sparse(stoicMatSmall');
model.c=ones(size(model.S,2),1);
model.lb=lbSmall;
model.ub=ubSmall;
model.rxns=IDSmall;
model.rxnNames=IDSmall;
revRxns=reversibilitySmall;

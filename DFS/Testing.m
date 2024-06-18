 preFileName=num2str(clusterID);
for i=1:3-length(preFileName)
    preFileName=strcat('0',preFileName);
end
TotalConfig=0;
Numconfig=zeros(692,1);
for i=1:692
    preFileName2=num2str(i);
    for j=1:4-length(preFileName2)
        preFileName2=strcat('0',preFileName2);
    end
    binaryFile=strcat('TFM_',preFileName,'_',preFileName2,'.bin');
    BinaryOutput = readBinaryCfgs(binaryFile, sum(revRxns));
    
    Numconfig(i)=size(BinaryOutput,1);
    TotalConfig=TotalConfig+Numconfig(i);
    fprintf(1,'File #: %d Number of configs: %d.\n',i,TotalConfig);
end
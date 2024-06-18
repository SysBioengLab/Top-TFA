function output = readBinaryCfgs(binaryFile, numRxns)
fileID = fopen(binaryFile);
elements = fread(fileID, inf, 'uint8');

% Find number of row for clusters with more than 8 rxns
if numRxns>8
    numCharsPerRow = ceil(numRxns/8);
    output=zeros(length(elements)/numCharsPerRow,numCharsPerRow);
    for i=1:length(elements)/numCharsPerRow
        preOutput=elements(numCharsPerRow*(i-1)+1:1:numCharsPerRow*(i));
        output(i,:)=preOutput';
    end
else
    output=elements;
end
fclose(fileID);
end


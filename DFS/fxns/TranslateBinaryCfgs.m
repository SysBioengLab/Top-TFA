function output = TranslateBinaryCfgs(BinaryOutput, numRxns)
numRows=size(BinaryOutput,1);
%Create the output array with zeros and ones
output = zeros(numRows,numRxns);
for i = 1:numRows
    rowPos=1;
    for j = 1:size(BinaryOutput,2)-1
        %Read the number in
        binString = dec2bin(BinaryOutput(i, j), 8);
        assert(length(binString) == 8);
        for k = 1:8
            output(i,rowPos) = str2double(binString(k)) * 2 - 1;
            rowPos = rowPos + 1;
        end        
    end
    %Read only the significant bits from the final number
    finalBits = dec2bin(BinaryOutput(i,end),8);
    extraBits=(8*size(BinaryOutput,2))-numRxns;
    finalBits = finalBits(extraBits+1:end);
    
    %THIS SHOULD IGNORE THE LEADING BITS e.g. ->(00)1101001
    for k = 1:8-extraBits
        output(i,rowPos) = str2double(finalBits(k))  * 2 - 1;
        rowPos = rowPos + 1;
    end
end




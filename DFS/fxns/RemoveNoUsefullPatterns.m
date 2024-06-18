function unfeasiblePattern=RemoveNoUsefullPatterns(FullunfeasiblePattern) 

unfeasiblePattern2=FullunfeasiblePattern;
unfeasiblePattern2(logical(unfeasiblePattern2))=1;
[counts,index]=sort(sum(unfeasiblePattern2,2));
unfeasiblePattern=FullunfeasiblePattern(index,:);
ListToRemove=[];
for i=max(counts==1)+1:size(FullunfeasiblePattern,1)
    LongPattern=unfeasiblePattern(i,:);
    for j=1:i-1
        ShortPattern=unfeasiblePattern(j,:);
        LongPattern(ShortPattern==0) = 0;
        if sum(LongPattern==ShortPattern)==length(LongPattern)
            ListToRemove=[ListToRemove; i];
            break
        end
       LongPattern=unfeasiblePattern(i,:); 
    end
end
unfeasiblePattern(ListToRemove,:)=[];
function withoutempty=DeleteEmpty(withempty)
for i = length(withempty):-1:1
    if isempty(withempty{i})
        withempty(i)=[];
    end
end
withoutempty=withempty;
    
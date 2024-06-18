function newWord = correct_TEXT(word)
% Replaces _ with \_
N = strfind(word,"_");

if isequal(N,[])
    newWord = word;
else
    newWord = strcat(word(1:N-1),'\',word(N:end));
end
end
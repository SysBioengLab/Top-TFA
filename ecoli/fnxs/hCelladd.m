function letter = hCelladd(li,ni)
% Takes the column li and adds ni. Returns the new column. Excel.

if size(li,2) == 1
    prev = 64;
    num = double(li) + ni;
else
    prev = double(li(1));
    num  = double(li(2)) + ni;
end

if num <= double('Z')
    if prev<65
        letter = char(num);
    else
        letter = char([prev num]);
    end
else
    while num > double('Z')
        prev = prev+1;
        num = num-26;
    end
    letter = char([prev num]);
end
end
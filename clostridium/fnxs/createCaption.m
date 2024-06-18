function caption = createCaption(s)

if size(s,2) <= 30
    caption = s;
elseif size(s,2) <= 60
    caption = [convertCharsToStrings(s(1:30)),convertCharsToStrings(s(31:end))];
elseif size(s,2) <= 90
    caption = [convertCharsToStrings(s(1:30)),convertCharsToStrings(s(31:60)),convertCharsToStrings(s(61:end))];
else
    caption = [convertCharsToStrings(s(1:30)),convertCharsToStrings(s(31:60)),convertCharsToStrings(s(61:90))];
end
function [newStyle] = fixDirectionalities(oldStyle)
newStyle = int16(zeros(size(oldStyle, 1), 3));
for i = 1:size(oldStyle, 2)
    indices = oldStyle{i} + 2;
    newStyle(i, indices) = 1;
end
end


function model = generateNH(model)
% Calculates the number of protons of mets in a model.
% Requires the field metFormulas.

NH = zeros(size(model.mets));

for i=1:size(model.metFormulas)
    molecule = model.metFormulas{i};
    n = strfind(molecule,'H') + 1;
    if ~isempty(n)
        if n > size(molecule,2)
        NH(i) = 1;
        else
            if isnan(str2double(molecule(n)))
                NH(i) = 1;
            else
                if n+1 > size(molecule,2)
                    NH(i) = str2double(molecule(n));
                elseif isnan(str2double(molecule(n+1)))
                    NH(i) = str2double(molecule(n));
                else
                    NH(i) = str2double(molecule(n:n+1));
                end
            end
        end
    end
end

model.NH = NH;
end
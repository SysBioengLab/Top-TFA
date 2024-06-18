function [string,frecuency]=FrecuencyString(cell)
cell=DeleteEmpty(cell);
cellu=unique(cell);
frecuency=zeros(length(cellu),1);
for j=1:length(cellu)
    frecuency(j)=length(strmatch(cellu(j),cell,'exact'));
end
string=cellu;

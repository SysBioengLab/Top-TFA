function Side=concatenateSpecies(v,s)
if isempty(v)
    Side={};
end
for i=1:length(v)
   if i==1;
       if s(i)==1
           Side=v(i);
       else
           sn=num2str(s(i));
           Side=strcat('(',sn,')',v(i));
       end
   else
        if s(i)==1
           Side=strcat(Side,'+',v(i));
       else
           sn=num2str(s(i));
           Side=strcat(Side,'+','(',sn,')',v(i));
        end
   end
end
        
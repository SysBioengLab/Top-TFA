function rxnLine = stoicrxnLine2(stoicMat,speciesListName)
rxnLine=cell(size(stoicMat,1),1);
for i = 1:size(stoicMat,1)
    hitReactants = speciesListName(find(stoicMat(i,:)<0));
    hitReactantsStoic = abs(stoicMat(i,stoicMat(i,:)<0));
    hitProducts = speciesListName(stoicMat(i,:)>0);
    hitProductsStoic = abs(stoicMat(i,stoicMat(i,:)>0));
    compr=cell(length(hitReactants),1);
    compp=cell(length(hitProducts),1);
    r=hitReactants;
    p=hitProducts;
    for j=1:length(hitReactants)
        rea=hitReactants{j};
        compr{j}=rea(end-1:end-1);
        r{j}=rea(1:end-3);
    end
    for j=1:length(hitProducts)
        prd=hitProducts{j};
        compp{j}=prd(end-1:end-1);
        p{j}=prd(1:end-3);
    end
    count=0;
    if length(hitReactants)==1
        count=0;
    else
        for j=2:length(hitReactants)
            j1=j-1;
            a=compr{j1};
            b=compr{j};
            count=count+strfind(a,b);
        end
    end
    for j=1:length(hitProducts)
        if j==1 && ~isempty(compr)
            a=compr{1};
            b=compp{j};
            count=count+strfind(a,b);
        elseif j>1;
            j1=j-1;
            a=compp{j1};
            b=compp{j};
            count=count+strfind(a,b);
        end
    end
   if count==length(hitReactants)+length(hitProducts)-1
       if count>0
           if i == size(stoicMat,1)
               LHS = concatenateSpecies(hitReactants,hitReactantsStoic);
               RHS = concatenateSpecies(p,hitProductsStoic);
           else
               LHS = concatenateSpecies(r,hitReactantsStoic);
               if ~isempty(LHS)
                   LHS = strcat('[',compr(1),']',LHS);
                   RHS = concatenateSpecies(p,hitProductsStoic);
               else
                   RHS = concatenateSpecies(hitProducts,hitProductsStoic);
               end
           end
       else
           LHS = concatenateSpecies(hitReactants,hitReactantsStoic);
           RHS = concatenateSpecies(hitProducts,hitProductsStoic);
       end
   else
       LHS = concatenateSpecies(hitReactants,hitReactantsStoic);
       RHS = concatenateSpecies(hitProducts,hitProductsStoic);
   end
   if isempty(LHS)
       LHS{1}=' ';
   end
   if isempty(RHS)
       RHS{1}=' ';
   end
       rxnLine(i)=strcat(LHS,'=',RHS);
       fprintf('%1.0f ',i);
end


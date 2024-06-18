function model=ForCytoscape(M,IDrxn,metsName)
model={};
counter=1;
for i=1:size(M,2)
    substrates=find(M(:,i)<0);
    products=find(M(:,i)>0);
    premodel=[];
    for j=1:length(substrates)
        if j==1
            premodel=[metsName(substrates(j)) num2str(counter) IDrxn(i)];
            counter=counter+1;
        else
            premodel=[premodel;metsName(substrates(j)) num2str(counter) IDrxn(i)];
            counter=counter+1;
        end
    end
    for j=1:length(products)
            premodel=[premodel;IDrxn(i) num2str(counter) metsName(products(j))];
            counter=counter+1;
    end
    if i==1
        model=premodel;
    else
        model=[model;premodel];
    end
end

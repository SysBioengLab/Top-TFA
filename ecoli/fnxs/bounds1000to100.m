function model = bounds1000to100(model)
% Changes limits from 1000 mmol/(gdcw h) to 100 mmol/(gdcw h).

model.ub(model.ub>100)  =  100;
model.lb(model.lb<-100) = -100;
end
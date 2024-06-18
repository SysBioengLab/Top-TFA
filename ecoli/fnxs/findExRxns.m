function exchange = findExRxns(model)

exchange = sum(abs(sign(model.S)),1) == 1;
exchange = exchange';

end
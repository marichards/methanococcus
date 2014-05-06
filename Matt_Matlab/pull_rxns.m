function rxns = pull_rxns(model,met_name)

%Pull out all reactions for a given metabolite NAME and print their formulas

%First find the metabolite
[~,idx]=intersect(model.metNames,met_name);

%Grab the reactions for that metabolite
rxns=model.rxns(find(model.S(idx,:)));

%Print the reactions
printRxnFormula(model,rxns);
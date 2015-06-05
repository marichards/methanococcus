function pullNADNADP(model,solution,workbook)

%Pull out the important biosynthetic and catabolic pathways 

%Do an FVA and FBA
[minFlux,maxFlux] = fluxVariability(model);

%Find all the reactions that have both NAD and NADP using nadRxns
rxns = nadRxns(model);

%Find them in the model
[rxns,idx]=intersect(model.rxns,rxns);

%Assemble all the info
names = model.rxnNames(idx);
genes = model.grRules(idx);
formulas = printRxnFormula(model,rxns,'False');
fluxes = num2cell(solution.x(idx));
lb = num2cell(model.lb(idx));
ub = num2cell(model.ub(idx));
minFluxes = num2cell(minFlux(idx));
maxFluxes = num2cell(maxFlux(idx));
A = [rxns,fluxes,minFluxes,maxFluxes,formulas,genes,names,lb,ub];
xlswrite(workbook,A,'NAD_NADP Couplets','A2');        

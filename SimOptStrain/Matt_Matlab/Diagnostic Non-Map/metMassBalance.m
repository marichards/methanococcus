function [rxns,amounts,directions,names] = metMassBalance(model,met,solution,threshold,printFlag)

%Take in a model, a metabolite, and an FBA solution.  Find all the major
%sinks and sources of that metabolite by 
%
%Inputs
%model - a COBRA model structure
%met - a metabolite of interest
%solution - an FBA solution of the model
%
%Outputs
%rxns - list of reactions that are a major sink or source for metabolite of
%interest
%amount - fluxes from FBA solution multiplied by reaction coefficients that correspond to 'rxns'  
%directions - either 'sink' or 'source' corresponding to each reaction

%Step 1: Find all the reactions for the metabolite
rxns = findRxnsFromMets(model,met);

%Step 2: Find the index of the metabolite for finding coefficients
[~,met_idx] = intersect(model.mets,met);

%Step 3: Find reaction indices
[rxns,rxn_idx] = intersect(model.rxns,rxns);

%%%Calculate the amount (coeff*flux)
%Step 4: Pull out fluxes
fluxes = solution.x(rxn_idx);

%Step 5: Pull out metabolite coefficients for rxns (transpose for
%dimensions and convert to non-sparse
coeffs = full(model.S(met_idx,rxn_idx))';

%Step 6: Multiply fluxes by coeffs to get amounts
amounts = fluxes.*coeffs;

%Step 7: Create a tuneable threshold percentage of maximum flux to grab
%"notable" fluxes for
if (nargin<4)
    threshold = 0.1;
end

%Step 8: Use threshold to create a binary vector of fluxes above it
above = abs(amounts)>=threshold*max(abs(amounts));

%Step 9: Pare down the rxns and amounts to only values above the threshold
rxns = rxns(above);
amounts = amounts(above);

%Step 10: Create a directions vector to store 'sink' or 'source'
directions = cell(length(amounts),1);

%Step 11: For each reaction, designate as a source or sink based on the
%amount
for i = 1:length(amounts)
    if amounts(i)>0
        directions{i} = 'source';
    else
        directions{i} = 'sink';
    end
end

%Step 12: Look at arguments; if no printFlag, make it 'false'
if (nargin<5)
    printFlag=false;
end

%Step 13: If printFlag is true, print out in columns
if printFlag
    fprintf('\n')
    for i=1:length(rxns)
        fprintf('%s\t%f\t%s\n',rxns{i},amounts(i),directions{i})
    end
    fprintf('\n')
end



%Step 14: Find the reaction names
[~,idx] = intersect(model.rxns,rxns);
names = model.rxnNames(idx);

end

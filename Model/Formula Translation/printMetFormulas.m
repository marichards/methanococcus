function printMetFormulas(model,rxn)
%
%Take in a model and reaction, print the formulas for all metabolites in
%the reaction
%
%Inputs
%model - a COBRA model
%rxn - a reaction in the COBRA model
%
%Outputs
%Prints the metabolites and formulas to the screen

%Find the metabolites from the reaction
mets = findMetsFromRxns(model,rxn);
%Clear some print space
fprintf('\n\nReaction metabolites and formulas\n');
%Run through the metabolites
for i=1:length(mets)
   %Find the metabolite in the model
   [~,idx] = intersect(model.mets,mets{i});
   %Print the metabolite and the formula
   fprintf('%s\t%s\n',model.mets{idx},model.metFormulas{idx})
    
end
    
end
function [model,unchanged] = addBounds(model)

%For all reactions, examine the lower and upper bounds, attach them
%accordingly

%Load the bounds
load('2014_10_13_bounds.mat');

%Intersect the model indices with the bounds ones (the dictionary)
[~,model_idx,dict_idx] = intersect(model.rxns,rxn_labs);

%Loop through the indices and put things in
for i = 1:length(model_idx)
    
    %Make the bounds match up
    model.lb(model_idx(i)) = bounds(dict_idx(i),1);
    model.ub(model_idx(i)) = bounds(dict_idx(i),2);
    
    %Change model.rev values too; right now all are reversible
    if model.lb(model_idx(i)) == 0
        
        %Change rev to 0
        model.rev(model_idx(i)) = 0;
    end
        
    
end

%Finally, return reactions that didn't show up
unchanged = setdiff(model.rxns,rxn_labs);

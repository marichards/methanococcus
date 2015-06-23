function model = addDG(model,rxn_list,dG_list)

%Written 02/28/2013 by Matthew Richards

%INPUTS
%model: a COBRA model, can have field model.freeEnergy but not required
%rxn_list: a reaction ID or list of reaction IDs for which free energy of
%formation values will be added (all reactions are exchanges!)
%dG_list: an value or array of values for free energy of formation,
%corresponding to the reaction IDs in rxn_list

%OUTPUTS
%model: a COBRA model with the field model.freeEnergy, which containts free
%energy of formation values

%Check to see if the freeEnergy field is in the model, if not, create it
if ~isfield(model,'freeEnergy')
    model.freeEnergy = zeros(length(model.rxns),1);
end


%Check to see if rxn_list and dG_list are there
if nargin >2

    %Check to be sure all reactions are in the model
    [rxns,rxn_idx] = intersect(rxn_list,model.rxns,'stable');
    if ~isequal(length(rxns),length(model.rxns))
        different = setdiff(rxn_list,model.rxns);
        for j = 1:length(different)
            fprintf('Reaction %s not in model\n',different{j});
            
        end
        
        warning('Above reactions are not in the model and have not been added')
        %Remove the different ones using intersect from before
        rxn_list = rxns;
        dG_list = dG_list(rxn_idx);
        
    end

    %Check to see if rxn_list is a string or a list
    if isa(rxn_list,'cell')
        %Check to see if rxn_list and dG_list are the same length
        if ~isequal(length(rxn_list),length(dG_list))
        error('Reaction list and free energy list must have same length');
        end
        %Add a list of values
        for i = 1:length(rxn_list)
            [~,idx] = intersect(model.rxns,rxn_list{i});
            model.freeEnergy(idx) = dG_list(i);     
        end

    %Add one value
    elseif isa(rxn_list,'char')
        %Check to make sure dG isn't an array
        if length(dG_list)==1
            [~,idx] = intersect(model.rxns,rxn_list);
            model.freeEnergy(idx) = dG_list;
        end
    %If it's not a list of strings or string, break things    
    else
        error('Reaction list must be a cell array or string');
    end
end



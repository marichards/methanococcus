function model = addDG(model,rxn_list,dG_list)

% Takes in a model and adds the freeEnergy field with the supplied values
% for the supplied exchange reactions. If the field already exists, it
% only adds the supplied values for the supplied exchange reactions.
%
% INPUT
% model: a COBRA Toolbox model, which can optionally have field
% model.freeEnergy but is not required possess it. 
% rxn_list: an exchange reaction ID or list of exchange reaction IDs for
% which free energy of formation values will be added
% dG_list: an value or array of values for free energy of formation,
% corresponding to the reaction IDs in rxn_list
%
% OUTPUT
% model: the supplied model with the field model.freeEnergy, which contains
% values for free energy of formation for supplied exchange reactions
%
% Matthew Richards, 02/28/2013


% Check to see if the freeEnergy field is in the model, if not, create it
if ~isfield(model,'freeEnergy')
    model.freeEnergy = zeros(length(model.rxns),1);
end


% Check to see if rxn_list and dG_list are there
if nargin >2
    % Check to be sure all reactions are in the model
    [rxns,rxn_idx] = intersect(rxn_list,model.rxns,'stable');
    if ~isequal(length(rxns),length(model.rxns))
        different = setdiff(rxn_list,model.rxns);
        if length(different)>0
            for j = 1:length(different)
            fprintf('Reaction %s not in model\n',different{j});
            end
            warning('Above reactions are not in the model and have not been added')
        end              
        % Remove the different ones using intersect from before
        rxn_list = rxns;
        dG_list = dG_list(rxn_idx);     
    end
    % Check to see if rxn_list is a string or a list
    if isa(rxn_list,'cell')
        % Check to see if rxn_list and dG_list are the same length
        if ~isequal(length(rxn_list),length(dG_list))
        error('Reaction list and free energy list must have same length');
        end
        % Add a list of values
        for i = 1:length(rxn_list)
            [~,idx] = intersect(model.rxns,rxn_list{i});
            model.freeEnergy(idx) = dG_list(i);     
        end

    % Add one value
    elseif isa(rxn_list,'char')
        % Check to make sure dG isn't an array
        if length(dG_list)==1
            [~,idx] = intersect(model.rxns,rxn_list);
            model.freeEnergy(idx) = dG_list;
        end
    % If it's not a list of strings or string, break things    
    else
        error('Reaction list must be a cell array or string');
    end
else 
    error('Must supply reactions and free energy values')
end
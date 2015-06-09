function createJSONModel(model,model_filename)

% Take in a COBRA model structure and write it to a JSON file for use in
% COBRApy or in Escher

% JSON format is as follows: 
% Reactions:
% {"reactions": [{"subsystem": "Pyruvate Metabolism", "name": "acetaldehyde
% dehydrogenase (acetylating)", "upper_bound": 1000.0, "lower_bound":
% -1000.0, "notes": {}, "metabolites": {"accoa_c": 1.0, "acald_c": -1.0,
% "coa_c": -1.0, "h_c": 1.0, "nad_c": -1.0, "nadh_c": 1.0}, "objective_coefficient": 0.0, "variable_kind": "continuous", "id": "ACALD", "gene_reaction_rule": "(b0351 or b1241)"},
% {"subsystem": "Transport, Extracellular", "name": "acetaldehyde reversible transport",
% "upper_bound": 1000.0, "lower_bound": -1000.0, "notes": {}, 
% "metabolites": {"acald_c": 1.0, "acald_e": -1.0}, "objective_coefficient": 0.0, "variable_kind": "continuous", "id": "ACALDt", "gene_reaction_rule": "s0001"},...
% End of Reactions and model:
% ": "TPI", "gene_reaction_rule": "b3919"}], "description": "Ecoli_core_model", "notes": {},
% Genes
% "genes": [{"name": "b1241", "id": "b1241"}, {"name": "b0351", "id": "b0351"}, 
% End of Genes and Metabolites:
% "b3919", "id": "b3919"}], "metabolites": [{"name": "H", "notes": "{}", "annotation": "{}", "_constraint_sense": "E", "charge": "0", "_bound": "0.0", "formula": "H", "compartment": "c", "id": "h_c"}, 
% {"name": "Nicotinamide-adenine-dinucleotide", "notes": "{}", "annotation": "{}", "_constraint_sense": "E", "charge": "0", "_bound": "0.0", "formula": "C21H26N7O14P2", "compartment": "c", "id": "nad_c"},
% End of Metabolites and Model:
% ": "c", "id": "s7p_c"}], "id": "Ecoli_core_model"}

% Written by Matt Richards, 6/3/2015
%%%%%%%%%%%%%%%%%%

% First, initiate the file:
file_id = fopen(sprintf('%s.json',model_filename),'w');

% And make a non-sparse matrix
A = full(model.S);

% PRINT THE REACTIONS
% Initiate the reactions section
fprintf(file_id,'{"reactions": [');
% Loop through the reactions; write things in chunks
for i=1:length(model.rxns)
    % Print out the subsytem, if it exists, the name, and upper bound
    fprintf(file_id,'{"subsystem": "%s", "name": "%s", "upper_bound": %f, ',...
        model.subSystems{i},model.rxnNames{i},model.ub(i));
    % Print out the lower bound and an empty field for notes, and header for
    % metabolites that we'll fill next
    fprintf(file_id,'"lower_bound": %f, "notes": {}, "metabolites": {',...
        model.lb(i));
    
    % Use A matrix to print out coefficients for metabolites
    % Find the things with coefficients in that reaction
    idx = find(model.S(:,i));
    
    % For each thing in the index 
    for j=1:length(idx)
        % Print out the name of the met and its coefficient
        fprintf(file_id,'"%s": %f',model.mets{idx(j)},A(idx(j),i));
        % If it's not the last thing in the list, print a comma
        if j<length(idx)
            fprintf(file_id,', ');
        % Otherwise, end the metabolites section of reactions
        else
            fprintf(file_id,'}, ');
        end
    end
    
    % Print the objective coefficient, variable kind, ID, and GPR
    fprintf(file_id,'"objective_coefficient": %f, "variable_kind": "continuous", ',...
        model.c(i));
    fprintf(file_id,'"id": "%s", "gene_reaction_rule": "%s"}',...
        model.rxns{i},model.grRules{i});
    % Decide what to print after the gene_reaction_rule
    if i<length(model.rxns)
        fprintf(file_id,', ');
    % Otherwise, end the genes section
    else
        fprintf(file_id,'], ');
    end
end

% FINISH THE REACTIONS SECTION AND PRINT THE MODEL DESCRIPTION AND EMPTY NOTES FIELD
fprintf(file_id,'"description": "%s", "notes": {}, ',model.description);

% PRINT THE GENES
fprintf(file_id,'"genes": [');
% Loop through the genes
for i = 1:length(model.genes)
    % For each, print the name and id (they're the same here)
    fprintf(file_id,'{"name": "%s", "id": "%s"}',...
        model.genes{i},model.genes{i});
    % If it's the not the last thing, print a comma
    if i<length(model.genes)
        fprintf(file_id,', ');
    % Otherwise, end the genes section
    else
        fprintf(file_id,'], "metabolites": [');
    end
end

% PRINT THE METABOLITES
% Loop through the list
for i = 1:length(model.mets)
    % For each, print out its name, notes(empty), annotation(empty),
    fprintf(file_id,'{"name": "%s", "notes": "{}", "annotation": "{}", ',...
        model.metNames{i});
    % _constraint_sense, charge,
    fprintf(file_id,'"_constraint_sense": "E", "charge": "%d", ',...
        model.metCharge(i));
    % _bound, formula, 
    fprintf(file_id,'"_bound": "0.0", "formula": "%s", ',...
        model.metFormulas{i});
    % Print correct compartment 
    if regexp(model.mets{i},'_c0')
        fprintf(file_id,'"compartment": "c", ');
    else
        fprintf(file_id,'"compartment": "e", ');
    end
    % Print the ID
    fprintf(file_id,'"id": "%s"}',model.mets{i});
    % If it's the not the last thing, print a comma
    if i<length(model.mets)
        fprintf(file_id,', ');
    % Otherwise, end the metabolites section
    else
        fprintf(file_id,'], ');
    end    
end

% PRINT THE MODEL ID
fprintf(file_id,'"id": "%s"}',model.id');

end
    
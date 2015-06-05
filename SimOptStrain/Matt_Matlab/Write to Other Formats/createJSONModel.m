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
    fprintf(file_id,'{"subsystem": "%s", "name": "%s", "upper_bound: %0.1f, ',...
        model.subSystems{i},model.rxnNames{i},model.ub(i));
    % Print out the lower bound and an empty field for notes, and header for
    % metabolites that we'll fill next
    fprintf(file_id,'"lower_bound": %0.1f, "notes": {}, "metabolites": {',...
        model.lb(i));
    
    % Use A matrix to print out coefficients for metabolites
    % Find the things with coefficients in that reaction
    idx = find(model.S(:,1));
    
    % For each thing in the index 
    for j=1:length(idx)
        % Print out the name of the met and its coefficient
        fprintf(file_id,'"%s": %0.1f',model.mets{idx(j)},A(idx(j),i));
        % If it's not the last thing in the list, print a comma
        if j<length(idx)
            fprintf(file_id,', ');
        else
            fprintf(file_id,'}, ');
        end
    end
    
    % Print the objective coefficient, variable kind, ID, and GPR
    fprintf(file_id,'"objective_coefficient": %f, "variable_kind": "continuous", ',...
        model.c(i));
    fprintf(file_id,'"id": "%s", gene_reaction_rule": "%s"}, ',...
        model.rxns{i},model.grRules{i});
    
end

% PRINT THE MODEL DESCRIPTION AND EMPTY NOTES FIELD
fprintf(file_id,'"description": "%s", "notes": {}, ',model.description);
    
%%%%%ALL THAT'S LEFT IS GENES AND METABOLITES%%%%%
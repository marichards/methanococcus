function createTSVModelFiles(model,rxns_filename,genes_filename)
%%
% Takes in a COBRA model and creates TSV files. One lists all reactions and
% associated metadata; the other lists all genes and associated reactions
%
% INPUT
% model: a COBRA Toolbox model structure
% rxns_filename: a supplied name for the reactions file created
% genes_filename: a supplied name for the genes file created
%
%Written by Matthew Richards 06/08/2016

%Create the rxns file
rxn_file_id = fopen(sprintf('%s.txt',rxns_filename),'w');
%Print the headers
fprintf(rxn_file_id,'Rxn_ID\tRxn_Name\tGPR\tSubsystem\tFormula\tReadable_Formula\n');
%Create a full index of reactions
rxn_idx = (1:length(model.rxns))';
%Find the biomass equation index using the objective
bio_idx = find(model.c);
%Take out the biomass index
rxn_idx = setdiff(rxn_idx,bio_idx,'stable');
%Add the biomass index back to the top
rxn_idx = [bio_idx;rxn_idx];

% Change genes to kbase IDs
%%%% No longer necessary as of 05/05
% model = changeGenesToKbase(model);

%Now we have the correct index order; create the correct fields
%Loop through the index
for i=1:length(rxn_idx)
    %Make sure to use the number in the index for grabbing each reaction
    %First grab the reaction formula
    formula = printRxnFormula(model,model.rxns{rxn_idx(i)},false);
    % And the human-readable version of the formula
    readable_formula = printRxnFormula(model,model.rxns{rxn_idx(i)},false,false,true);

    %Now we have everything; print it out according to the headers
    fprintf(rxn_file_id,'%s\t%s\t%s\t%s\t%s\t%s\t\n',...
        model.rxns{rxn_idx(i)},model.rxnNames{rxn_idx(i)}...
        ,model.grRules{rxn_idx(i)},model.subSystems{rxn_idx(i)},...
        formula{1},readable_formula{1});
end

% Create the genes file
gene_file_id = fopen(sprintf('%s.txt',genes_filename),'w');
% Print the headers
fprintf(gene_file_id,'Gene\tRxns\n');

% For each gene in the model:
for i = 1:length(model.genes)
    % Find the list of associated reactions and put them in a cell array
    rxns = struct2cell(findRxnsFromGenes(model,model.genes{i}));
    % Find the size of it, aka the number of reactions
    [m,~] = size(rxns{1,1});
    % Now print things
    fprintf(gene_file_id,'%s\t',model.genes{i});
    
    % Catch 1-reaction cases
    if m ==1 
        fprintf(gene_file_id,'%s',rxns{1,1}{1}{1});
    else
        % Loop through the associated reactions and print them
        for j = 1:m-1
            fprintf(gene_file_id,'%s;',rxns{1,1}{j}{1});
        end
        fprintf(gene_file_id,'%s',rxns{1,1}{m}{1});
    end
    
    % Go to the next line
    fprintf(gene_file_id,'\n');
end


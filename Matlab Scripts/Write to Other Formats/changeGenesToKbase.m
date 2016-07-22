function model = changeGenesToKbase(model)

% Note: this takes just under a minute to run
% Correct mapping: load the inverted gene dictionary container.Map
load('gene_dictionary.mat');

% For each gene in the model, check if it's in the keys
for i = 1:length(model.genes)
    if ismember(model.genes{i},keys(gene_dict))
        % If it is, then replace it with its value; else leave it
        model.genes{i} = gene_dict(model.genes{i});
    end
end

% Do the same for grRules, but use regular expressions
% This takes about 50 seconds
% Grab all IDs in reverse numerical order
all_kbase_ids = keys(gene_dict);
for i = 1:length(all_kbase_ids)
    % Translate the Kbase ID into the correct patter
    pattern = regexptranslate('escape',all_kbase_ids{i});
    % Nest in a loop through all grRules
    for j = 1:length(model.grRules)
        % Replace any hits with the value in the dictionary
        model.grRules{j} = regexprep(model.grRules{j},pattern,...
            gene_dict(all_kbase_ids{i}));
    end
end


% Outdated function; use dictionary method instead
% % First change the genes themselves
% for i = 1:length(model.genes)
%     %Catch 3-zero cases
%     if regexp(model.genes{i},'mmp000')
%         model.genes{i} = regexprep(model.genes{i},'mmp000','kb\|g.575.peg.');
%     end
%     %Catch 2-zero cases
%     if regexp(model.genes{i},'mmp00')
%         model.genes{i} = regexprep(model.genes{i},'mmp00','kb\|g.575.peg.');
%     end
%     %Catch 1-zero cases
%     if regexp(model.genes{i},'mmp0')
%         model.genes{i} = regexprep(model.genes{i},'mmp0','kb\|g.575.peg.');
%     end    
%           
% 	%Catch no-zero cases
% 	model.genes{i} = regexprep(model.genes{i},'mmp','kb\|g.575.peg.');
% end
% 
% % Do the same for grRules
% for i = 1:length(model.grRules)
%     %Catch 3-zero cases
%     if regexp(model.grRules{i},'mmp000')
%         model.grRules{i} = regexprep(model.grRules{i},'mmp000','kb\|g.575.peg.');
%     end
%     %Catch 2-zero cases
%     if regexp(model.grRules{i},'mmp00')
%         model.grRules{i} = regexprep(model.grRules{i},'mmp00','kb\|g.575.peg.');
%     end
%     %Catch 1-zero cases
%     if regexp(model.grRules{i},'mmp0')
%         model.grRules{i} = regexprep(model.grRules{i},'mmp0','kb\|g.575.peg.');
%     end              
% 	%Catch no-zero cases
% 	model.grRules{i} = regexprep(model.grRules{i},'mmp','kb\|g.575.peg.');
% end
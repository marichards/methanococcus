function model = addChEBIIDs(model)

% Go through the model and for anything with a KEGG ID, translate to a
% ChEBI ID using the obtained list of cross references


% First, create a new field in the model for ChEBI IDs
model.metChEBIID = cell(length(model.mets),1);

% Load the chebi/kegg translation dictionary
load('2015_09_30_chebi_to_kegg.mat')
% Go through the list and find things with KEGG IDs
for i = 1:length(model.mets)
if ~isempty(model.metKEGGID{i})
    % If they have KEGG IDs, then pull them out and split them on the "|"
    % symbol that may or may not be there
    keggs = strsplit(model.metKEGGID{i},'|');
    % Now I have a list of KEGG IDs. And it's a cell. For each of these,
    % find the corresponding ChEBI ID and add these

    % For all but the end, add a "|"
    for j = 1:length(keggs)-1
        % Find the corresponding ChEBI ID
        [~,kegg_idx] = intersect(kegg_ids,keggs{j}) ;
        % Add it with the bar to the correct index
        model.metChEBIID{i} = strcat(model.metChEBIID{i},...
            chebi_ids{kegg_idx},'|');        
    end
            
    % For the end one, add just the ID
    [~,kegg_idx] = intersect(kegg_ids,keggs{end});
    model.metChEBIID{i} = strcat(model.metChEBIID{i},chebi_ids{kegg_idx});
end
end
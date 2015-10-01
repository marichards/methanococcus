function model = addChEBIIDs(model)

% Go through the model and for anything with a KEGG ID, translate to a
% ChEBI ID using the obtained list of cross references


% First, create a new field in the model for ChEBI IDs
model.metChEBIID = cell(length(model.mets),1);

% Go through the list and find things with KEGG IDs
for i = 1:length(model.mets)
if ~isempty(model.metKEGGID{i})
    % If they have KEGG IDs, then pull them out and split them on the "|"
    % symbol that may or may not be there
    keggs = strplit(model.metKEGGID{i},'\');
    % Now I have a list of KEGG IDs. For each of these, find the
    % corresponding ChEBI ID and add these
    
    % For all but the end, add a "|"
    for j = 1:length(keggs)-1
        
    end
    
    % For the end one, add just the ID
[~,idx] = intersect(kegg_ids,model.metKEGGID{i});
chebi = chebi_ids(idx)
end
end
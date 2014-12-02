function matches = findMetabolomicsMets(model)

%Note: KEGG IDs alone aren't perfect; some may not match because of
%stereochemistry (e.g. L-Phenylalanine instead of Phenylalanine)

%Read in the metabolomics KEGG IDs
[~,text] = xlsread('Targeted Metabolomics - Aqueous 8-2014.xls','IDs & Pathways','B2:B172');

%We're going to have to actually look using a string match (regex)
%Create an empy array
matches = {};
%Loop through the items of the text
for i=1:length(text)
    %For each, search through and find a match
    for j=1:length(model.metKEGGID)
        %Match it
        if regexp(model.metKEGGID{j},text{i})
            %Add it to the matches
            matches = [matches;model.mets{j}];
            break
        end
    end
end
function extMets = countExternalMets(model)

%Go through model and count all the external metabolites
%
%Input: a COBRA model
%
%Output: all the external metabolites
%

%Create an array to write to
extMets = {};

%Loop through metabolites
for i=1:length(model.mets)
   
    %Check if there's an e0
    if regexp(model.mets{i},'e0')
        %Add it
        extMets=[extMets; model.mets{i}];
    end
    
    
end



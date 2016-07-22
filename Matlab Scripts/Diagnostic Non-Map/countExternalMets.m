function extMets = countExternalMets(model,external_compartment)

% Goes through COBRA model and grabs all the external metabolites
%
% INPUT
% model: a COBRA Toolbox model structure
%
% OPTIONAL INPUT
% external_compartment: notation for the extracelluar compartment used in
% the supplied model. (Default = 'e0')
%
% OUTPUT 
% extMets: all the external metabolites
%
% Matthew Richards, 09/24/2015

% Check if there's a different external compartment specified
if nargin < 2
    external_compartment = 'e0';
end

% Create an array to write to
extMets = {};

% Loop through metabolites
for i=1:length(model.mets)
   
    % Check if there's an e0
    if regexp(model.mets{i},external_compartment)
        % Add it if there is
        extMets=[extMets; model.mets{i}];
    end
    
    
end



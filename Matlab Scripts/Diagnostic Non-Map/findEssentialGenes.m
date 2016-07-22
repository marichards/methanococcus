function essentials = findEssentialGenes(model,genes,cutoff)

% Based on a cutoff for percentage of growth considered "lethal", finds the
% essential genes in the given list of genes
%
% INPUT
% model: a COBRA Toolbox model structure
%
% OPTIONAL INPUT
% genes: a list of genes in the supplied model to evaluate for essentiality
% (Default = model.genes)
% cutoff: a supplied cutoff for the ratio of knockout growth to wild-type
% growth that signifies lethality. Scales from 0-1. (Default = 0)
% 
% OUTPUT
% essentials: subset of supplied "genes" that are essential for growth
% based upon the supplied cutoff
%
% Matthew Richards, 09/24/2014


% Qualify if there isn't a cutoff; make it 0
if (nargin<3)
    cutoff = 0;
end

% If there's no genes given, then just use all of them
if (nargin<2)
    genes = model.genes;
end

% Run the single gene deletion
grRatio = singleGeneDeletion(model,'FBA',genes);

% Pull out things where the growth ratio is below the cutoff
essentials = genes(grRatio<=cutoff);
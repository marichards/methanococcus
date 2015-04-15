function [missing,model] = findMissingBiomassMetsMOD(model)

% Take in a model that doesn't grow and find the metabolite(s) that are
% missing from the biomass. 

% First pull out the list of all metabolites required in the biomass
% Pull out the biomass index
idx = find(model.c);
% Find the indices of metabolites in biomass reaction that are less than 0
% (the reactants)
required_mets = model.mets(model.S(:,idx)<0);

% Take out all of the required mets from the biomass; this should allow
% growth no matter what

% Find the biomass reaction index
[~,bio_idx] = intersect(model.rxns,'biomass0');

% Find index of all metabolites from the biomass
[~,idx] = intersect(model.mets,required_mets);
% But save them first
coeffs = model.S(idx,bio_idx);
% Then remove them
model.S(idx,bio_idx) = 0;

% Now try adding them one at a time. Shouldn't take any randomization at
% all, so just do it in order

% Create an array to keep the missing mets in
missing = {};
% Now loop through each metabolite
for i=1:length(required_mets)
    % Try putting it back into the biomass
    model.S(idx(i),bio_idx) = coeffs(i);
    % Test for growth
    solution = optimizeCbModel(model);
    % If it doesn't grow 
    if solution.f <1e-10
        % Then remove it from the biomass and keep it in an array
        model.S(idx(i),bio_idx) = 0;
        missing = [missing;required_mets{i}];
    end
end




   
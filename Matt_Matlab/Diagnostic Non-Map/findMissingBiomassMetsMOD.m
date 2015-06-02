function [missing,model] = findMissingBiomassMetsMOD(model)

% Take in a model that doesn't grow and find the metabolite(s) that are
% missing from the biomass. 

% First pull out the list of all metabolites required in the biomass
% Pull out the biomass index
idx = find(model.c);
% Find the indices of metabolites in biomass reaction that are less than 0
% (the reactants)
required_mets = model.mets(model.S(:,idx)<0);

% Find the biomass reaction index
[~,bio_idx] = intersect(model.rxns,'biomass0');

% Find index of all metabolites from the biomass
[~,idx] = intersect(model.mets,required_mets);
% But save them first
coeffs = model.S(idx,bio_idx);

% Now try KOing 1 biomass component at a time and testing for growth
% Create an array to keep the missing mets in
for i=1:length(required_mets)
    %Create a test model
    test = model;
    test.S(idx(i),bio_idx) = 0;
    % Test for growth
    solution = optimizeCbModel(model);
    % If it grows, break the loop
    if solution.f > 1e-5
        solved = 1;
        missing = required_mets{i};
        model = test;
        break
    else
        solved = 0;
        
    end
end

if ~solved
    for i = 1:length(required_mets)-1
        
        for j = (i+1):length(required_mets)
            % Make a model with 2 knockouts
            test = model;
            test.S(idx(i),bio_idx) = 0;
            test.S(idx(j),bio_idx) = 0;
            solution = optimizeCbModel(model);
    % If it grows, break the loop
    if solution.f > 1e-5
        solved = 1;
        missing = {required_mets{i},required_mets{j}};
        model = test;
        break
    else
        solved = 0;
    end
        end
    end
end
        

% Function end
end



   
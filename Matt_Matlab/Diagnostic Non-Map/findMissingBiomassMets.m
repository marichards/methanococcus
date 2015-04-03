function [missing,model] = findMissingBiomassMets(model)

% Take in a model that doesn't grow and find the metabolite(s) that are
% missing from the biomass. 

% First pull out the list of all metabolites required in the biomass
% Pull out the biomass index
idx = find(model.c);
% Find the indices of metabolites in biomass reaction that are less than 0
% (the reactants)
required_mets = model.mets(model.S(:,idx)<0);

% Test out adding one exchange at a time, until it works

%for i=1:length(required_mets)
    
%   new_model = addReaction(model,sprintf('%s sink',required_mets{i}),...
%        sprintf('%s <=> ',required_mets{i}));
   
%    solution = optimizeCbModel(new_model);
    
%    if solution.f > 1e-10
        
%        missing_met = required_mets{i};
%        fprintf('Missing biomass component is %s',missing_met);
        
%        break
%    else
%        missing_met = 0;
%    end
%end


% There may be more than 1....if there's no missing_met, then return that
%there's more than once missing
%if ~missing_met
%    fprintf('More than one metabolite is missing');
%end

% Try the opposite way: add exchanges for all of them and then try taking them away one at a time
% Loop through all required biomass metabolites
added_exchanges = {};
for i=1:length(required_mets)
   % Add an exchange for each one 
   model = addReaction(model,sprintf('%s sink',required_mets{i}),...
        sprintf('%s <=> ',required_mets{i}));
   %Catalogue them as I go along
   added_exchanges = [added_exchanges; sprintf('%s sink',required_mets{i})];
end

% Now try taking them away one at a time, doing that until I can't take one
% away
% First, make an array to catalogue the reactions we MUST keep:
missing = {};
% Now start a while loop; we want to keep taking out reactions until
% there's no more in our list of added exchanges
while length(added_exchanges) > 0
    % For each iteration, take a random number between 1 and the length of
    % added_exchanges
    idx = ceil(length(added_exchanges)*rand);
    % Remove that reaction
    test_model = removeRxns(model,added_exchanges{idx});
    % Now test for growth
    solution = optimizeCbModel(test_model);
    % If it still grows, then the reaction isn't necessary, so toss it
    if solution.f > 1e-10
        % Make the model the new one
        model = test_model;
    else
        %If it doesn't grow, it's necessary; put it in missing and keep the
        %model
        missing=[missing;added_exchanges{idx}];
    end
    %Remove the reaction from the queue
    added_exchanges = setdiff(added_exchanges,added_exchanges{idx});
end


   
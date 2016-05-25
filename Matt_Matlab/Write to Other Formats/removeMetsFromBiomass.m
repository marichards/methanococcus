function model = removeMetsFromBiomass(model,mets)


% Find biomass index
bio_idx = find(model.c);
% Cycle through mets
for i = 1:length(mets)
    [~,met_idx] = intersect(model.mets,mets{i});
    model.S(met_idx,bio_idx) = 0;
end
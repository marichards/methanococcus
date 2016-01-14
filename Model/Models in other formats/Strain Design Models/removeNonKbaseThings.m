function model = removeNonKbaseThings(model)

% As a pre-processing step, remove compounds that aren't in Kbase from the
% media. Also substitute the arbitrary compound "cpd15000" for dG, as it's
% also an issue

% First, remove the 3 problematic exchanges
model = removeRxns(model,{'EX_NAC[c0]','EX_Membrane_lipid[c0]','EX_Flagellin[e0]'});

% Now remove archaellin from the biomass
[~,met_idx] = intersect(model.mets,'ARCN[e0]');
[~,bio_idx] = intersect(model.rxns,'biomass0');
model.S(met_idx,bio_idx) = 0;

% Substitute info for dG
[~,met_idx] = intersect(model.mets,'dG');
if met_idx
    model.mets{met_idx} = 'cpd15000[c0]';
end
% Do 
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

%%New Part
% Find non-Seed compounds
logs = regexp(model.mets,'cpd[0-9]{5}');
non_Seed = {};
for i=1:length(logs)
if isempty(logs{i})
non_Seed = [non_Seed;model.mets{i}];
end
end

% Find reactions that use those compounds
rxns = setdiff(findRxnsFromMets(model,non_Seed),'biomass0');
model = removeRxns(model,rxns);

% Remove all extra things from biomass
to_remove = {'cpd03425[c0])','cpd02817[c0])','cpd00895[c0]','cpd02246[c0]',...
    'cpd00643[c0]','cpd00649[c0]','cpd15868[c0]','cpd00131[c0]','cpd16579[c0]',...
    'SATARCHL[c0]','SATARCHLS[c0]'};

for i=1:length(to_remove)
    [~,met_idx] = intersect(model.mets,to_remove{i});
    [~,bio_idx] = intersect(model.rxns,'biomass0');
    model.S(met_idx,bio_idx) = 0;
end

% Remove more things from the biomass
to_remove = {'cpd00161_c0','cpd00033_c0','cpd02817_c0','cpd11493_c0',...
    'cpd00042_c0','cpd00166_c0'};

for i=1:length(to_remove)
    [~,met_idx] = intersect(model.mets,to_remove{i});
    [~,bio_idx] = intersect(model.rxns,'biomass0');
    model.S(met_idx,bio_idx) = 0;
end

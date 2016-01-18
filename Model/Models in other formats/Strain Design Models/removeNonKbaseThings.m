function model = removeNonKbaseThings(model)

% As a pre-processing step, remove compounds that aren't in Kbase from the
% media. Also substitute the arbitrary compound "cpd15000" for dG, as it's
% also an issue

% First, remove the 3 problematic exchanges
model = removeRxns(model,{'EX_NAC[c0]','EX_Membrane_lipid[c0]','EX_Flagellin[e0]'});

% Now remove archaellin from the biomass
model = removeMetsFromBiomass(model,{'ARCN[e0]'});

% Substitute info for dG
%[~,met_idx] = intersect(model.mets,'dG');
%if met_idx
%    model.mets{met_idx} = 'cpd15000[c0]';
%end

% Remove dG from model and GIBBS rxn
if ismember('dG',model.mets)
    model = removeMetabolites(model,'dG',false);
    model = removeRxns(model,'GIBBS_kJ/GDW');
end

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

model = removeMetsFromBiomass(model,to_remove);

% Remove more things from the biomass
to_remove = {'cpd00161[c0]','cpd00033[c0]','cpd02817[c0]',...
    'cpd00042[c0]','cpd03425[c0]'};

model = removeMetsFromBiomass(model,to_remove);

% Compile list of biomass compounds that are okay (35 of them right now)
kbase_ok = {'cpd11416[c0]','cpd00063[c0]','cpd10515[c0]','cpd00001[c0]',...%4/35
    'cpd00002[c0]','cpd00067[c0]','cpd17041[c0]','cpd00009[c0]',...%8/35
    'cpd00254[c0]','cpd00099[c0]','cpd10516[c0]','cpd00205[c0]',...%12/35
    'cpd00149[c0]','cpd17043[c0]','cpd00008[c0]','cpd17042[c0]',...%16/35
    'cpd00035[c0]','cpd00062[c0]','cpd00264[c0]','cpd00051[c0]',...%20/35
    'cpd00069[c0]','cpd00356[c0]','cpd00060[c0]','cpd00054[c0]',...%24/35
    'cpd00065[c0]','cpd00039[c0]','cpd00132[c0]','cpd00118[c0]',...%28/35
    'cpd00041[c0]','cpd00156[c0]','cpd00357[c0]','cpd00012[c0]',...%32/35
    'cpd00056[c0]','cpd00066[c0]','cpd00016[c0]'};%35/35

% Use this list to find all remaining biomass things the aren't okay and
% remove them
all_biomass = findMetsFromRxns(model,'biomass0');
model = removeMetsFromBiomass(model,setdiff(all_biomass,kbase_ok));
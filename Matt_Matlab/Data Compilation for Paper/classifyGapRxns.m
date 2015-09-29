function [stats,rxns] = classifyGapRxns(model)

% For all of the gapfilling reactions in a model, pull out the subsystems
% associated with them and return the results as a table of statistics and
% as the reactions themselves
%
% INPUT
% model: a COBRA Toolbox model structure
%
% OUTPUT
% stats: a table listing all subsystems of gapfilled (no-gene) internal
% reactions, the number of reactions belonging to each subsystem, and the
% percentage of gapfilled reactions in that subsystem. 
% rxns: a list of the gapfilled internal reactions in the model
%
% Matthew Richards, 09/28/2015

% First, pull out the set of reaction without genes and don't include
% transport reactions. Exchanges are already excluded
rxns = setdiff(findRxnsWOGenes(model),findTransRxnsMOD(model));

% Now find those reactions in the model
[rxns,idx] = intersect(model.rxns,rxns);

% Pull out the subsystems
subs = model.subSystems(idx);

% Put them in table form
stats = tabulate(subs);


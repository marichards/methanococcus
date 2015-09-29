function model = switchToAlanine(model)

% Switches M.maripaludis model to use L-alanine as the nitrogen source
% for growth. Turns on NH3 uptake and turns off NH3/N2 uptakes
%
% INPUT
% model: the M. maripaludis model, a COBRA Toolbox model structure
%
% OUTPUT
% model: a reconfigured version of the model, set to grow with L-alanine as
% the nitrogen source
%
% Matthew Richards, 09/29/2015


% Turn down ammonia
model = changeRxnBounds(model,'EX_cpd00013[e0]',0,'l');

% Turn up alanine
model = changeRxnBounds(model,'EX_cpd00035[e0]',-1000,'l');

% Turn down nitrogen
model = changeRxnBounds(model,'EX_cpd00528[e0]',0,'l');
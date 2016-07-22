function model = switchToN2(model)

% Switches M.maripaludis model to use N2 as the nitrogen source for growth,
% or "nitrogen fixing" conditions. Turns on N2 uptake and turns off
% NH3/alanine uptakes
%
% INPUT
% model: the M. maripaludis model, a COBRA Toolbox model structure
%
% OUTPUT
% model: a reconfigured version of the model, set to grow with N2 as the
% nitrogen source
%
% Matthew Richards, 09/29/2015


% Turn down ammonia
model = changeRxnBounds(model,'EX_cpd00013[e0]',0,'l');

% Turn down alanine
model = changeRxnBounds(model,'EX_cpd00035[e0]',0,'l');

% Turn up nitrogen
model = changeRxnBounds(model,'EX_cpd00528[e0]',-1000,'l');
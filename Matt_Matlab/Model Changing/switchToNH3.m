function model = switchToNH3(model)

% Switches M.maripaludis model to use NH3 (ammonia) as the nitrogen source
% for growth. Turns on NH3 uptake and turns off N2/alanine uptakes
%
% INPUT
% model: the M. maripaludis model, a COBRA Toolbox model structure
%
% OUTPUT
% model: a reconfigured version of the model, set to grow with NH3 as the
% nitrogen source
%
% Matthew Richards, 09/29/2015


% Turn up ammonia
model = changeRxnBounds(model,'EX_cpd00013[e0]',-1000,'l');

% Turn down alanine
model = changeRxnBounds(model,'EX_cpd00035[e0]',0,'l');

% Turn down nitrogen
model = changeRxnBounds(model,'EX_cpd00528[e0]',0,'l');
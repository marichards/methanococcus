function model = switchToH2(model)

% Switches M.maripaludis model to H2 as the main electron source. Also
% ensures that acetate and formate are turned off and protons are not
% supplied
%
% INPUT
% model: the M. maripaludis model, a COBRA Toolbox model structure
%
% OUTPUT
% model: a reconfigured version of the model, set to grow with H2 as the
% electron donor
%
% Matthew Richards, 09/29/2015


% Turn H2 to 45 %Test doing it at 1000
model = changeRxnBounds(model,'EX_cpd11640[e0]',-1000,'l');

% Turn off formate
model = changeRxnBounds(model,'EX_cpd00047[e0]',0,'l');

%Turn off H+ uptake
model = changeRxnBounds(model,'EX_cpd00067[e0]',0,'l');

% Be sure to turn off the acetate input, just in case
model = changeRxnBounds(model,'EX_cpd00029[e0]',0,'l');

% Set a bound on methane
model = changeRxnBounds(model,'EX_cpd01024[e0]',30,'b');

function model = switchToFormateOnly(model)

% Switches M.maripaludis model from H2 to formate as the main electron
% source. Also ensures that acetate is turned off and protons are supplied
%
% INPUT
% model: the M. maripaludis model, a COBRA Toolbox model structure
%
% OUTPUT
% model: a reconfigured version of the model, set to grow with formate as
% the electron donor
%
% Matthew Richards, 09/29/2015


% Turn H2 off
model = changeRxnBounds(model,'EX_cpd11640[e0]',0,'l');

% Turn up formate
model = changeRxnBounds(model,'EX_cpd00047[e0]',-1000,'l');

%Turn on H+ or else it gets no growth
model = changeRxnBounds(model,'EX_cpd00067[e0]',-1000,'l');

% Be sure to turn off the acetate input, just in case
% Turn on acetate
model = changeRxnBounds(model,'EX_cpd00029[e0]',0,'l');

% Set a bound on methane
model = setMethaneSecretion(model,50);
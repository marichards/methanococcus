function model = switchToAcetateOnly(model)

% Switches M.maripaludis model to use H2 as the main electron source and
% supplement with acetate. Also ensures that formate is turned off and
% protons are not supplied
%
% INPUT
% model: the M. maripaludis model, a COBRA Toolbox model structure
%
% OUTPUT
% model: a reconfigured version of the model, set to grow with H2 as the
% main electron donor and acetate as a supplementary source
%
% Matthew Richards, 09/29/2015


% Turn H2 off
model = changeRxnBounds(model,'EX_cpd11640[e0]',0,'l');

% Turn off formate
model = changeRxnBounds(model,'EX_cpd00047[e0]',0,'l');

%Turn off H+ uptake
model = changeRxnBounds(model,'EX_cpd00067[e0]',0,'l');

% Turn acetate on
model = changeRxnBounds(model,'EX_cpd00029[e0]',-1000,'l');

% Set a bound on methane
model = changeRxnBounds(model,'EX_cpd01024[e0]',46,'b');

% Check if model is specific ferredoxins or not and set bound on Eha and
% Ehb in either case
if ismember('Eha/Ehb',model.rxns)
    % If not, then set bounds on Eha/Ehb
    model = changeRxnBounds(model,'Eha/Ehb',4.6,'u');
    model = changeRxnBounds(model,'Eha/Ehb',-4.6,'l');
else
    % If it is, then sent on both Eha and Ehb
    model = changeRxnBounds(model,'Eha',4.6,'u');
    model = changeRxnBounds(model,'Eha',-4.6,'l');
    model = changeRxnBounds(model,'Ehb',4.6,'u');
    model = changeRxnBounds(model,'Ehb',-4.6,'l');
end
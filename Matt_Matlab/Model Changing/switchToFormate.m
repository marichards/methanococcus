function model = switchToFormate(model)

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

% Turn on acetate
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
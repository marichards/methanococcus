function model = setMethaneSecretion(model,secretion_rate)

% Sets methane secretion rate to the specified value, constraining the
% model to produce that much methane. Also restricts Eha/Ehb to 10% of the
% methane secretion rate value, so that they can only play an anaplerotic
% role. 
%
% INPUT
% model: the M. maripaludis model, a COBRA Toolbox model structure
% secretion_rate: desired methane secretion rate for model simulation. This
% rate serves as the main model constraint in all default simulations
%
% OUTPUT
% model: the M. maripaludis model set to the supplied methane secretion
% rate
%
% Matthew Richards, 12/21/2015

% Change the methane secretion bound directly
model = changeRxnBounds(model,'EX_cpd01024[e0]',secretion_rate,'b');

% Now, check if model is specific ferredoxins or not and set bound on Eha and
% Ehb in either case
if ismember('Eha/Ehb',model.rxns)
    % If not, then set bounds on Eha/Ehb
    model = changeRxnBounds(model,'Eha/Ehb',secretion_rate/10,'u');
    model = changeRxnBounds(model,'Eha/Ehb',-secretion_rate/10,'l');
else
    % If it is, then sent on both Eha and Ehb
    model = changeRxnBounds(model,'Eha',secretion_rate/10,'u');
    model = changeRxnBounds(model,'Eha',-secretion_rate/10,'l');
    model = changeRxnBounds(model,'Ehb',secretion_rate/10,'u');
    model = changeRxnBounds(model,'Ehb',-secretion_rate/10,'l');
end
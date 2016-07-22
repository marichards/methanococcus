function model = removeEhaBounds(model)

% The default M. maripaludis model has a bound on Eha/Ehb preventing more
% than 10% the methane flux from flowing through. This function removes
% that bound and allows Eha/Ehb to run at full capacity, up to -1000/1000
%
% INPUT
% model: the M. maripaludis model, a COBRA toolbox model structure,
% presumably with Eha/Ehb constrained to 10% the methane flux. 
%
% OUTPUT
% model: the M. mariapludis model with no 10% constraint on the Eha/Ehb
% reaction(s)

% Check if model is specific ferredoxins or not and unset bounds on Eha and
% Ehb in either case
if ismember('Eha/Ehb',model.rxns)
    % If not, then set bounds on Eha/Ehb
    model = changeRxnBounds(model,'Eha/Ehb',1000,'u');
    model = changeRxnBounds(model,'Eha/Ehb',-1000,'l');
else
    % If it is, then sent on both Eha and Ehb
    model = changeRxnBounds(model,'Eha',1000,'u');
    model = changeRxnBounds(model,'Eha',-1000,'l');
    model = changeRxnBounds(model,'Ehb',1000,'u');
    model = changeRxnBounds(model,'Ehb',-1000,'l');
end
function [qs,dGs] = analyzeQSensitivity(model)
%%
% Takes in the M. maripaludis model and does a sensitivity analysis of
% equilibrium quotient, without addressing specific concentrations

% INPUT
% model: the M. maripaludis model, a COBRA toolbox model structure
%
% OUTPUT
% qs: an array of equilibrium quotient values used for the sensitivity
% analysis
% dGs: the array of free energy predictions corresponding to the "qs" array
%
% Matthew Richards, 10/13/2015
%%

% Set the model to uptake methane and make methanol

% First add in methanol pathway
% Add in the reaction for converting methanol to methyl-coM
model = addReaction(model,{'rxn10568[c0]','Methanol: coenzyme M methyltransferase'},...
    'cpd00116[c0] + cpd02246[c0] <=> cpd02438[c0] + cpd00001[c0]');

% Add in the methanol uptake
model = addReaction(model,{'rxn10570[c0]','Methanol diffusion'},...
    'cpd00116[e0] <=> cpd00116[c0]');
model = addReaction(model,{'EX_cpd00116[e0]','EX_Methanol[e0]'},...
    'cpd00116[e0] <=> ');

% Add metabolite info for methanol
[~,idx] = intersect(model.mets,'cpd00116[c0]');
model.metNames{idx} = 'Methanol[c0]';
model.metCharge(idx)=0;
model.metFormulas{idx}='CH4O';
[~,idx] = intersect(model.mets,'cpd00116[e0]');
model.metNames{idx} = 'Methanol[e0]';
model.metCharge(idx)=0;
model.metFormulas{idx}='CH4O';

% Add it to the free energy
[~,idx] = intersect(model.rxns,'EX_cpd00116[e0]');
model.freeEnergy(idx) = -0.0302;
[~,idx] = intersect(model.rxns,'rxn10570[c0]');
model.freeEnergy(idx) = 0;

% Change bounds to function in reverse
% Change bounds such that methane goes IN instead of OUT
model = changeRxnBounds(model,'EX_cpd01024[e0]',-1000,'l');
model = changeRxnBounds(model,'EX_cpd01024[e0]',0,'u');
% Force methanol to come OUT
model = changeRxnBounds(model,'EX_cpd00116[e0]',1,'l');
model = changeRxnBounds(model,'EX_cpd00116[e0]',1,'u');
% Turn off Hydrogen input and let it come out
model = changeRxnBounds(model,'EX_cpd11640[e0]',0,'l');
model = changeRxnBounds(model,'EX_cpd11640[e0]',1000,'u');
% Try turning off proton output, else it pours out. 
model = changeRxnBounds(model,'EX_cpd00067[e0]',0,'b');

% Be sure formate/acetate and N2/Alanine are off, and that NH3 is on
% Turn down formate
model = changeRxnBounds(model,'EX_cpd00047[e0]',0,'l');
% Turn down acetate
model = changeRxnBounds(model,'EX_cpd00029[e0]',0,'l');
% Make certain that NH3 is on, N2 is off, Alanine is off
% Turn up ammonia
model = changeRxnBounds(model,'EX_cpd00013[e0]',-1000,'l');
% Turn down alanine
model = changeRxnBounds(model,'EX_cpd00035[e0]',0,'l');
% Turn down nitrogen
model = changeRxnBounds(model,'EX_cpd00528[e0]',0,'l');

% Turn off FWD too, so we can't grow on CO2
model = changeRxnBounds(model,'rxn11938[c0]',0,'b');

% Add in the H2 -> Fd_red reactions 
model = addReaction(model,{'rxn05759[c0]','Reduced ferredoxin:H+ oxidoreductase'},...
    '2 cpd00067[c0] + cpd11620[c0] <=> cpd11640[c0] + cpd11621[c0]');

% Add a row to freeEnergy to make dimensions correct
[~,idx] = intersect(model.rxns,'rxn05759[c0]');
model.freeEnergy(idx) = 0;

% Remove the bounds on Eha
model = removeEhaBounds(model);

% Rather than setting any substrate reactions, take in the generic model
% free energy and calculate the standard free energies based on fluxes

% Simulate growth to get flux distributions
solution = optimizeCbModel(model,[],'one');

% Calculate standard free energy
% First, find all the exchange reaction indices
exc_idx = findExcRxns(model);

% Use it to multiply the fluxes by free energies and get standard dG
gibbs_0 = sum(solution.x(exc_idx).*model.freeEnergy(exc_idx));

% Set standard values
R = 8.314e-6; %kJ/mmol*K
T = 310; %K

% Create a varied array of Q
qs = logspace(-200,0,200);

dGs = gibbs_0 + R*T*log(qs);

% Find value of Q that crosses 0
q_0 = exp(-gibbs_0/(R*T));


% Plot a figure showing Q vs dG
figure
semilogx(qs,dGs)
xlabel('Equilibrium Quotient (Q)')
ylabel('Free Energy Produced (kJ/gDCW/h)')
hold on
plot(q_0,0,'ro')


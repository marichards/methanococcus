function [design_2,solution1,design_3,solution2] = tryManualDesigns(model)

% Test out the manual strain designs from Dr. John Leigh


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

% Change bounds to function in reverse
% Change bounds such that methane goes IN instead of OUT
model = changeRxnBounds(model,'EX_cpd01024[e0]',-10,'l');
model = changeRxnBounds(model,'EX_cpd01024[e0]',0,'u');
% Change methanol to come OUT instead of IN
model = changeRxnBounds(model,'EX_cpd00116[e0]',0,'l');
model = changeRxnBounds(model,'EX_cpd00116[e0]',1000,'u');
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

% Let Eha/Ehb run at full blast
model = removeEhaBounds(model);

%%%%%%%%%%%
% All models will have the above things in common, now come the differences

% Design#1: Add in a reduction pathway that uses H2 directly; need to clone
% in a hydrogenase too


% Design#2: Add in H2 -> Fd reaction plus a reduction pathway (e.g. nitrate)
% then require that free energy be negative

design_2 = addReaction(model,{'rxn05759[c0]','Reduced ferredoxin:H+ oxidoreductase'},...
    '2 cpd00067[c0] + cpd11620[c0] <=> cpd11640[c0] + cpd11621[c0]');

% Turn off Fwd, which will ruin everything if left on
design_2 = changeRxnBounds(design_2,'rxn11938[c0]',0,'b');

% This completes the metabolite picture, but not the thermodynamic picture
% Do so using nitrate/nitrite

% Add outlet for nitrite (nitrate is already there)
design_2 = addReaction(design_2,{'rxn05626[c0]',...
    'Nitrite transport out via proton antiport'},...
    'cpd00067[e0] + cpd00075[c0] <=> cpd00067[c0] + cpd00075[e0]');
design_2 = addReaction(design_2,{'EX_cpd00075[e0]','EX_Nitrite[e0]'},...
    'cpd00075[e0] <=> ');
% Add metabolite info for Nitrite
[~,idx] = intersect(design_2.mets,'cpd00075[c0]');
design_2.metNames{idx} = 'Nitrite[c0]';
design_2.metCharge(idx)=-1;
design_2.metFormulas{idx}='NO2';
[~,idx] = intersect(design_2.mets,'cpd00075[e0]');
design_2.metNames{idx} = 'Nitrite[e0]';
design_2.metCharge(idx)=-1;
design_2.metFormulas{idx}='NO2';


% Add it to the free energy
[~,idx] = intersect(design_2.rxns,'EX_cpd00075[e0]');
design_2.freeEnergy(idx) = -0.0503;

% Add a Nitrate-> Nitrite reduction
% Ferredoxin version
design_2 = addReaction(design_2,{'rxn05894[c0]',...
    'Nitrite:ferredoxin oxidoreductase'},...
    'cpd00001[c0] + cpd00075[c0] + 2 cpd11621[c0] <=> 2 cpd00067[c0] + cpd00209[c0] + 2 cpd11620[c0]');

% NAD version
design_2 = addReaction(design_2,{'rxn00571[c0]',...
    'NADH:nitrate oxidoreductase'},...
    'cpd00001[c0] + cpd00075[c0] + cpd00003[c0] <=> cpd00067[c0] + cpd00209[c0] + cpd00004[c0]');

% NADP version
design_2 = addReaction(design_2,{'rxn00572[c0]',...
    'NADPH:nitrate oxidoreductase'},...
    'cpd00001[c0] + cpd00075[c0] + cpd00005[c0] <=> cpd00067[c0] + cpd00209[c0] + cpd00006[c0]');

% Non-mass balance version
design_2 = addReaction(design_2,'madeup','cpd00075[c0] <=> cpd00209[c0]');

% Add free energy for these reactions too (it's 0)
[~,idx] = intersect(design_2.rxns,'rxn05894[c0]');
design_2.freeEnergy(idx) = 0;
[~,idx] = intersect(design_2.rxns,'rxn00571[c0]');
design_2.freeEnergy(idx) = 0;
[~,idx] = intersect(design_2.rxns,'rxn00572[c0]');
design_2.freeEnergy(idx) = 0;
[~,idx] = intersect(design_2.rxns,'madeup');
design_2.freeEnergy(idx) = 0;

% Allow nitrate uptake/secretion infinitely
design_2 = changeRxnBounds(design_2,'EX_cpd00209[c0]',-inf,'l');
design_2 = changeRxnBounds(design_2,'EX_cpd00209[c0]',inf,'u');
% Also allow nitrite uptake/secretion infinitely
design_2 = changeRxnBounds(design_2,'EX_cpd00075[c0]',-inf,'l');
design_2 = changeRxnBounds(design_2,'EX_cpd00075[c0]',inf,'u');


% Specify substrate reactions and concentrations as 1 mM
substrate_rxns = {'EX_cpd00116[e0]','EX_cpd11640[e0]','EX_cpd01024[e0]',...
    'EX_cpd00075[e0]','EX_cpd00209[e0]'};
concentrations = [1 1 1 1 1];

% Solve by maximizing biomass, with negative free energy not allowed
[solution1,gibbs_flux,design_2] = optimizeThermoModel(design_2,substrate_rxns...
    ,concentrations,310,'EX_cpd00001[e0]',true);

% Find the reaction indices
[~,h2_idx]  = intersect(design_2.rxns,'EX_cpd11640[e0]');
[~,meoh_idx] = intersect(design_2.rxns,'EX_cpd00116[e0]');
[~,co2_idx] = intersect(design_2.rxns,'EX_cpd00011[e0]');
[~,ch4_idx] = intersect(design_2.rxns,'EX_cpd01024[e0]');
[~,h2o_idx] = intersect(design_2.rxns,'EX_cpd00001[e0]');
[~,nh3_idx] = intersect(design_2.rxns,'EX_cpd00013[e0]');
[~,po4_idx] = intersect(design_2.rxns,'EX_cpd00009[e0]');
[~,no3_idx] = intersect(design_2.rxns,'EX_cpd00209[e0]');
[~,no2_idx] = intersect(design_2.rxns,'EX_cpd00075[e0]');

if solution1.f > 0 
% Print the biomass flux
fprintf('\n\nBiomass flux: %f\n\n',solution1.f);

% Print the reaction fluxes
fprintf('Methanol flux: %f\n',solution1.x(meoh_idx))
fprintf('CO2 flux: %f\n',solution1.x(co2_idx))
fprintf('H2 flux: %f\n',solution1.x(h2_idx))
fprintf('H2O flux: %f\n',solution1.x(h2o_idx))
fprintf('CH4 flux: %f\n',solution1.x(ch4_idx))
fprintf('NH3 flux: %f\n',solution1.x(nh3_idx))
fprintf('PO4 flux: %f\n',solution1.x(po4_idx))
fprintf('NO3 flux: %f\n',solution1.x(no3_idx))
fprintf('NO2 flux: %f\n',solution1.x(no2_idx))

% Print the Gibbs flux
fprintf('\nPredicted Free Energy Generation: %f kJ/gDCW\n\n',gibbs_flux)
end


% Design#3: Use Fumarate --> Succinate conversion to reform disulfide; has
% no way to account for fd though! [Rxn is 'Tfr' and goes correct way
% already]

% Turn off Fwd, which will ruin everything if left on
design_3 = changeRxnBounds(model,'rxn11938[c0]',0,'b');

% Try putting in the ferredoxin->H2 conversion again (need a source of
% reduced ferredoxin for energy conservation)
design_3 = addReaction(design_3,{'rxn05759[c0]','Reduced ferredoxin:H+ oxidoreductase'},...
    '2 cpd00067[c0] + cpd11620[c0] <=> cpd11640[c0] + cpd11621[c0]');



% Add uptake/transport for fumarate and succinate
% First fumarate
design_3 = addReaction(design_3,{'EX_cpd00106[e0]','EX_Fumarate[e0]'},...
    'cpd00106[e0] <=> ');
design_3 = addReaction(design_3,{'rxn11013[c0]','EX_fum_e'},...
    'cpd00106[e0] <=> cpd00106[c0]');
% Then succinate
design_3 = addReaction(design_3,{'EX_cpd00036[e0]','EX_Succinate[e0]'},...
    'cpd00036[e0] <=> ');
design_3 = addReaction(design_3,{'rxn10952[c0]','EX_succ_e'},...
    'cpd00036[e0] <=> cpd00036[c0]');

% Add metabolite info for Both external metabolites
[~,idx] = intersect(design_3.mets,'cpd00106[e0]');
design_3.metNames{idx} = 'Fumarate[e0]';
design_3.metCharge(idx)=-2;
design_3.metFormulas{idx}='C4H2O4';
[~,idx] = intersect(design_3.mets,'cpd00036[e0]');
design_3.metNames{idx} = 'Succinate[e0]';
design_3.metCharge(idx)=-1;
design_3.metFormulas{idx}='C4H4O4';

% Add free energy for these reactions too (it's 0 for transports)
[~,idx] = intersect(design_3.rxns,'EX_cpd00106[e0]');
design_3.freeEnergy(idx) = -0.5367;
[~,idx] = intersect(design_3.rxns,'rxn11013[c0]');
design_3.freeEnergy(idx) = 0;
[~,idx] = intersect(design_3.rxns,'EX_cpd00036[e0]');
design_3.freeEnergy(idx) = -0.5358;
[~,idx] = intersect(design_3.rxns,'rxn10952[c0]');
design_3.freeEnergy(idx) = 0;

% Allow fumarate uptake/secretion infinitely
design_3 = changeRxnBounds(design_3,'EX_cpd00106[c0]',-inf,'l');
design_3 = changeRxnBounds(design_3,'EX_cpd00106[c0]',inf,'u');
% Also allow succinate uptake/secretion infinitely
design_3 = changeRxnBounds(design_3,'EX_cpd00036[c0]',-inf,'l');
design_3 = changeRxnBounds(design_3,'EX_cpd00036[c0]',inf,'u');

% Make sure Tfr is reversible
design_3 = changeRxnBounds(design_3,'Tfr',-1000,'l');

% Specify substrate reactions and concentrations as 1 mM
substrate_rxns = {'EX_cpd00116[e0]','EX_cpd11640[e0]','EX_cpd01024[e0]',...
    'EX_cpd00106[e0]','EX_cpd00036[e0]'};
concentrations = [1 1 1 1 1];

% Turn off Hdr's
%design_3 = changeRxnBounds(design_3,{'HdrABC','Hdr_formate'},0,'b');

% Solve by maximizing biomass, with negative free energy not allowed
[solution2,gibbs_flux,design_3] = optimizeThermoModel(design_3,substrate_rxns...
    ,concentrations,310,'EX_cpd00001[e0]',true);

% Find the reaction indices
[~,h2_idx]  = intersect(design_3.rxns,'EX_cpd11640[e0]');
[~,meoh_idx] = intersect(design_3.rxns,'EX_cpd00116[e0]');
[~,co2_idx] = intersect(design_3.rxns,'EX_cpd00011[e0]');
[~,ch4_idx] = intersect(design_3.rxns,'EX_cpd01024[e0]');
[~,h2o_idx] = intersect(design_3.rxns,'EX_cpd00001[e0]');
[~,nh3_idx] = intersect(design_3.rxns,'EX_cpd00013[e0]');
[~,po4_idx] = intersect(design_3.rxns,'EX_cpd00009[e0]');
[~,succ_idx] = intersect(design_3.rxns,'EX_cpd00036[e0]');
[~,fum_idx] = intersect(design_3.rxns,'EX_cpd00106[e0]');

if solution2.f > 0 
% Print the biomass flux
fprintf('\n\nBiomass flux: %f\n\n',solution2.f);

% Print the reaction fluxes
fprintf('Methanol flux: %f\n',solution2.x(meoh_idx))
fprintf('CO2 flux: %f\n',solution2.x(co2_idx))
fprintf('H2 flux: %f\n',solution2.x(h2_idx))
fprintf('H2O flux: %f\n',solution2.x(h2o_idx))
fprintf('CH4 flux: %f\n',solution2.x(ch4_idx))
fprintf('NH3 flux: %f\n',solution2.x(nh3_idx))
fprintf('PO4 flux: %f\n',solution2.x(po4_idx))
fprintf('Succinate flux: %f\n',solution2.x(succ_idx))
fprintf('Fumarate flux: %f\n',solution2.x(fum_idx))

% Print the Gibbs flux
fprintf('\nPredicted Free Energy Generation: %f kJ/gDCW\n\n',gibbs_flux)
end



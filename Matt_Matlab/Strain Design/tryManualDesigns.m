function [design_1,design_2,design_3,design_4] = tryManualDesigns(model)

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

% Add Fd <=> H2 Redox reaction
model = addReaction(model,{'rxn05759[c0]','Reduced ferredoxin:H+ oxidoreductase'},...
    '2 cpd00067[c0] + cpd11620[c0] <=> cpd11640[c0] + cpd11621[c0]');

% Turn off Fwd, which will ruin everything if left on by allowing CO2 to be
% metabolized
model = changeRxnBounds(model,'rxn11938[c0]',0,'b');

%%%%%%%%%%%
% All models will have the above things in common, now come the differences

% Design#1: Sulfate Reduction
design_1 = reduceSulfate(model);

% Design#2: Nitrate Reduction
design_2 = reduceNitrate(model);

% Design#3: Manganese Reduction
design_3 = reduceManganese(model);

% Design#4: Iron Reduction
design_4 = reduceIron(model);




end

function reduceSulfate(model)

end

function reduceNitrate(model)
% Design#2: Add in H2 -> Fd reaction plus a reduction pathway (e.g. nitrate)
% then require that free energy be negative

model = addReaction(model,{'rxn05759[c0]','Reduced ferredoxin:H+ oxidoreductase'},...
    '2 cpd00067[c0] + cpd11620[c0] <=> cpd11640[c0] + cpd11621[c0]');

% Turn off Fwd, which will ruin everything if left on
model = changeRxnBounds(model,'rxn11938[c0]',0,'b');

% This completes the metabolite picture, but not the thermodynamic picture
% Do so using nitrate/nitrite

% Add outlet for nitrite (nitrate is already there)
model = addReaction(model,{'rxn05626[c0]',...
    'Nitrite transport out via proton antiport'},...
    'cpd00067[e0] + cpd00075[c0] <=> cpd00067[c0] + cpd00075[e0]');
model = addReaction(model,{'EX_cpd00075[e0]','EX_Nitrite[e0]'},...
    'cpd00075[e0] <=> ');
% Add metabolite info for Nitrite
[~,idx] = intersect(model.mets,'cpd00075[c0]');
model.metNames{idx} = 'Nitrite[c0]';
model.metCharge(idx)=-1;
model.metFormulas{idx}='NO2';
[~,idx] = intersect(model.mets,'cpd00075[e0]');
model.metNames{idx} = 'Nitrite[e0]';
model.metCharge(idx)=-1;
model.metFormulas{idx}='NO2';

% Add it to the free energy
[~,idx] = intersect(model.rxns,'EX_cpd00075[e0]');
model.freeEnergy(idx) = -0.0503;

% Add a Nitrate-> Nitrite reduction
% Ferredoxin version
model = addReaction(model,{'rxn05894[c0]',...
    'Nitrite:ferredoxin oxidoreductase'},...
    'cpd00001[c0] + cpd00075[c0] + 2 cpd11621[c0] <=> 2 cpd00067[c0] + cpd00209[c0] + 2 cpd11620[c0]');

% NAD version
%design_2 = addReaction(design_2,{'rxn00571[c0]',...
%    'NADH:nitrate oxidoreductase'},...
%    'cpd00001[c0] + cpd00075[c0] + cpd00003[c0] <=> cpd00067[c0] + cpd00209[c0] + cpd00004[c0]');

% NADP version
%design_2 = addReaction(design_2,{'rxn00572[c0]',...
%    'NADPH:nitrate oxidoreductase'},...
%    'cpd00001[c0] + cpd00075[c0] + cpd00005[c0] <=> cpd00067[c0] + cpd00209[c0] + cpd00006[c0]');

% Add free energy for these reactions too (it's 0)
[~,idx] = intersect(model.rxns,'rxn05894[c0]');
model.freeEnergy(idx) = 0;
% [~,idx] = intersect(model.rxns,'rxn00571[c0]');
% model.freeEnergy(idx) = 0;
% [~,idx] = intersect(model.rxns,'rxn00572[c0]');
% model.freeEnergy(idx) = 0;
% [~,idx] = intersect(model.rxns,'madeup');
% model.freeEnergy(idx) = 0;

% Allow nitrate uptake/secretion infinitely
model = changeRxnBounds(model,'EX_cpd00209[c0]',-inf,'l');
model = changeRxnBounds(model,'EX_cpd00209[c0]',inf,'u');
% Also allow nitrite uptake/secretion infinitely
model = changeRxnBounds(model,'EX_cpd00075[c0]',-inf,'l');
model = changeRxnBounds(model,'EX_cpd00075[c0]',inf,'u');
% Do the same for reaction(s) I just added
model = changeRxnBounds(model,'rxn05894[c0]',-inf,'l');
model = changeRxnBounds(model,'rxn05894[c0]',inf,'u');
%design_2 = changeRxnBounds(design_2,'rxn00571[c0]',-inf,'l');
%design_2 = changeRxnBounds(design_2,'rxn00571[c0]',inf,'u');
%design_2 = changeRxnBounds(design_2,'rxn00572[c0]',-inf,'l');
%design_2 = changeRxnBounds(design_2,'rxn00572[c0]',inf,'u');

% Specify substrate reactions and concentrations as 1 mM
substrate_rxns = {'EX_cpd00116[e0]','EX_cpd11640[e0]','EX_cpd01024[e0]',...
    'EX_cpd00075[e0]','EX_cpd00209[e0]'};
concentrations = [1 1 1 1 1];

% Solve by maximizing biomass, with negative free energy not allowed
[solution1,gibbs_flux,model] = optimizeThermoModel(model,substrate_rxns...
    ,concentrations,310,'EX_cpd00001[e0]',true);

% Find the reaction indices
[~,h2_idx]  = intersect(model.rxns,'EX_cpd11640[e0]');
[~,meoh_idx] = intersect(model.rxns,'EX_cpd00116[e0]');
[~,co2_idx] = intersect(model.rxns,'EX_cpd00011[e0]');
[~,ch4_idx] = intersect(model.rxns,'EX_cpd01024[e0]');
[~,h2o_idx] = intersect(model.rxns,'EX_cpd00001[e0]');
[~,nh3_idx] = intersect(model.rxns,'EX_cpd00013[e0]');
[~,po4_idx] = intersect(model.rxns,'EX_cpd00009[e0]');
[~,no3_idx] = intersect(model.rxns,'EX_cpd00209[e0]');
[~,no2_idx] = intersect(model.rxns,'EX_cpd00075[e0]');

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
end


function reduceManganese(model)

end

function reduceIron(model)

end

% Design#3: Use Fumarate --> Succinate conversion to reform disulfide; has
% no way to account for fd though! [Rxn is 'Tfr' and goes correct way
% already]

% Turn off Fwd, which will ruin everything if left on
%model = changeRxnBounds(model,'rxn11938[c0]',0,'b');

% Try putting in the ferredoxin->H2 conversion again (need a source of
% reduced ferredoxin for energy conservation)
%model = addReaction(model,{'rxn05759[c0]','Reduced ferredoxin:H+ oxidoreductase'},...
%     '2 cpd00067[c0] + cpd11620[c0] <=> cpd11640[c0] + cpd11621[c0]');



% Add uptake/transport for fumarate and succinate
% First fumarate
%model = addReaction(model,{'EX_cpd00106[e0]','EX_Fumarate[e0]'},...
%    'cpd00106[e0] <=> ');
%model = addReaction(model,{'rxn11013[c0]','EX_fum_e'},...
%    'cpd00106[e0] <=> cpd00106[c0]');
% Then succinate
%model = addReaction(model,{'EX_cpd00036[e0]','EX_Succinate[e0]'},...
%    'cpd00036[e0] <=> ');
%model = addReaction(model,{'rxn10952[c0]','EX_succ_e'},...
%    'cpd00036[e0] <=> cpd00036[c0]');

% Add metabolite info for Both external metabolites
%[~,idx] = intersect(model.mets,'cpd00106[e0]');
%model.metNames{idx} = 'Fumarate[e0]';
%model.metCharge(idx)=-2;
%model.metFormulas{idx}='C4H2O4';
%[~,idx] = intersect(model.mets,'cpd00036[e0]');
%model.metNames{idx} = 'Succinate[e0]';
%model.metCharge(idx)=-1;
%model.metFormulas{idx}='C4H4O4';

% Add free energy for these reactions too (it's 0 for transports)
% [~,idx] = intersect(model.rxns,'EX_cpd00106[e0]');
% model.freeEnergy(idx) = -0.5367;
% [~,idx] = intersect(model.rxns,'rxn11013[c0]');
% model.freeEnergy(idx) = 0;
% [~,idx] = intersect(model.rxns,'EX_cpd00036[e0]');
% model.freeEnergy(idx) = -0.5358;
% [~,idx] = intersect(model.rxns,'rxn10952[c0]');
% model.freeEnergy(idx) = 0;

% % Allow fumarate uptake/secretion infinitely
% model = changeRxnBounds(model,'EX_cpd00106[c0]',-inf,'l');
% model = changeRxnBounds(model,'EX_cpd00106[c0]',inf,'u');
% % Also allow succinate uptake/secretion infinitely
% model = changeRxnBounds(model,'EX_cpd00036[c0]',-inf,'l');
% model = changeRxnBounds(model,'EX_cpd00036[c0]',inf,'u');
% 
% % Make sure Tfr is reversible
% model = changeRxnBounds(model,'Tfr',-1000,'l');
% 
% % Specify substrate reactions and concentrations as 1 mM
% substrate_rxns = {'EX_cpd00116[e0]','EX_cpd11640[e0]','EX_cpd01024[e0]',...
%     'EX_cpd00106[e0]','EX_cpd00036[e0]'};
% concentrations = [1 1 1 1 1];
% 
% % Turn off Hdr's
% %design_3 = changeRxnBounds(design_3,{'HdrABC','Hdr_formate'},0,'b');
% 
% % Solve by maximizing biomass, with negative free energy not allowed
% [solution2,gibbs_flux,model] = optimizeThermoModel(model,substrate_rxns...
%     ,concentrations,310,'EX_cpd00001[e0]',true);
% 
% % Find the reaction indices
% [~,h2_idx]  = intersect(model.rxns,'EX_cpd11640[e0]');
% [~,meoh_idx] = intersect(model.rxns,'EX_cpd00116[e0]');
% [~,co2_idx] = intersect(model.rxns,'EX_cpd00011[e0]');
% [~,ch4_idx] = intersect(model.rxns,'EX_cpd01024[e0]');
% [~,h2o_idx] = intersect(model.rxns,'EX_cpd00001[e0]');
% [~,nh3_idx] = intersect(model.rxns,'EX_cpd00013[e0]');
% [~,po4_idx] = intersect(model.rxns,'EX_cpd00009[e0]');
% [~,succ_idx] = intersect(model.rxns,'EX_cpd00036[e0]');
% [~,fum_idx] = intersect(model.rxns,'EX_cpd00106[e0]');
% 
% if solution2.f > 0 
% % Print the biomass flux
% fprintf('\n\nBiomass flux: %f\n\n',solution2.f);
% 
% % Print the reaction fluxes
% fprintf('Methanol flux: %f\n',solution2.x(meoh_idx))
% fprintf('CO2 flux: %f\n',solution2.x(co2_idx))
% fprintf('H2 flux: %f\n',solution2.x(h2_idx))
% fprintf('H2O flux: %f\n',solution2.x(h2o_idx))
% fprintf('CH4 flux: %f\n',solution2.x(ch4_idx))
% fprintf('NH3 flux: %f\n',solution2.x(nh3_idx))
% fprintf('PO4 flux: %f\n',solution2.x(po4_idx))
% fprintf('Succinate flux: %f\n',solution2.x(succ_idx))
% fprintf('Fumarate flux: %f\n',solution2.x(fum_idx))
% 
% % Print the Gibbs flux
% fprintf('\nPredicted Free Energy Generation: %f kJ/gDCW\n\n',gibbs_flux)
% end



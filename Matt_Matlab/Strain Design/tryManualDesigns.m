function [solution_1,solution_2,solution_3,solution_4,...
    gibbs_flux_1,gibbs_flux_2,gibbs_flux_3,gibbs_flux_4,...
    model_1,model_2,model_3,model_4] = tryManualDesigns(model)

% Test out the manual strain designs from Dr. John Leigh


% First add in methanol pathway
model = addReverseMethanolPathway(model);

% Add Fd <=> H2 Redox reaction
model = addReaction(model,{'rxn05759[c0]','Reduced ferredoxin:H+ oxidoreductase'},...
    '2 cpd00067[c0] + cpd11620[c0] <=> cpd11640[c0] + cpd11621[c0]');

% Turn off Fwd, which will ruin everything if left on by allowing CO2 to be
% metabolized
model = changeRxnBounds(model,'rxn11938[c0]',0,'b');

%%%%%%%%%%%
% All models will have the above things in common, now come the differences

% Simulate all models
[solution_1,gibbs_flux_1,model_1] = reduceSulfate(model);
[solution_2,gibbs_flux_2,model_2] = reduceNitrate(model);
[solution_3,gibbs_flux_3,model_3] = reduceManganese(model);
[solution_4,gibbs_flux_4,model_4] = reduceIron(model);

end

function [solution,gibbs_flux,model] = reduceSulfate(model)

% Reactions from original model:
% EX_cpd00048_e0	Sulfate_e0 	<=>		
% rxn05651_c0	H_e0 + Sulfate_e0 	<=>	H_c0 + Sulfate_c0 
% rxn00379_c0	H_c0 + ATP_c0 + Sulfate_c0 	->	PPi_c0 + APS_c0
% rxn05256_c0	trdrd_c0 + APS_c0 	->	trdox_c0 + Sulfite_c0 + AMP_c0
% rxn00623_c0	4.000000 H_c0 + Sulfite_c0 + 3.000000 NADPH_c0 	->	3.000000 H2O_c0 + 3.000000 NADP_c0 + H2S_c0

% Compound IDs:
% Sulfate: cpd00048
% APS: cpd00193
% Sulfite: cpd00081
% H2S: cpd00239
% trdrd: cpd11421
% trdox: cpd11420
% NADPH: cpd00005
% NADP: cpd00006
% PPi: cpd00012
% 

% Add the uptake and inlet for sulfate
model = addReaction(model,{'EX_cpd00048[e0]','EX_Sulfate_e0'},...
    'cpd00048[e0] <=> ');
model = addReaction(model,{'rxn05651[c0]','Sulfate transport via proton symport'},...
    'cpd00048[e0] + cpd00067[e0] <=> cpd00048[c0] + cpd00067[c0]');
% Add the Sulfate --> APS pathway
model = addReaction(model,{'rxn00379[c0]','Sulfate reductase'},...
    'cpd00067[c0] + cpd00002[c0] + cpd00048[c0] -> cpd00012[c0] + cpd00193[c0]');
% Add the APS --> Sulfite pathway
model = addReaction(model,{'rxn05256[c0]','APS reductase'},...
    'cpd11421[c0] + cpd00193[c0] -> cpd11420[c0] + cpd00081[c0] + cpd00018[c0]');

% Add free energies for all reactions (0 for all but the exchange)
% Add it to the free energy
[~,idx] = intersect(model.rxns,'EX_cpd00048[e0]');
model.freeEnergy(idx) = -0.7649;
[~,idx] = intersect(model.rxns,'rxn05651[c0]');
model.freeEnergy(idx) = 0;
[~,idx] = intersect(model.rxns,'rxn00379[c0]');
model.freeEnergy(idx) = 0;
[~,idx] = intersect(model.rxns,'rxn05256[c0]');
model.freeEnergy(idx) = 0;

% Allow sulfate uptake/secretion infinitely
model = changeRxnBounds(model,'EX_cpd00048[e0]',-inf,'l');
model = changeRxnBounds(model,'EX_cpd00048[e0]',inf,'u');

% OPTION: Add one, comment the other

% Option 1, add Sulfite -> H2S
model = addReaction(model,{'rxn00623[c0]','Sulfite reductase'},...
    '4 cpd00067[c0] + cpd00081[c0] + 3 cpd00005[c0] -> 3 cpd00001[c0] + 3 cpd00006[c0] + cpd00239[c0]');
[~,idx] = intersect(model.rxns,'rxn00623[c0]');
model.freeEnergy(idx) = 0;
% Also allow sulfide uptake/secretion infinitely
model = changeRxnBounds(model,'EX_cpd00239[e0]',-inf,'l');
model = changeRxnBounds(model,'EX_cpd00239[e0]',inf,'u');


% % Option 2, add Sulfite outlet
% model = addReaction(model,{'rxn10604[c0]','Sulfate transport exchange via proton symport'},...
%     'cpd00067[e0] + cpd00081[e0] <=> cpd00067[c0] + cpd00081[c0]');
% model = addReaction(model,{'EX_cpd00081[e0]','EX_Sulfite_e0'},...
%     'cpd00081[e0] <=> ');
% % Add free energies
% [~,idx] = intersect(model.rxns,'EX_cpd00048[e0]');
% model.freeEnergy(idx) = -0.5046;
% [~,idx] = intersect(model.rxns,'rxn10604[c0]');
% model.freeEnergy(idx) = 0;
% % Also allow sulfite uptake/secretion infinitely
% model = changeRxnBounds(model,'EX_cpd00081[e0]',10,'l');
% model = changeRxnBounds(model,'EX_cpd00081[e0]',10,'u');


% Specify substrate reactions and concentrations as 1 mM
substrate_rxns = {'EX_cpd00116[e0]','EX_cpd11640[e0]','EX_cpd01024[e0]',...%,...
    'EX_cpd00048[e0]','EX_cpd00239[e0]'};
concentrations = [1 1 1 1e10 0.001];

% Solve by maximizing biomass, with negative free energy not allowed
[solution,gibbs_flux,model] = optimizeThermoModel(model,substrate_rxns...
    ,concentrations,310,'EX_cpd00001[e0]',true);

% Find the reaction indices
[~,h2_idx]  = intersect(model.rxns,'EX_cpd11640[e0]');
[~,meoh_idx] = intersect(model.rxns,'EX_cpd00116[e0]');
[~,co2_idx] = intersect(model.rxns,'EX_cpd00011[e0]');
[~,ch4_idx] = intersect(model.rxns,'EX_cpd01024[e0]');
[~,h2o_idx] = intersect(model.rxns,'EX_cpd00001[e0]');
[~,nh3_idx] = intersect(model.rxns,'EX_cpd00013[e0]');
[~,po4_idx] = intersect(model.rxns,'EX_cpd00009[e0]');
[~,h2s_idx] = intersect(model.rxns,'EX_cpd00239[e0]');
[~,so4_idx] = intersect(model.rxns,'EX_cpd00048[e0]');


if solution.f > 0 
% Print the biomass flux
fprintf('\n\nBiomass flux: %f\n\n',solution.f);

% Print the reaction fluxes
fprintf('Methanol flux: %f\n',solution.x(meoh_idx))
fprintf('CO2 flux: %f\n',solution.x(co2_idx))
fprintf('H2 flux: %f\n',solution.x(h2_idx))
fprintf('H2O flux: %f\n',solution.x(h2o_idx))
fprintf('CH4 flux: %f\n',solution.x(ch4_idx))
fprintf('NH3 flux: %f\n',solution.x(nh3_idx))
fprintf('PO4 flux: %f\n',solution.x(po4_idx))
fprintf('H2S flux: %f\n',solution.x(h2s_idx))
fprintf('SO4 flux: %f\n',solution.x(so4_idx))

% Print the Overall Reaction
fprintf('Predicted Overall Reaction: CH4 + %0.2f CO2 + %0.2f SO4 --> %0.2f CH3OH + %0.2f H2S+ %0.2f H2O',...
    (solution.x(co2_idx)/solution.x(ch4_idx)),(solution.x(so4_idx)/solution.x(ch4_idx)),...
    -(solution.x(meoh_idx)/solution.x(ch4_idx)),-(solution.x(h2s_idx)/solution.x(ch4_idx)),...
    -(solution.x(h2o_idx)/solution.x(ch4_idx)));

% Print the Gibbs flux
fprintf('\nPredicted Free Energy Generation: %f kJ/gDCW\n\n',gibbs_flux)
end

end

function [solution,gibbs_flux,model] = reduceNitrate(model)
% Design#2: Add in H2 -> Fd reaction plus a reduction pathway (e.g. nitrate)
% then require that free energy be negative

% Add outlet for nitrite (nitrate is already there)
model = addReaction(model,{'rxn05626[c0]',...
    'Nitrite transport out via proton antiport'},...
    'cpd00067[e0] + cpd00075[c0] <=> cpd00067[c0] + cpd00075[e0]');
% Try adding ABC-transporter
% model = addReaction(model,'test_rxn',...
% 'cpd00001[c0] + cpd00002[c0] + cpd00075[e0] <=> cpd00067[c0] + cpd00009[c0] + cpd00008[c0] + cpd00075[c0]');

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

% Allow nitrate uptake/secretion infinitely
model = changeRxnBounds(model,'EX_cpd00209[e0]',-inf,'l');
model = changeRxnBounds(model,'EX_cpd00209[e0]',inf,'u');
% Also allow nitrite uptake/secretion infinitely
model = changeRxnBounds(model,'EX_cpd00075[e0]',-inf,'l');
model = changeRxnBounds(model,'EX_cpd00075[e0]',inf,'u');

% Add a Nitrate-> Nitrite reduction
% Ferredoxin version
% model = addReaction(model,{'rxn05894[c0]',...
%     'Nitrite:ferredoxin oxidoreductase'},...
%     'cpd00001[c0] + cpd00075[c0] + 2 cpd11621[c0] <=> 2 cpd00067[c0] + cpd00209[c0] + 2 cpd11620[c0]');
% % Add free energy and bounds for these reactions too (it's 0 and inf/-inf)
% [~,idx] = intersect(model.rxns,'rxn05894[c0]');
% model.freeEnergy(idx) = 0;
% model = changeRxnBounds(model,'rxn05894[c0]',-inf,'l');
% model = changeRxnBounds(model,'rxn05894[c0]',inf,'u');

% NAD version
model = addReaction(model,{'rxn00571[c0]',...
   'NADH:nitrate oxidoreductase'},...
   'cpd00001[c0] + cpd00075[c0] + cpd00003[c0] <=> cpd00067[c0] + cpd00209[c0] + cpd00004[c0]');
% Add free energy and bounds for these reactions too (it's 0 and inf/-inf)
[~,idx] = intersect(model.rxns,'rxn00571[c0]');
model.freeEnergy(idx) = 0;
model = changeRxnBounds(model,'rxn00571[c0]',-inf,'l');
model = changeRxnBounds(model,'rxn00571[c0]',inf,'u');
% 
% % NADP version
% model = addReaction(model,{'rxn00572[c0]',...
%    'NADPH:nitrate oxidoreductase'},...
%    'cpd00001[c0] + cpd00075[c0] + cpd00006[c0] <=> cpd00067[c0] + cpd00209[c0] + cpd00005[c0]');
% % Add free energy and bounds for these reactions too (it's 0 and inf/-inf)
% [~,idx] = intersect(model.rxns,'rxn00572[c0]');
% model.freeEnergy(idx) = 0;
% model = changeRxnBounds(model,'rxn00572[c0]',-inf,'l');
% model = changeRxnBounds(model,'rxn00572[c0]',inf,'u');

% Made up version; NO3 + H2 -> NO2 + H2O
% model = addReaction(model,'madeup',...
%     'cpd00075[c0] <=> cpd00209[c0] + 2 cpd00001[c0]');
% [~,idx] = intersect(model.rxns,'madeup');
% model.freeEnergy(idx) = 0;
% model = changeRxnBounds(model,'madeup',-inf,'l');
% model = changeRxnBounds(model,'madeup',inf,'u');
% % Force flux
% model = changeRxnBounds(model,'madeup',10,'b');

% Add different nitrate transporter
model = addReaction(model,{'madeup2',...
    'Nitrate transport out via proton antiport'},...
    'cpd00067[e0] + cpd00209[c0] <=> cpd00067[c0] + cpd00209[e0]');

% Specify substrate reactions and concentrations as 1 mM
substrate_rxns = {'EX_cpd00116[e0]','EX_cpd11640[e0]','EX_cpd01024[e0]',...
    'EX_cpd00075[e0]','EX_cpd00209[e0]'};
concentrations = [1 1 1 1 1];

% Solve by maximizing biomass, with negative free energy not allowed
[solution,gibbs_flux,model] = optimizeThermoModel(model,substrate_rxns...
    ,concentrations,310,'EX_cpd00001[e0]',true);

% Find the reaction indices
[~,h2_idx]  = intersect(model.rxns,'EX_cpd11640[e0]');
[~,meoh_idx] = intersect(model.rxns,'EX_cpd00116[e0]');
[~,co2_idx] = intersect(model.rxns,'EX_cpd00011[e0]');
[~,ch4_idx] = intersect(model.rxns,'EX_cpd01024[e0]');
[~,h2o_idx] = intersect(model.rxns,'EX_cpd00001[e0]');
[~,nh3_idx] = intersect(model.rxns,'EX_cpd00013[e0]');
[~,po4_idx] = intersect(model.rxns,'EX_cpd00009[e0]');
[~,h2s_idx] = intersect(model.rxns,'EX_cpd00239[e0]');
[~,no3_idx] = intersect(model.rxns,'EX_cpd00209[e0]');
[~,no2_idx] = intersect(model.rxns,'EX_cpd00075[e0]');

if solution.f > 0 
% Print the biomass flux
fprintf('\n\nBiomass flux: %f\n\n',solution.f);

% Print the reaction fluxes
fprintf('Methanol flux: %f\n',solution.x(meoh_idx))
fprintf('CO2 flux: %f\n',solution.x(co2_idx))
fprintf('H2 flux: %f\n',solution.x(h2_idx))
fprintf('H2O flux: %f\n',solution.x(h2o_idx))
fprintf('CH4 flux: %f\n',solution.x(ch4_idx))
fprintf('NH3 flux: %f\n',solution.x(nh3_idx))
fprintf('PO4 flux: %f\n',solution.x(po4_idx))
fprintf('H2S flux: %f\n',solution.x(h2s_idx))
fprintf('NO3 flux: %f\n',solution.x(no3_idx))
fprintf('NO2 flux: %f\n',solution.x(no2_idx))

% Print the Overall Reaction
fprintf('Predicted Overall Reaction: CH4 + %0.2f CO2 + %0.2f NO3 --> %0.2f CH3OH + %0.2f NO2+ %0.2f H2O',...
    (solution.x(co2_idx)/solution.x(ch4_idx)),(solution.x(no3_idx)/solution.x(ch4_idx)),...
    -(solution.x(meoh_idx)/solution.x(ch4_idx)),-(solution.x(no2_idx)/solution.x(ch4_idx)),...
    -(solution.x(h2o_idx)/solution.x(ch4_idx)));

% Print the Gibbs flux
fprintf('\nPredicted Free Energy Generation: %f kJ/gDCW\n\n',gibbs_flux)
end

end


function [solution,gibbs_flux,model] = reduceManganese(model)

% Manganese already has an exchange, but MnO2 does not; give it one, plus a
% transporter for both and a redox reaction w/NAD

% Add exchange for MnO2, with free energy
model = addReaction(model,'EX_cpd17028[e0]','cpd17028[e0] <=> ');
[~,idx] = intersect(model.rxns,'EX_cpd17028[e0]');
model.freeEnergy(idx) = -0.0503;

% Specify substrate reactions and concentrations as 1 mM
substrate_rxns = {'EX_cpd00116[e0]','EX_cpd11640[e0]','EX_cpd01024[e0]'};
concentrations = [1 1 1];

% Solve by maximizing biomass, with negative free energy not allowed
[solution,gibbs_flux,model] = optimizeThermoModel(model,substrate_rxns...
    ,concentrations,310,'EX_cpd00001[e0]',true);

end

function [solution,gibbs_flux,model] = reduceIron(model)


% Iron forms are both already in the model; add an oxidoreductase

% Potential problems: ABC transporters, but both types of iron have the
% same type so we should be okay.

% % NAD version (note that there's no water produced here, so it's not quite
% % as favorable as I'd like. However, this is the only option I have
% % available
% model = addReaction(model,{'rxn00068[c0]',...
%    'NADH:Fe2+ oxidoreductase'},...
%     '2 cpd10515[c0] + cpd00067[c0] + cpd00003[c0] <=> 2 cpd10516[c0] + cpd00004[c0]');
% % Add free energy and bounds for these reactions too (it's 0 and inf/-inf)
% [~,idx] = intersect(model.rxns,'rxn00068[c0]');
% model.freeEnergy(idx) = 0;
% model = changeRxnBounds(model,'rxn00068[c0]',-inf,'l');
% model = changeRxnBounds(model,'rxn00068[c0]',inf,'u');

% % NAD version (note that there's no water produced here, so it's not quite
% % as favorable as I'd like. However, this is the only option I have
% % available
% model = addReaction(model,{'altered_nad',...
%    'NADH:Fe2+ oxidoreductase'},...
%     'cpd00067[c0] + cpd00001[c0] + 2 cpd10515[c0] + cpd00003[c0] <=>  2 cpd10516[c0] + cpd00004[c0]');
% % Add free energy and bounds for these reactions too (it's 0 and inf/-inf)
% [~,idx] = intersect(model.rxns,'altered_nad');
% model.freeEnergy(idx) = 0;
% model = changeRxnBounds(model,'altered_nad',-inf,'l');
% model = changeRxnBounds(model,'altered_nad',inf,'u');


% % Try adding an NAD Source/Sink
%  model = addReaction(model,'NAD_sink','cpd00003[c0] <=>');
% model = addReaction(model,'NADH_sink','cpd00004[c0] <=>');


% NADP version (not an actual Kbase reaction)
model = addReaction(model,{'madeup',...
   'NADPH:Fe2+ oxidoreductase'},...
    '2 cpd10515[c0] + cpd00067[c0] + cpd00006[c0] <=> 2 cpd10516[c0] + cpd00005[c0]');
% Add free energy and bounds for these reactions too (it's 0 and inf/-inf)
[~,idx] = intersect(model.rxns,'madeup');
model.freeEnergy(idx) = 0;
model = changeRxnBounds(model,'madeup',-inf,'l');
model = changeRxnBounds(model,'madeup',inf,'u');


% % % Try a madeup version
% model = addReaction(model,'madeup','cpd10515[e0] <=> cpd10516[e0]');
% % Add free energy and bounds for these reactions too (it's 0 and inf/-inf)
% [~,idx] = intersect(model.rxns,'madeup');
% model.freeEnergy(idx) = 0;
% model = changeRxnBounds(model,'madeup',-inf,'l');
% model = changeRxnBounds(model,'madeup',inf,'u');



% Allow iron in/out as much as necessary; also turn on the transporters to
% infinite bounds
% Allow Fe3+ uptake/secretion infinitely
model = changeRxnBounds(model,'EX_cpd10516[e0]',-inf,'l');
model = changeRxnBounds(model,'EX_cpd10516[e0]',inf,'u');
% Allow Fe2+ uptake/secretion infinitely
model = changeRxnBounds(model,'EX_cpd10515[e0]',-inf,'l');
model = changeRxnBounds(model,'EX_cpd10515[e0]',inf,'u');
% Allow Fe3+ transport infinitely
model = changeRxnBounds(model,'rxn05195[c0]',-inf,'l');
model = changeRxnBounds(model,'rxn05195[c0]',inf,'u');
% Allow Fe2+ transport infinitely
model = changeRxnBounds(model,'rxn05555[c0]',-inf,'l');
model = changeRxnBounds(model,'rxn05555[c0]',inf,'u');


% Specify substrate reactions and concentrations as 1 mM
substrate_rxns = {'EX_cpd00116[e0]','EX_cpd11640[e0]','EX_cpd01024[e0]',...
    'EX_cpd10515[e0]','EX_cpd10516[e0]'};
concentrations = [1 1 1 1 1];

% Solve by maximizing biomass, with negative free energy not allowed
[solution,gibbs_flux,model] = optimizeThermoModel(model,substrate_rxns...
    ,concentrations,310,'EX_cpd00001[e0]',true);

% Find the reaction indices
[~,h2_idx]  = intersect(model.rxns,'EX_cpd11640[e0]');
[~,meoh_idx] = intersect(model.rxns,'EX_cpd00116[e0]');
[~,co2_idx] = intersect(model.rxns,'EX_cpd00011[e0]');
[~,ch4_idx] = intersect(model.rxns,'EX_cpd01024[e0]');
[~,h2o_idx] = intersect(model.rxns,'EX_cpd00001[e0]');
[~,nh3_idx] = intersect(model.rxns,'EX_cpd00013[e0]');
[~,po4_idx] = intersect(model.rxns,'EX_cpd00009[e0]');
[~,h2s_idx] = intersect(model.rxns,'EX_cpd00239[e0]');
[~,fe3_idx] = intersect(model.rxns,'EX_cpd10516[e0]');
[~,fe2_idx] = intersect(model.rxns,'EX_cpd10515[e0]');

if solution.f > 0 
% Print the biomass flux
fprintf('\n\nBiomass flux: %f\n\n',solution.f);

% Print the reaction fluxes
fprintf('Methanol flux: %f\n',solution.x(meoh_idx))
fprintf('CO2 flux: %f\n',solution.x(co2_idx))
fprintf('H2 flux: %f\n',solution.x(h2_idx))
fprintf('H2O flux: %f\n',solution.x(h2o_idx))
fprintf('CH4 flux: %f\n',solution.x(ch4_idx))
fprintf('NH3 flux: %f\n',solution.x(nh3_idx))
fprintf('PO4 flux: %f\n',solution.x(po4_idx))
fprintf('H2S flux: %f\n',solution.x(h2s_idx))
fprintf('Fe3 flux: %f\n',solution.x(fe3_idx))
fprintf('Fe2 flux: %f\n',solution.x(fe2_idx))

% Print the Overall Reaction
fprintf('Predicted Overall Reaction: CH4 + %0.2f CO2 + %0.2f Fe3+ --> %0.2f CH3OH + %0.2f Fe2+ %0.2f H2O',...
    (solution.x(co2_idx)/solution.x(ch4_idx)),(solution.x(fe3_idx)/solution.x(ch4_idx)),...
    -(solution.x(meoh_idx)/solution.x(ch4_idx)),-(solution.x(fe2_idx)/solution.x(ch4_idx)),...
    -(solution.x(h2o_idx)/solution.x(ch4_idx)));

% Print the Gibbs flux
fprintf('\nPredicted Free Energy Generation: %f kJ/gDCW\n\n',gibbs_flux)
end

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

function model = addReverseMethanolPathway(model)

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
model = changeRxnBounds(model,'EX_cpd01024[e0]',-inf,'l');
model = changeRxnBounds(model,'EX_cpd01024[e0]',0,'u');
% Change methanol to come OUT instead of IN
model = changeRxnBounds(model,'EX_cpd00116[e0]',10,'b');

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

end

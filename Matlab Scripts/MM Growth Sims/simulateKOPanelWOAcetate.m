function simulateKOPanelWOAcetate(model)

% For the M. maripaludis S2 model, simulate the model for known gene KO
% experiments to get predictions and compare predictions to reality. Do all
% possible KOs for all 4 conditions, not just the replicates of
% experiments. In this specialized case, also remove acetate from the media
% to represent what we'd predict without it. 
%
% INPUT
% model: the M. maripaludis model, a COBRA Toolbox model structure
%
% Matthew Richards, 09/24/2015


% Alteration on 05/26/2015: Add MCC and accuracy calculations
% Create TP/TN/FP/FN metrics to fill up on appropriate things
tp = 0; tn = 0; fp = 0; fn = 0;

% Set GAPOR off to begin with
model = changeRxnBounds(model,'rxn07191[c0]',0,'b');
% Make sure model is set to H2
model = switchToH2Only(model);

% H2-CO2 simulations
fprintf('================================\nGrowth on H2 + CO2\n================================');

% First simulate Wild-type growth
solution = optimizeCbModel(model,[],'one');
fprintf('\nWild-Type Growth: %0.2f\n\n',solution.f);
wt_growth = solution.f;

% Simulate Hmd KO (MMP0127)
ko_model = deleteModelGenes(model,'MMP0127',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-Hmd Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (non-lethal)
if solution.f/wt_growth >= 0.1
    % Then add to True Positive
    tp = tp+1;
else
    % Then add to False Negative
    fn = fn+1;
end

% Simulate Mtd KO (MMP0372)
ko_model = deleteModelGenes(model,'MMP0372',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-Mtd Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (non-lethal)
if solution.f/wt_growth >= 0.1
    % Then add to True Positive
    tp = tp+1;
else
    % Then add to False Negative
    fn = fn+1;
end

% Simulate FrcA KO (MMP0820)
ko_model = deleteModelGenes(model,'MMP0820',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FrcA Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (non-lethal)
if solution.f/wt_growth >= 0.1
    % Then add to True Positive
    tp = tp+1;
else
    % Then add to False Negative
    fn = fn+1;
end

% Simulate FruA KO (MMP1382)
ko_model = deleteModelGenes(model,'MMP1382',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FruA Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (non-lethal)
if solution.f/wt_growth >= 0.1
    % Then add to True Positive
    tp = tp+1;
else
    % Then add to False Negative
    fn = fn+1;
end

% Simulate FrcA-FruA double KO (MMP0820 and MMP1382)
ko_model = deleteModelGenes(model,{'MMP0820','MMP1382'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FrcA-FruA Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (non-lethal)
if solution.f/wt_growth >= 0.1
    % Then add to True Positive
    tp = tp+1;
else
    % Then add to False Negative
    fn = fn+1;
end

% Simulate  VhuAU-VhcA triple KO (MMP1694, MMP1693, MMP0823)
ko_model = deleteModelGenes(model,{'MMP0680','MMP1694','MMP1693','MMP0823'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-VhuAU-VhcA Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (non-lethal)
if solution.f/wt_growth >= 0.1
    % Then add to True Positive
    tp = tp+1;
else
    % Then add to False Negative
    fn = fn+1;
end

% Simulate  HdrB2 KO (MMP1053)
ko_model = deleteModelGenes(model,{'MMP0680','MMP1053'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-HdrB2 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  FdhA1 KO (MMP1298)
ko_model = deleteModelGenes(model,{'MMP1298'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA1 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  FdhA2 KO (MMP0138)
ko_model = deleteModelGenes(model,{'MMP0138'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA2 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  FdhA1-FdhA2 double KO (MMP1298 and MMP0138)
ko_model = deleteModelGenes(model,{'MMP1298','MMP0138'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA1-FdhA2 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  FdhA2B2 KO (MMP0138 and MMP0139)
ko_model = deleteModelGenes(model,{'MMP0680','MMP0138','MMP0139'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA2B2 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate EhbF KO (MMP1628)
ko_model = deleteModelGenes(model,'MMP1628',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-EhbF Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (non-lethal)
if solution.f/wt_growth >= 0.1
    % Then add to True Positive
    tp = tp+1;
else
    % Then add to False Negative
    fn = fn+1;
end

% Simulate 3H2ase KOs of frcAGB,fruAGB,hmd
% (MMP0820, MMP0818, MMP817, MMP1382, MMP1384, MMP1385, and MMP0127)
ko_model = deleteModelGenes(model,...
    {'MMP0680','MMP0820','MMP0818','MMP0817','MMP1382','MMP1384','MMP1385','MMP0127'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-3H2ase Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (non-lethal)
if solution.f/wt_growth >= 0.1
    % Then add to True Positive
    tp = tp+1;
else
    % Then add to False Negative
    fn = fn+1;
end

% Simulate 5H2ase KOs of frcA,fruA,hmd,vhuAU,vhcA
% (MMP0820, MMP1382,MMP0127,MMP1694,MMP1693,MMP0823)
ko_model = deleteModelGenes(model,...
    {'MMP0680','MMP0820','MMP1382','MMP0127','MMP1694','MMP1693','MMP0823'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-5H2ase Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (lethal)
if solution.f/wt_growth >= 0.1
    % Then add to False Positive
    fp = fp+1;
else
    % Then add to False Negative
    tn = tn+1;
end

% Simulate 6H2ase KOs of frcA,fruA,hmd,vhuAU,vhcA,ehbN
% (MMP0820, MMP1382, MMP0127,MMP1694,MMP1693,MMP0823,MMP1153)
ko_model = deleteModelGenes(model,...
    {'MMP0680','MMP0820','MMP1382','MMP0127','MMP1694','MMP1693','MMP0823','MMP1153'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-6H2ase Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (lethal)
if solution.f/wt_growth >= 0.1
    % Then add to False Positive
    fp = fp+1;
else
    % Then add to False Negative
    tn = tn+1;
end

% Simulate 6H2ase-cdh KOs of frcAGB,fruAGB,hmd,vhuAU,vhcA,ehbN, and cdh WITH CO supp
% (MMP0820, MMP0818, MMP817, MMP1382, MMP1384, MMP1385, MMP0127,MMP1694,MMP1693,MMP0823,MMP1153,MMP0983-0995)
ko_model = deleteModelGenes(ko_model,...
    {'MMP0983','MMP0984','MMP0985'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-6H2ase-cdh Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Turn on GAPOR
model = changeRxnBounds(model,'rxn07191[c0]',-1000,'l');
model = changeRxnBounds(model,'rxn07191[c0]',1000,'u');

% Simulate 6H2ase_supp KOs of frcAGB,fruAGB,hmd,vhuAU,vhcA,ehbN
% (MMP0820, MMP1382, MMP0127,MMP1694,MMP1693,MMP0823,MMP1153)
ko_model = deleteModelGenes(model,...
    {'MMP0680','MMP0820','MMP1382','MMP0127','MMP1694','MMP1693','MMP0823','MMP1153'}...
    ,0);

solution = optimizeCbModel(ko_model,[],'one');
fprintf('-6H2ase_supp Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate 7H2ase_supp KOs of frcAGB,fruAGB,hmd,vhuAU,vhcA,ehbN,ehaNO
% (MMP0820, MMP1382, MMP0127,MMP1694,MMP1693,MMP0823,MMP1153,MMP1461,MMP1462)
ko_model = deleteModelGenes(model,...
    {'MMP0680','MMP0820','MMP1382','MMP0127','MMP1694','MMP1693','MMP0823','MMP1153','MMP1461','MMP1462'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-7H2ase Growth Ratio: %0.2f\n',solution.f/wt_growth);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Formate simulations
fprintf('\n================================\nGrowth on Formate\n================================');
model = switchToFormateOnly(model);

% Set GAPOR off to begin with
model = changeRxnBounds(model,'rxn07191[c0]',0,'b');

% First simulate Wild-type growth
solution = optimizeCbModel(model,[],'one');
fprintf('\nWild-Type Growth: %0.2f\n\n',solution.f);
wt_growth = solution.f;

% Simulate Hmd KO (MMP0127)
ko_model = deleteModelGenes(model,'MMP0127',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-Hmd Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (non-lethal)
if solution.f/wt_growth >= 0.1
    % Then add to True Positive
    tp = tp+1;
else
    % Then add to False Negative
    fn = fn+1;
end

% Simulate Mtd KO (MMP0372)
ko_model = deleteModelGenes(model,'MMP0372',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-Mtd Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (non-lethal)
if solution.f/wt_growth >= 0.1
    % Then add to True Positive
    tp = tp+1;
else
    % Then add to False Negative
    fn = fn+1;
end

% Simulate FrcA KO (MMP0820)
ko_model = deleteModelGenes(model,'MMP0820',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FrcA Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (non-lethal)
if solution.f/wt_growth >= 0.1
    % Then add to True Positive
    tp = tp+1;
else
    % Then add to False Negative
    fn = fn+1;
end

% Simulate FruA KO (MMP1382)
ko_model = deleteModelGenes(model,'MMP1382',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FruA Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (non-lethal)
if solution.f/wt_growth >= 0.1
    % Then add to True Positive
    tp = tp+1;
else
    % Then add to False Negative
    fn = fn+1;
end

% Simulate FrcA-FruA double KO (MMP0820 and MMP1382)
ko_model = deleteModelGenes(model,{'MMP0820','MMP1382'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FrcA-FruA Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (non-lethal)
if solution.f/wt_growth >= 0.1
    % Then add to True Positive
    tp = tp+1;
else
    % Then add to False Negative
    fn = fn+1;
end

% Simulate  VhuAU-VhcA triple KO (MMP1694, MMP1693, MMP0823)
ko_model = deleteModelGenes(model,{'MMP0680','MMP1694','MMP1693','MMP0823'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-VhuAU-VhcA Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (non-lethal)
if solution.f/wt_growth >= 0.1
    % Then add to True Positive
    tp = tp+1;
else
    % Then add to False Negative
    fn = fn+1;
end

% Simulate  HdrB2 KO (MMP1053)
ko_model = deleteModelGenes(model,{'MMP0680','MMP1053'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-HdrB2 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (non-lethal)
if solution.f/wt_growth >= 0.1
    % Then add to True Positive
    tp = tp+1;
else
    % Then add to False Negative
    fn = fn+1;
end

% Simulate  FdhA1 KO (MMP1298)
ko_model = deleteModelGenes(model,{'MMP1298'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA1 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (non-lethal)
if solution.f/wt_growth >= 0.1
    % Then add to True Positive
    tp = tp+1;
else
    % Then add to False Negative
    fn = fn+1;
end

% Simulate  FdhA2 KO (MMP0138)
ko_model = deleteModelGenes(model,{'MMP0138'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA2 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (non-lethal)
if solution.f/wt_growth >= 0.1
    % Then add to True Positive
    tp = tp+1;
else
    % Then add to False Negative
    fn = fn+1;
end

% Simulate  FdhA1-FdhA2 double KO (MMP1298 and MMP0138)
ko_model = deleteModelGenes(model,{'MMP1298','MMP0138'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA1-FdhA2 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (lethal)
if solution.f/wt_growth >= 0.1
    % Then add to False Positive
    fp = fp+1;
else
    % Then add to False Negative
    tn = tn+1;
end

% Simulate  FdhA2B2 KO (MMP0138 and MMP0139)
ko_model = deleteModelGenes(model,{'MMP0680','MMP0138','MMP0139'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA2B2 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (non-lethal)
if solution.f/wt_growth >= 0.1
    % Then add to True Positive
    tp = tp+1;
else
    % Then add to False Negative
    fn = fn+1;
end

% Simulate EhbF KO (MMP1628)
ko_model = deleteModelGenes(model,'MMP1628',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-EhbF Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate 3H2ase KOs of frcAGB,fruAGB,hmd
% (MMP0820, MMP0818, MMP817, MMP1382, MMP1384, MMP1385, and MMP0127)
ko_model = deleteModelGenes(model,...
    {'MMP0680','MMP0820','MMP0818','MMP0817','MMP1382','MMP1384','MMP1385','MMP0127'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-3H2ase Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (non-lethal)
if solution.f/wt_growth >= 0.1
    % Then add to True Positive
    tp = tp+1;
else
    % Then add to False Negative
    fn = fn+1;
end

% Simulate 5H2ase KOs of frcA,fruA,hmd,vhuAU,vhcA
% (MMP0820, MMP1382,MMP0127,MMP1694,MMP1693,MMP0823)
ko_model = deleteModelGenes(model,...
    {'MMP0680','MMP0820','MMP1382','MMP0127','MMP1694','MMP1693','MMP0823'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-5H2ase Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (lethal)
if solution.f/wt_growth >= 0.1
    % Then add to False Positive
    fp = fp+1;
else
    % Then add to True Negative
    tn = tn+1;
end

% Simulate 6H2ase KOs of frcAGB,fruAGB,hmd,vhuAU,vhcA,ehbN
% (MMP0820, MMP1382, MMP0127,MMP1694,MMP1693,MMP0823,MMP1153)
ko_model = deleteModelGenes(model,...
    {'MMP0680','MMP0820','MMP1382','MMP0127','MMP1694','MMP1693','MMP0823','MMP1153'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-6H2ase Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (lethal)
if solution.f/wt_growth >= 0.1
    % Then add to False Positive
    fp = fp+1;
else
    % Then add to True Negative
    tn = tn+1;
end

% Simulate 6H2ase-cdh KOs of frcAGB,fruAGB,hmd,vhuAU,vhcA,ehbN, and cdh WITH CO supp
% (MMP0820, MMP0818, MMP817, MMP1382, MMP1384, MMP1385, MMP0127,MMP1694,MMP1693,MMP0823,MMP1153,MMP0983-0995)
ko_model = deleteModelGenes(ko_model,...
    {'MMP0983','MMP0984','MMP0985'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-6H2ase-cdh Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Turn on GAPOR
model = changeRxnBounds(model,'rxn07191[c0]',-1000,'l');
model = changeRxnBounds(model,'rxn07191[c0]',1000,'u');

% Simulate 6H2ase_supp KOs of frcAGB,fruAGB,hmd,vhuAU,vhcA,ehbN
% (MMP0820, MMP1382, MMP0127,MMP1694,MMP1693,MMP0823,MMP1153)
ko_model = deleteModelGenes(model,...
    {'MMP0680','MMP0820','MMP1382','MMP0127','MMP1694','MMP1693','MMP0823','MMP1153'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-6H2ase_supp Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (non-lethal)
if solution.f/wt_growth >= 0.1
    % Then add to True Positive
    tp = tp+1;
else
    % Then add to False Negative
    fn = fn+1;
end

% Simulate 7H2ase_supp KOs of frcAGB,fruAGB,hmd,vhuAU,vhcA,ehbN,ehaNO
% (MMP0820, MMP1382, MMP0127,MMP1694,MMP1693,MMP0823,MMP1153,MMP1461,MMP1462)
ko_model = deleteModelGenes(model,...
    {'MMP0680','MMP0820','MMP1382','MMP0127','MMP1694','MMP1693','MMP0823','MMP1153','MMP1461','MMP1462'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-7H2ase Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (non-lethal)
if solution.f/wt_growth >= 0.1
    % Then add to True Positive
    tp = tp+1;
else
    % Then add to False Negative
    fn = fn+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Formate plus H2 simulations
fprintf('\n================================\nGrowth on Formate + H2\n================================');
model = changeRxnBounds(model,'EX_cpd11640[e0]',-1000,'l');

% Set GAPOR off to begin with
model = changeRxnBounds(model,'rxn07191[c0]',0,'b');

% First simulate Wild-type growth
solution = optimizeCbModel(model,[],'one');
fprintf('\nWild-Type Growth: %0.2f\n\n',solution.f);
wt_growth = solution.f;

% Simulate Hmd KO (MMP0127)
ko_model = deleteModelGenes(model,'MMP0127',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-Hmd Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate Mtd KO (MMP0372)
ko_model = deleteModelGenes(model,'MMP0372',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-Mtd Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate FrcA KO (MMP0820)
ko_model = deleteModelGenes(model,'MMP0820',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FrcA Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate FruA KO (MMP1382)
ko_model = deleteModelGenes(model,'MMP1382',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FruA Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate FrcA-FruA double KO (MMP0820 and MMP1382)
ko_model = deleteModelGenes(model,{'MMP0820','MMP1382'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FrcA-FruA Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  VhuAU-VhcA triple KO (MMP1694, MMP1693, MMP0823)
ko_model = deleteModelGenes(model,{'MMP0680','MMP1694','MMP1693','MMP0823'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-VhuAU-VhcA Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  HdrB2 KO (MMP1053)
ko_model = deleteModelGenes(model,{'MMP0680','MMP1053'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-HdrB2 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  FdhA1 KO (MMP1298)
ko_model = deleteModelGenes(model,{'MMP1298'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA1 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  FdhA2 KO (MMP0138)
ko_model = deleteModelGenes(model,{'MMP0138'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA2 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  FdhA1-FdhA2 double KO (MMP1298 and MMP0138)
ko_model = deleteModelGenes(model,{'MMP1298','MMP0138'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA1-FdhA2 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  FdhA2B2 KO (MMP0138 and MMP0139)
ko_model = deleteModelGenes(model,{'MMP0680','MMP0138','MMP0139'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA2B2 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate EhbF KO (MMP1628)
ko_model = deleteModelGenes(model,'MMP1628',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-EhbF Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate 3H2ase KOs of frcAGB,fruAGB,hmd
% (MMP0820, MMP0818, MMP817, MMP1382, MMP1384, MMP1385, and MMP0127)
ko_model = deleteModelGenes(model,...
    {'MMP0680','MMP0820','MMP0818','MMP0817','MMP1382','MMP1384','MMP1385','MMP0127'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-3H2ase Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate 5H2ase KOs of frcA,fruA,hmd,vhuAU,vhcA
% (MMP0820, MMP1382,MMP0127,MMP1694,MMP1693,MMP0823)
ko_model = deleteModelGenes(model,...
    {'MMP0680','MMP0820','MMP1382','MMP0127','MMP1694','MMP1693','MMP0823'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-5H2ase Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (non-lethal)
if solution.f/wt_growth >= 0.1
    % Then add to True Positive
    tp = tp+1;
else
    % Then add to False Negative
    fn = fn+1;
end

% Simulate 6H2ase KOs of frcAGB,fruAGB,hmd,vhuAU,vhcA,ehbN
% (MMP0820, MMP1382, MMP0127,MMP1694,MMP1693,MMP0823,MMP1153)
ko_model = deleteModelGenes(model,...
    {'MMP0680','MMP0820','MMP1382','MMP0127','MMP1694','MMP1693','MMP0823','MMP1153'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-6H2ase Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (non-lethal)
if solution.f/wt_growth >= 0.1
    % Then add to True Positive
    tp = tp+1;
else
    % Then add to False Negative
    fn = fn+1;
end

% Simulate 6H2ase-cdh KOs of frcAGB,fruAGB,hmd,vhuAU,vhcA,ehbN, and cdh WITH CO supp
% (MMP0820, MMP0818, MMP817, MMP1382, MMP1384, MMP1385, MMP0127,MMP1694,MMP1693,MMP0823,MMP1153,MMP0983-0995)
ko_model = deleteModelGenes(ko_model,...
    {'MMP0983','MMP0984','MMP0985'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-6H2ase-cdh Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Turn on GAPOR
model = changeRxnBounds(model,'rxn07191[c0]',-1000,'l');
model = changeRxnBounds(model,'rxn07191[c0]',1000,'u');

% Simulate 6H2ase_supp KOs of frcAGB,fruAGB,hmd,vhuAU,vhcA,ehbN
% (MMP0820, MMP1382, MMP0127,MMP1694,MMP1693,MMP0823,MMP1153)
ko_model = deleteModelGenes(model,...
    {'MMP0680','MMP0820','MMP1382','MMP0127','MMP1694','MMP1693','MMP0823','MMP1153'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-6H2ase_supp Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate 7H2ase_supp KOs of frcAGB,fruAGB,hmd,vhuAU,vhcA,ehbN,ehaNO
% (MMP0820, MMP1382, MMP0127,MMP1694,MMP1693,MMP0823,MMP1153,MMP1461,MMP1462)
ko_model = deleteModelGenes(model,...
    {'MMP0680','MMP0820','MMP1382','MMP0127','MMP1694','MMP1693','MMP0823','MMP1153','MMP1461','MMP1462'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-7H2ase Growth Ratio: %0.2f\n',solution.f/wt_growth);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Formate plus CO simulations
fprintf('\n================================\nGrowth on Formate + CO\n================================');
model = changeRxnBounds(model,'EX_cpd11640[e0]',0,'l');
model = changeRxnBounds(model,'EX_cpd00204[e0]',-1000,'l');

% Set GAPOR off to begin with
model = changeRxnBounds(model,'rxn07191[c0]',0,'b');

% First simulate Wild-type growth
solution = optimizeCbModel(model,[],'one');
fprintf('\nWild-Type Growth: %0.2f\n\n',solution.f);
wt_growth = solution.f;

% Simulate Hmd KO (MMP0127)
ko_model = deleteModelGenes(model,'MMP0127',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-Hmd Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate Mtd KO (MMP0372)
ko_model = deleteModelGenes(model,'MMP0372',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-Mtd Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate FrcA KO (MMP0820)
ko_model = deleteModelGenes(model,'MMP0820',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FrcA Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate FruA KO (MMP1382)
ko_model = deleteModelGenes(model,'MMP1382',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FruA Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate FrcA-FruA double KO (MMP0820 and MMP1382)
ko_model = deleteModelGenes(model,{'MMP0820','MMP1382'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FrcA-FruA Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  VhuAU-VhcA triple KO (MMP1694, MMP1693, MMP0823)
ko_model = deleteModelGenes(model,{'MMP0680','MMP1694','MMP1693','MMP0823'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-VhuAU-VhcA Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  HdrB2 KO (MMP1053)
ko_model = deleteModelGenes(model,{'MMP0680','MMP1053'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-HdrB2 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  FdhA1 KO (MMP1298)
ko_model = deleteModelGenes(model,{'MMP1298'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA1 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  FdhA2 KO (MMP0138)
ko_model = deleteModelGenes(model,{'MMP0138'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA2 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  FdhA1-FdhA2 double KO (MMP1298 and MMP0138)
ko_model = deleteModelGenes(model,{'MMP1298','MMP0138'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA1-FdhA2 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  FdhA2B2 KO (MMP0138 and MMP0139)
ko_model = deleteModelGenes(model,{'MMP0680','MMP0138','MMP0139'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA2B2 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate EhbF KO (MMP1628)
ko_model = deleteModelGenes(model,'MMP1628',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-EhbF Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate 3H2ase KOs of frcAGB,fruAGB,hmd
% (MMP0820, MMP0818, MMP817, MMP1382, MMP1384, MMP1385, and MMP0127)
ko_model = deleteModelGenes(model,...
    {'MMP0680','MMP0820','MMP0818','MMP0817','MMP1382','MMP1384','MMP1385','MMP0127'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-3H2ase Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate 5H2ase KOs of frcA,fruA,hmd,vhuAU,vhcA
% (MMP0820, MMP1382,MMP0127,MMP1694,MMP1693,MMP0823)
ko_model = deleteModelGenes(model,...
    {'MMP0680','MMP0820','MMP1382','MMP0127','MMP1694','MMP1693','MMP0823'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-5H2ase Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate 6H2ase KOs of frcAGB,fruAGB,hmd,vhuAU,vhcA,ehbN
% (MMP0820, MMP1382, MMP0127,MMP1694,MMP1693,MMP0823,MMP1153)
ko_model = deleteModelGenes(model,...
    {'MMP0680','MMP0820','MMP1382','MMP0127','MMP1694','MMP1693','MMP0823','MMP1153'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-6H2ase Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (non-lethal)
if solution.f/wt_growth >= 0.1
    % Then add to True Positive
    tp = tp+1;
else
    % Then add to False Negative
    fn = fn+1;
end

% Simulate 6H2ase-cdh KOs of frcAGB,fruAGB,hmd,vhuAU,vhcA,ehbN, and cdh WITH CO supp
% (MMP0820, MMP0818, MMP817, MMP1382, MMP1384, MMP1385, MMP0127,MMP1694,MMP1693,MMP0823,MMP1153,MMP0983-0995)
ko_model = deleteModelGenes(ko_model,...
    {'MMP0983','MMP0984','MMP0985'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-6H2ase-cdh Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Compare to experimental result (lethal)
if solution.f/wt_growth >= 0.1
    % Then add to False Positive
    fp = fp+1;
else
    % Then add to True Negative
    tn = tn+1;
end

% Turn on GAPOR
model = changeRxnBounds(model,'rxn07191[c0]',-1000,'l');
model = changeRxnBounds(model,'rxn07191[c0]',1000,'u');

% Simulate 6H2ase_supp KOs of frcAGB,fruAGB,hmd,vhuAU,vhcA,ehbN
% (MMP0820, MMP1382, MMP0127,MMP1694,MMP1693,MMP0823,MMP1153)
ko_model = deleteModelGenes(model,...
    {'MMP0680','MMP0820','MMP1382','MMP0127','MMP1694','MMP1693','MMP0823','MMP1153'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-6H2ase_supp Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate 7H2ase_supp KOs of frcAGB,fruAGB,hmd,vhuAU,vhcA,ehbN,ehaNO
% (MMP0820, MMP1382, MMP0127,MMP1694,MMP1693,MMP0823,MMP1153,MMP1461,MMP1462)
ko_model = deleteModelGenes(model,...
    {'MMP0680','MMP0820','MMP1382','MMP0127','MMP1694','MMP1693','MMP0823','MMP1153','MMP1461','MMP1462'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-7H2ase Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Addition on 5/26/2015: Add MCC and accuracy calculation
total = tp+tn+fp+fn;
% First calculate total accuracy
fprintf('\nGene Knockout Accuracy: %0.1f%%(%d/%d)\n',100*(tp+tn)/total,(tp+tn),total);

% Next, calculate MCC
mcc = (tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));
fprintf('Matthews Correlation Coefficient: %0.2f\n',mcc)

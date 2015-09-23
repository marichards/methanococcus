function simulateKOPanel(model)

% For the M. maripaludis S2 model, simulate the model for known gene KO
% experiments to get predictions and compare predictions to reality. Do all
% possible KOs for all 4 conditions

% Alteration on 05/26/2015: Add MCC and accuracy calculations
% Create TP/TN/FP/FN metrics to fill up on appropriate things
tp = 0; tn = 0; fp = 0; fn = 0;

% H2-CO2 simulations
fprintf('================================\nGrowth on H2 + CO2\n================================');

% First simulate Wild-type growth
solution = optimizeCbModel(model,[],'one');
fprintf('\nWild-Type Growth: %0.2f\n\n',solution.f);
wt_growth = solution.f;

% Simulate Hmd KO (mmp0127)
ko_model = deleteModelGenes(model,'mmp0127',0);
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

% Simulate Mtd KO (mmp0372)
ko_model = deleteModelGenes(model,'mmp0372',0);
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

% Simulate FrcA KO (mmp0820)
ko_model = deleteModelGenes(model,'mmp0820',0);
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

% Simulate FruA KO (mmp1382)
ko_model = deleteModelGenes(model,'mmp1382',0);
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

% Simulate FrcA-FruA double KO (mmp0820 and mmp1382)
ko_model = deleteModelGenes(model,{'mmp0820','mmp1382'},0);
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

% Simulate  VhuAU-VhcA triple KO (mmp1694, mmp1693, mmp0823)
ko_model = deleteModelGenes(model,{'mmp0680','mmp1694','mmp1693','mmp0823'},0);
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

% Simulate  HdrB2 KO (mmp1053)
ko_model = deleteModelGenes(model,{'mmp0680','mmp1053'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-HdrB2 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  FdhA1 KO (mmp1298)
ko_model = deleteModelGenes(model,{'mmp1298'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA1 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  FdhA2 KO (mmp0138)
ko_model = deleteModelGenes(model,{'mmp0138'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA2 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  FdhA1-FdhA2 double KO (mmp1298 and mmp0138)
ko_model = deleteModelGenes(model,{'mmp1298','mmp0138'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA1-FdhA2 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  FdhA2B2 KO (mmp0138 and mmp0139)
ko_model = deleteModelGenes(model,{'mmp0680','mmp0138','mmp0139'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA2B2 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate EhbF KO (mmp1628)
ko_model = deleteModelGenes(model,'mmp1628',0);
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
% (mmp0820, mmp0818, mmp817, mmp1382, mmp1384, mmp1385, and mmp0127)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp0818','mmp0817','mmp1382','mmp1384','mmp1385','mmp0127'}...
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
% (mmp0820, mmp1382,mmp0127,mmp1694,mmp1693,mmp0823)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp1382','mmp0127','mmp1694','mmp1693','mmp0823'}...
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
% (mmp0820, mmp1382, mmp0127,mmp1694,mmp1693,mmp0823,mmp1153)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp1382','mmp0127','mmp1694','mmp1693','mmp0823','mmp1153'}...
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
% (mmp0820, mmp0818, mmp817, mmp1382, mmp1384, mmp1385, mmp0127,mmp1694,mmp1693,mmp0823,mmp1153,mmp0983-0995)
ko_model = deleteModelGenes(ko_model,...
    {'mmp0983','mmp0984','mmp0985'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-6H2ase-cdh Growth Ratio: %0.2f\n',solution.f/wt_growth);


% Simulate 6H2ase_supp KOs of frcAGB,fruAGB,hmd,vhuAU,vhcA,ehbN
% (mmp0820, mmp1382, mmp0127,mmp1694,mmp1693,mmp0823,mmp1153)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp1382','mmp0127','mmp1694','mmp1693','mmp0823','mmp1153'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-6H2ase_supp Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate 7H2ase_supp KOs of frcAGB,fruAGB,hmd,vhuAU,vhcA,ehbN,ehaNO
% (mmp0820, mmp1382, mmp0127,mmp1694,mmp1693,mmp0823,mmp1153,mmp1461,mmp1462)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp1382','mmp0127','mmp1694','mmp1693','mmp0823','mmp1153','mmp1461','mmp1462'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-7H2ase Growth Ratio: %0.2f\n',solution.f/wt_growth);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Formate simulations
fprintf('\n================================\nGrowth on Formate\n================================');
model = switchToFormate(model);

% First simulate Wild-type growth
solution = optimizeCbModel(model,[],'one');
fprintf('\nWild-Type Growth: %0.2f\n\n',solution.f);
wt_growth = solution.f;

% Simulate Hmd KO (mmp0127)
ko_model = deleteModelGenes(model,'mmp0127',0);
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

% Simulate Mtd KO (mmp0372)
ko_model = deleteModelGenes(model,'mmp0372',0);
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

% Simulate FrcA KO (mmp0820)
ko_model = deleteModelGenes(model,'mmp0820',0);
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

% Simulate FruA KO (mmp1382)
ko_model = deleteModelGenes(model,'mmp1382',0);
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

% Simulate FrcA-FruA double KO (mmp0820 and mmp1382)
ko_model = deleteModelGenes(model,{'mmp0820','mmp1382'},0);
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

% Simulate  VhuAU-VhcA triple KO (mmp1694, mmp1693, mmp0823)
ko_model = deleteModelGenes(model,{'mmp0680','mmp1694','mmp1693','mmp0823'},0);
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

% Simulate  HdrB2 KO (mmp1053)
ko_model = deleteModelGenes(model,{'mmp0680','mmp1053'},0);
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

% Simulate  FdhA1 KO (mmp1298)
ko_model = deleteModelGenes(model,{'mmp1298'},0);
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

% Simulate  FdhA2 KO (mmp0138)
ko_model = deleteModelGenes(model,{'mmp0138'},0);
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

% Simulate  FdhA1-FdhA2 double KO (mmp1298 and mmp0138)
ko_model = deleteModelGenes(model,{'mmp1298','mmp0138'},0);
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

% Simulate  FdhA2B2 KO (mmp0138 and mmp0139)
ko_model = deleteModelGenes(model,{'mmp0680','mmp0138','mmp0139'},0);
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

% Simulate EhbF KO (mmp1628)
ko_model = deleteModelGenes(model,'mmp1628',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-EhbF Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate 3H2ase KOs of frcAGB,fruAGB,hmd
% (mmp0820, mmp0818, mmp817, mmp1382, mmp1384, mmp1385, and mmp0127)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp0818','mmp0817','mmp1382','mmp1384','mmp1385','mmp0127'}...
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
% (mmp0820, mmp1382,mmp0127,mmp1694,mmp1693,mmp0823)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp1382','mmp0127','mmp1694','mmp1693','mmp0823'}...
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
% (mmp0820, mmp1382, mmp0127,mmp1694,mmp1693,mmp0823,mmp1153)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp1382','mmp0127','mmp1694','mmp1693','mmp0823','mmp1153'}...
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
% (mmp0820, mmp0818, mmp817, mmp1382, mmp1384, mmp1385, mmp0127,mmp1694,mmp1693,mmp0823,mmp1153,mmp0983-0995)
ko_model = deleteModelGenes(ko_model,...
    {'mmp0983','mmp0984','mmp0985'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-6H2ase-cdh Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate 6H2ase_supp KOs of frcAGB,fruAGB,hmd,vhuAU,vhcA,ehbN
% (mmp0820, mmp1382, mmp0127,mmp1694,mmp1693,mmp0823,mmp1153)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp1382','mmp0127','mmp1694','mmp1693','mmp0823','mmp1153'}...
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
% (mmp0820, mmp1382, mmp0127,mmp1694,mmp1693,mmp0823,mmp1153,mmp1461,mmp1462)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp1382','mmp0127','mmp1694','mmp1693','mmp0823','mmp1153','mmp1461','mmp1462'}...
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
model = changeRxnBounds(model,'EX_cpd11640[e0]',-45,'l');

% First simulate Wild-type growth
solution = optimizeCbModel(model,[],'one');
fprintf('\nWild-Type Growth: %0.2f\n\n',solution.f);
wt_growth = solution.f;

% Simulate Hmd KO (mmp0127)
ko_model = deleteModelGenes(model,'mmp0127',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-Hmd Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate Mtd KO (mmp0372)
ko_model = deleteModelGenes(model,'mmp0372',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-Mtd Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate FrcA KO (mmp0820)
ko_model = deleteModelGenes(model,'mmp0820',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FrcA Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate FruA KO (mmp1382)
ko_model = deleteModelGenes(model,'mmp1382',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FruA Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate FrcA-FruA double KO (mmp0820 and mmp1382)
ko_model = deleteModelGenes(model,{'mmp0820','mmp1382'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FrcA-FruA Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  VhuAU-VhcA triple KO (mmp1694, mmp1693, mmp0823)
ko_model = deleteModelGenes(model,{'mmp0680','mmp1694','mmp1693','mmp0823'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-VhuAU-VhcA Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  HdrB2 KO (mmp1053)
ko_model = deleteModelGenes(model,{'mmp0680','mmp1053'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-HdrB2 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  FdhA1 KO (mmp1298)
ko_model = deleteModelGenes(model,{'mmp1298'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA1 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  FdhA2 KO (mmp0138)
ko_model = deleteModelGenes(model,{'mmp0138'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA2 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  FdhA1-FdhA2 double KO (mmp1298 and mmp0138)
ko_model = deleteModelGenes(model,{'mmp1298','mmp0138'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA1-FdhA2 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  FdhA2B2 KO (mmp0138 and mmp0139)
ko_model = deleteModelGenes(model,{'mmp0680','mmp0138','mmp0139'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA2B2 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate EhbF KO (mmp1628)
ko_model = deleteModelGenes(model,'mmp1628',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-EhbF Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate 3H2ase KOs of frcAGB,fruAGB,hmd
% (mmp0820, mmp0818, mmp817, mmp1382, mmp1384, mmp1385, and mmp0127)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp0818','mmp0817','mmp1382','mmp1384','mmp1385','mmp0127'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-3H2ase Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate 5H2ase KOs of frcA,fruA,hmd,vhuAU,vhcA
% (mmp0820, mmp1382,mmp0127,mmp1694,mmp1693,mmp0823)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp1382','mmp0127','mmp1694','mmp1693','mmp0823'}...
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
% (mmp0820, mmp1382, mmp0127,mmp1694,mmp1693,mmp0823,mmp1153)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp1382','mmp0127','mmp1694','mmp1693','mmp0823','mmp1153'}...
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
% (mmp0820, mmp0818, mmp817, mmp1382, mmp1384, mmp1385, mmp0127,mmp1694,mmp1693,mmp0823,mmp1153,mmp0983-0995)
ko_model = deleteModelGenes(ko_model,...
    {'mmp0983','mmp0984','mmp0985'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-6H2ase-cdh Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate 6H2ase_supp KOs of frcAGB,fruAGB,hmd,vhuAU,vhcA,ehbN
% (mmp0820, mmp1382, mmp0127,mmp1694,mmp1693,mmp0823,mmp1153)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp1382','mmp0127','mmp1694','mmp1693','mmp0823','mmp1153'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-6H2ase_supp Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate 7H2ase_supp KOs of frcAGB,fruAGB,hmd,vhuAU,vhcA,ehbN,ehaNO
% (mmp0820, mmp1382, mmp0127,mmp1694,mmp1693,mmp0823,mmp1153,mmp1461,mmp1462)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp1382','mmp0127','mmp1694','mmp1693','mmp0823','mmp1153','mmp1461','mmp1462'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-7H2ase Growth Ratio: %0.2f\n',solution.f/wt_growth);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Formate plus CO simulations
fprintf('\n================================\nGrowth on Formate + CO\n================================');
model = changeRxnBounds(model,'EX_cpd11640[e0]',0,'l');
model = changeRxnBounds(model,'EX_cpd00204[e0]',-45,'l');

% First simulate Wild-type growth
solution = optimizeCbModel(model,[],'one');
fprintf('\nWild-Type Growth: %0.2f\n\n',solution.f);
wt_growth = solution.f;

% Simulate Hmd KO (mmp0127)
ko_model = deleteModelGenes(model,'mmp0127',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-Hmd Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate Mtd KO (mmp0372)
ko_model = deleteModelGenes(model,'mmp0372',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-Mtd Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate FrcA KO (mmp0820)
ko_model = deleteModelGenes(model,'mmp0820',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FrcA Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate FruA KO (mmp1382)
ko_model = deleteModelGenes(model,'mmp1382',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FruA Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate FrcA-FruA double KO (mmp0820 and mmp1382)
ko_model = deleteModelGenes(model,{'mmp0820','mmp1382'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FrcA-FruA Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  VhuAU-VhcA triple KO (mmp1694, mmp1693, mmp0823)
ko_model = deleteModelGenes(model,{'mmp0680','mmp1694','mmp1693','mmp0823'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-VhuAU-VhcA Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  HdrB2 KO (mmp1053)
ko_model = deleteModelGenes(model,{'mmp0680','mmp1053'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-HdrB2 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  FdhA1 KO (mmp1298)
ko_model = deleteModelGenes(model,{'mmp1298'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA1 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  FdhA2 KO (mmp0138)
ko_model = deleteModelGenes(model,{'mmp0138'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA2 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  FdhA1-FdhA2 double KO (mmp1298 and mmp0138)
ko_model = deleteModelGenes(model,{'mmp1298','mmp0138'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA1-FdhA2 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate  FdhA2B2 KO (mmp0138 and mmp0139)
ko_model = deleteModelGenes(model,{'mmp0680','mmp0138','mmp0139'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA2B2 Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate EhbF KO (mmp1628)
ko_model = deleteModelGenes(model,'mmp1628',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-EhbF Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate 3H2ase KOs of frcAGB,fruAGB,hmd
% (mmp0820, mmp0818, mmp817, mmp1382, mmp1384, mmp1385, and mmp0127)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp0818','mmp0817','mmp1382','mmp1384','mmp1385','mmp0127'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-3H2ase Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate 5H2ase KOs of frcA,fruA,hmd,vhuAU,vhcA
% (mmp0820, mmp1382,mmp0127,mmp1694,mmp1693,mmp0823)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp1382','mmp0127','mmp1694','mmp1693','mmp0823'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-5H2ase Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate 6H2ase KOs of frcAGB,fruAGB,hmd,vhuAU,vhcA,ehbN
% (mmp0820, mmp1382, mmp0127,mmp1694,mmp1693,mmp0823,mmp1153)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp1382','mmp0127','mmp1694','mmp1693','mmp0823','mmp1153'}...
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
% (mmp0820, mmp0818, mmp817, mmp1382, mmp1384, mmp1385, mmp0127,mmp1694,mmp1693,mmp0823,mmp1153,mmp0983-0995)
ko_model = deleteModelGenes(ko_model,...
    {'mmp0983','mmp0984','mmp0985'},0);
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

% Simulate 6H2ase_supp KOs of frcAGB,fruAGB,hmd,vhuAU,vhcA,ehbN
% (mmp0820, mmp1382, mmp0127,mmp1694,mmp1693,mmp0823,mmp1153)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp1382','mmp0127','mmp1694','mmp1693','mmp0823','mmp1153'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-6H2ase_supp Growth Ratio: %0.2f\n',solution.f/wt_growth);

% Simulate 7H2ase_supp KOs of frcAGB,fruAGB,hmd,vhuAU,vhcA,ehbN,ehaNO
% (mmp0820, mmp1382, mmp0127,mmp1694,mmp1693,mmp0823,mmp1153,mmp1461,mmp1462)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp1382','mmp0127','mmp1694','mmp1693','mmp0823','mmp1153','mmp1461','mmp1462'}...
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

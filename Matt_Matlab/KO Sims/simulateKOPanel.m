function simulateKOPanel(model)

% For the M. maripaludis S2 model, simulate the model for known gene KO
% experiments to get predictions and compare predictions to reality

% H2-CO2 simulations
fprintf('================================\nGrowth on H2 + CO2\n================================');

% First simulate Wild-type growth
solution = optimizeCbModel(model,[],'one');
fprintf('\nWild-Type Growth: %f\n\n',solution.f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hendrickson et al. (2008) J. Bacteriology:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulate Hmd KO (mmp0127)
ko_model = deleteModelGenes(model,'mmp0127',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-Hmd Growth: %f\n',solution.f);

% Simulate Mtd KO (mmp0372)
ko_model = deleteModelGenes(model,'mmp0372',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-Mtd Growth: %f\n',solution.f);

% Simulate FrcA KO (mmp0820)
ko_model = deleteModelGenes(model,'mmp0820',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FrcA Growth: %f\n',solution.f);

% Simulate FruA KO (mmp1382)
ko_model = deleteModelGenes(model,'mmp1382',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FruA Growth: %f\n',solution.f);

% Simulate FrcA-FruA double KO (mmp0820 and mmp1382)
ko_model = deleteModelGenes(model,{'mmp0820','mmp1382'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FrcA-FruA Growth: %f\n',solution.f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Costa et al. (2010) PNAS: MM901 (-mmp0680)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulate  VhuAU-VhcA triple KO (mmp1694, mmp1693, mmp0823)
ko_model = deleteModelGenes(model,{'mmp0680','mmp1694','mmp1693','mmp0823'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-VhuAU-VhcA Growth: %f\n',solution.f);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Porat et al. (2006) J. Bacteriology:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulate EhbF KO (mmp1628)
ko_model = deleteModelGenes(model,'mmp1628',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-EhbF Growth: %f\n',solution.f);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Costa et al. (2013) Mbio: MM901 (-mmp0680)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulate 3H2ase KOs of frcAGB,fruAGB,hmd
% (mmp0820, mmp0818, mmp817, mmp1382, mmp1384, mmp1385, and mmp0127)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp0818','mmp0817','mmp1382','mmp1384','mmp1385','mmp0127'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-3H2ase Growth: %f\n',solution.f);

% Simulate 5H2ase KOs of frcA,fruA,hmd,vhuAU,vhcA
% (mmp0820, mmp1382,mmp0127,mmp1694,mmp1693,mmp0823)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp1382','mmp0127','mmp1694','mmp1693','mmp0823'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-5H2ase Growth: %f\n',solution.f);

% Simulate 6H2ase KOs of frcAGB,fruAGB,hmd,vhuAU,vhcA,ehbN
% (mmp0820, mmp1382, mmp0127,mmp1694,mmp1693,mmp0823,mmp1153)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp1382','mmp0127','mmp1694','mmp1693','mmp0823','mmp1153'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-6H2ase Growth: %f\n',solution.f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Formate simulations
fprintf('\n================================\nGrowth on Formate\n================================');
model = switchToFormate(model);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Costa et al. (2010) PNAS: MM901 (-mmp0680)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First simulate Wild-type growth
solution = optimizeCbModel(model,[],'one');
fprintf('\nWild-Type Growth: %f\n\n',solution.f);

% Simulate  VhuAU-VhcA triple KO (mmp1694, mmp1693, mmp0823)
ko_model = deleteModelGenes(model,{'mmp0680','mmp1694','mmp1693','mmp0823'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-VhuAU-VhcA Growth: %f\n',solution.f);

% Simulate  HdrB2 KO (mmp1053)
ko_model = deleteModelGenes(model,{'mmp0680','mmp1053'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-HdrB2 Growth: %f\n',solution.f);

% Simulate  FdhA2B2 KO (mmp0138 and mmp0139)
ko_model = deleteModelGenes(model,{'mmp0680','mmp0138','mmp0139'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA2B2 Growth: %f\n',solution.f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lupa et al. (2008) App. and Env. Microbiol.:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulate  FdhA1 KO (mmp1298)
ko_model = deleteModelGenes(model,{'mmp1298'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA1 Growth: %f\n',solution.f);

% Simulate  FdhA2 KO (mmp0138)
ko_model = deleteModelGenes(model,{'mmp0138'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA2 Growth: %f\n',solution.f);

% Simulate  FdhA1-FdhA2 double KO (mmp1298 and mmp0138)
ko_model = deleteModelGenes(model,{'mmp1298','mmp0138'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FdhA1-FdhA2 Growth: %f\n',solution.f);

% Simulate Hmd KO (mmp0127)
ko_model = deleteModelGenes(model,'mmp0127',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-Hmd Growth: %f\n',solution.f);

% Simulate Mtd KO (mmp0372)
ko_model = deleteModelGenes(model,'mmp0372',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-Mtd Growth: %f\n',solution.f);

% Simulate FrcA KO (mmp0820)
ko_model = deleteModelGenes(model,'mmp0820',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FrcA Growth: %f\n',solution.f);

% Simulate FruA KO (mmp1382)
ko_model = deleteModelGenes(model,'mmp1382',0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FruA Growth: %f\n',solution.f);

% Simulate FrcA-FruA double KO (mmp0820 and mmp1382)
ko_model = deleteModelGenes(model,{'mmp0820','mmp1382'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-FrcA-FruA Growth: %f\n',solution.f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lie et al. (2012) PNAS: MM901 (-mmp0680)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulate 3H2ase KOs of frcAGB,fruAGB,hmd
% (mmp0820, mmp0818, mmp817, mmp1382, mmp1384, mmp1385, and mmp0127)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp0818','mmp0817','mmp1382','mmp1384','mmp1385','mmp0127'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-3H2ase Growth: %f\n',solution.f);

% Simulate 5H2ase KOs of frcA,fruA,hmd,vhuAU,vhcA
% (mmp0820, mmp1382, mmp0127,mmp1694,mmp1693,mmp0823)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp1382','mmp0127','mmp1694','mmp1693','mmp0823'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-5H2ase Growth: %f\n',solution.f);

% Simulate 6H2ase KOs of frcA,fruA,hmd,vhuAU,vhcA,ehbN
% (mmp0820,mmp1382,mmp0127,mmp1694,mmp1693,mmp0823,mmp1153)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp1382','mmp0127','mmp1694','mmp1693','mmp0823','mmp1153'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-6H2ase Growth: %f\n',solution.f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Costa et al. (2013) Mbio: MM901 (-mmp0680)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulate 6H2ase KOs of frcAGB,fruAGB,hmd,vhuAU,vhcA,ehbN WITH CO supp
% (mmp0820, mmp0818, mmp817, mmp1382, mmp1384, mmp1385, mmp0127,mmp1694,mmp1693,mmp0823,mmp1153)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp1382','mmp0127','mmp1694','mmp1693','mmp0823','mmp1153'}...
    ,0);
ko_model = changeRxnBounds(ko_model,'EX_cpd00204_e0',-10,'l');
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-6H2ase + 10 mmol/gDCW/h CO Growth: %f\n',solution.f);

% Simulate 6H2ase-cdh KOs of frcAGB,fruAGB,hmd,vhuAU,vhcA,ehbN, and cdh WITH CO supp
% (mmp0820, mmp0818, mmp817, mmp1382, mmp1384, mmp1385, mmp0127,mmp1694,mmp1693,mmp0823,mmp1153,mmp0983-0995)
ko_model = deleteModelGenes(ko_model,...
    {'mmp0983','mmp0984','mmp0985'},0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-6H2ase-cdh + 10 mmol/gDCW/h CO Growth: %f\n',solution.f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Formate plus H2 simulations
fprintf('\n================================\nGrowth on Formate + H2\n================================');
model = changeRxnBounds(model,'Ex_cpd11640_c0',-45,'l');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lie et al. (2012) PNAS: MM901 (-mmp0680)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First simulate Wild-type growth
solution = optimizeCbModel(model,[],'one');
fprintf('\nWild-Type Growth: %f\n\n',solution.f);

% Simulate 5H2ase KOs of frcA,fruA,hmd,vhuAU,vhcA
% (mmp0820, mmp1382, mmp0127,mmp1694,mmp1693,mmp0823)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp1382','mmp0127','mmp1694','mmp1693','mmp0823'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-5H2ase Growth: %f\n',solution.f);

% Simulate 6H2ase KOs of frcAGB,fruAGB,hmd,vhuAU,vhcA,ehbN
% (mmp0820, mmp0818, mmp817, mmp1382, mmp1384, mmp1385, mmp0127,mmp1694,mmp1693,mmp0823,mmp1153)
ko_model = deleteModelGenes(model,...
    {'mmp0680','mmp0820','mmp1382','mmp0127','mmp1694','mmp1693','mmp0823','mmp1153'}...
    ,0);
solution = optimizeCbModel(ko_model,[],'one');
fprintf('-6H2ase Growth: %f\n',solution.f);






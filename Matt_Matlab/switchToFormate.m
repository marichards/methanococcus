function model = switchToFormate(model)

%Switches M.maripaludis model from H2 to formate 

%Note: Formate is HCO2
%Turn off H2 Input
model = changeRxnBounds(model,'Ex_cpd11640_c0',0,'l');
%Turn on Formate Input
model = changeRxnBounds(model,'EX_cpd00047_e0',-10,'l');

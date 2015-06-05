function model = uptakeMethaneEatSulfate(model)

%Add in the reaction for converting methanol to methyl-coM
model = addReaction(model,'Methanol_to_MCoM',...
    'Methanol_c0 + CoM_c0 + H_c0 <=> Methyl_CoM_c0 + H2O_c0');

%Add in the methanol uptake; this is now connected to the methanogenesis
%pathway
model = addReaction(model,'EX_methanol_c0',...
    'Methanol_c0 <=> ');

%Turn off the MFR
%model = changeRxnBounds(model,'rxn11938_c0',0,'b');

%Change the media to SUPPLY methane

%Overall methanol reaction(s): 
%H2O + SO4 + 5 CH4 = H2S + 5 CH3OH

% From paper:
%CH4 + SO4(-2) -> HCO3- + HS- + H2O
%CH3 + 4 NO3- -> CO2 + 4 NO2- + 2 H2O

%Add in sulfur supply
model = addReaction(model,{'rxn05651_c0','sulfate transport in via proton symport c0'},...
    'H_e0 + Sulfate_e0 <=> H_c0 + Sulfate_c0');
model = addReaction(model,{'EX_cpd00048_e0','EX_Sulfate_e0'},...
    'Sulfate_e0 <=> ');

%Change bounds such that methane goes IN instead of OUT
model = changeRxnBounds(model,'Ex_cpd01024_c0',-1000,'l');
model = changeRxnBounds(model,'Ex_cpd01024_c0',0,'u');
%Change methanol to come OUT instead of IN
model = changeRxnBounds(model,'EX_methanol_c0',0,'l');
model = changeRxnBounds(model,'EX_methanol_c0',1000,'u');

%Turn off Hydrogen input and let it come out
model = changeRxnBounds(model,'Ex_cpd11640_c0',0,'l');
model = changeRxnBounds(model,'Ex_cpd11640_c0',1000,'u');

end
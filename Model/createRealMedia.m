function model = createRealMedia(model)

%Take in the Methanococcus model and remove unnecessary exchanges
%reactions, add necessary exchanges that are missing, and turn off uptake
%on others
%
%Inputs
%model - the Methanococcus metabolic model with rich H2/CO2 media
%
%Outputs
%model - the new Methanococcus metabolic model with proper H2/CO2 media

%Take away a bunch of reactions as media by setting their lower bounds to 0

%Non-essential (KO)
%EX Urea e0	-1000
model = removeRxns(model,'EX_cpd00073_e0');
%EX fe3 e0	-1000
model = removeRxns(model,'EX_cpd10516_e0');
%EX Spermine e0	-1000
model = removeRxns(model,'EX_cpd00558_e0');
%EX Nitrate e0	-1000
model = removeRxns(model,'EX_cpd00209_e0');
%EX BET e0	-1000
model = removeRxns(model,'EX_cpd00540_e0');
%EX Cytosine e0	-1000
model = removeRxns(model,'EX_cpd00307_e0');
%EX glycogenn-1 c0	-1000
model = removeRxns(model,'EX_cpd15302_c0');
%EX Uracil e0	-1000
model = removeRxns(model,'EX_cpd00092_e0');
%EX ddca e0	-1000
model = removeRxns(model,'EX_cpd01741_e0');
%EX Dephospho-CoA e0	-1000
model = removeRxns(model,'EX_cpd00655_e0');
%EX Cobinamide e0	-1000 %Necessary for Colamide?
model = removeRxns(model,'EX_cpd03422_e0');
%EX Oxidized glutathione e0	-1000
model = removeRxns(model,'EX_cpd00111_e0');

%Turn off uptake, don't take out
%EX H e0	-1000
model = changeRxnBounds(model,'EX_cpd00067_e0',0,'l');

%Metals in the biomass that must be removed
%EX Zn2 e0	-1000
%model = changeRxnBounds(model,'EX_cpd00034_e0',0,'l');
%EX Cu2 e0	-1000
%model = changeRxnBounds(model,'EX_cpd00058_e0',0,'l');
%EX Mn2 e0	-1000
%model = changeRxnBounds(model,'EX_cpd00030_e0',0,'l');



%Non-Essential but maybe should be?
%EX Nicotinamide ribonucleotide e0	-1000
model = changeRxnBounds(model,'EX_cpd00355_e0',0,'l');
%EX octadecenoate e0	-1000
%Taking this out screws up methane production...why?
%model = changeRxnBounds(model,'EX_cpd15269_e0',0,'l');




%Currently Essential
%EX Co2 e0	-1000
%Co2_c0, Peptidoglycan_polymer_n_subunits_c0, Calomide_c0, ACP_c0
%model = changeRxnBounds(model,'EX_cpd00149_e0',0,'l');
%EX Thiamin e0	-1000
%TPP_c0, Peptidoglycan_polymer_n_subunits_c0, Calomide_c0, ACP_c0
%model = changeRxnBounds(model,'EX_cpd00305_e0',0,'l');

%%%%%%%%%%%%%%%
%Add bicarbonate uptake





function model = switchToSpecificFd(model)

%Take in the MM model, switch it so that it uses specific ferredoxins
%instead of general ones
%
%Input
%model: the Methanococcus maripaludis S2 model with general Fds
%
%Output
%model: the Methanococcus maripaludis S2 model with specific Fds

%Commands from before
%Group EhA, HdrABC, and rxn11938_c0 together by giving them a specific
%ferredoxin
[~,idx] = intersect(model.rxns,'Eha/Ehb');
model.rxns{idx} = 'Eha';
model = addReaction(model,{'Eha','Eha'},...
    'Fdox*1_c0 + 2.000000 Na_e0 + H2_c0 <=>	2.000000 H_c0 + Fdred*1_c0 + 2.000000 Na_c0');
model = addReaction(model,{'HdrABC','Heterodisulfide reductase'},...
    'Fdox*1_c0 + CoM-S-S-CoB_c0 + 2.000000 H2_c0 ->	2.000000 H_c0 + Fdred*1_c0 + CoM_c0 + HTP_c0'); 
[~,idx] = intersect(model.rxns,'rxn11938_c0');
name = model.rxnNames{idx};
model = addReaction(model,{'rxn11938_c0',name},...
    'H2O_c0 + Fdox*1_c0 + Formylmethanofuran_c0	<=>	H_c0 + CO2_c0 + Fdred*1_c0 + Methanofuran_c0');

%Add EhB, indolepyruvate oxidoreductase with new specific ferredoxin; both are dead for now	
model = addReaction(model,'Ehb',...
    'Fdox*2_c0 + 2.000000 Na_e0 + H2_c0 <=>	2.000000 H_c0 + Fdred*2_c0 + 2.000000 Na_c0');
model = addReaction(model,{'rxn10561_c0','Indolepyruvate ferredoxin oxidoreductase'},...
    'Indole-3-pyruvate_c0 + Fdox*2_c0 + CoA_c0 <=> S-2-(indol-3-yl)acetyl-CoA_c0 + CO2_c0 + Fdred*2_c0 + H_c0');
%Give the CODH, rxn05938, and rxn05939 the same ferredoxin
[~,idx] = intersect(model.rxns,'CODH');
name = model.rxnNames{idx};	
model = addReaction(model,{'CODH',name},...
    'CO2_c0 + CoA_c0 + 2 H_c0 + Fdred*2_c0 + 5-Methyl-H4MPT_c0 <=> Acetyl-CoA_c0 + Fdox*2_c0 + H2O_c0 + H4MPT_c0');
[~,idx] = intersect(model.rxns,'rxn05938_c0');
name = model.rxnNames{idx};
model = addReaction(model,{'rxn05938_c0',name},...
    'H_c0 + CO2_c0 + Acetyl-CoA_c0 + Fdred*2_c0 <=>	CoA_c0 + Fdox*2_c0 + Pyruvate_c0');
[~,idx] = intersect(model.rxns,'rxn05939_c0');
name = model.rxnNames{idx};
model = addReaction(model,{'rxn05939_c0','2-oxoglutarate synthase'},...
    'H_c0 + CO2_c0 + Fdred*2_c0 + Succinyl-CoA_c0 	<=>	CoA_c0 + 2-Oxoglutarate_c0 + Fdox*2_c0');
%Add the Formate:Hdr
model = addReaction(model,{'Hdr_formate','Formate-utilizing heterodisulfide reductase'},...
    'CoM-S-S-CoB_c0 + 2 Formate_c0 + Fdox*1_c0 -> 2 CO2_c0 + CoM_c0 + HTP_c0 + Fdred*1_c0');

	
%Add charges for last 4 metabolites (4 ferredoxins)
%    'Fdox*1_c0'
model.metCharge(end-3)=6;
%    'Fdred*1_c0'
model.metCharge(end-2)=4;
%    'Fdox*2_c0'
model.metCharge(end-1)=6;
%    'Fdred*2_c0'
model.metCharge(end)=4;


function model = switchToSpecificFd(model)

% Take in the M. maripaludis model, and switches the ferredoxins on 3 Wolfe
% Cycle reactions to specific ferredoxins instead of general ones. Also
% adds reactions that allows specific ferredoxins to convert to general
% ones, but not visa versa.
%
% Input
% model: the M. maripaludis model, a COBRA Toolbox model structure (with
% general ferredoxins)
%
% Output
% model: the M. maripaludis model with specific Fds in 3 reactions and 2
% reactions that allow specific ferredoxins to substitute for general
% ferredoxins
%
% Matthew Richards, 09/29/2015


% Group EhA, HDRs, and rxn11938[c0] together by giving them a specific
% ferredoxin
% First, break off the former Eha/Ehb as just Ehb
[~,idx] = intersect(model.rxns,'Eha/Ehb');
model.rxns{idx} = 'Ehb';
model.rxnNames{idx} = 'Ehb';

% Add Eha
model = addReaction(model,{'Eha','Eha'},...
    'Fdox*1[c0] + 2.000000 cpd00971[e0] + cpd11640[c0] <=> 2.000000 cpd00067[c0] + Fdred*1[c0] + 2.000000 cpd00971[c0]');
[~,idx] = intersect(model.rxns,'Eha');
model.freeEnergy(idx) = 0;

% Add specific ferredoxins for Hdrs and rxn11938 (Fwd)
[~,idx] = intersect(model.rxns,'HdrABC');
name = model.rxnNames{idx};
model = addReaction(model,{'HdrABC',name},...
    'Fdox*1[c0] + cpd02935[c0] + 2.000000 cpd11640[c0] -> 2.000000 cpd00067[c0] + Fdred*1[c0] + cpd02246[c0] + cpd02817[c0]'); 
[~,idx] = intersect(model.rxns,'rxn11938[c0]');
name = model.rxnNames{idx};
model = addReaction(model,{'rxn11938[c0]',name},...
    'cpd00001[c0] + Fdox*1[c0] + cpd00735[c0]	<=>	cpd00067[c0] + cpd00011[c0] + Fdred*1[c0] + cpd00643[c0]');
[~,idx] = intersect(model.rxns,'Hdr_formate');
name = model.rxnNames{idx};
model = addReaction(model,{'Hdr_formate',name},...
    'cpd02935[c0] + 2 cpd00047[c0] + Fdox*1[c0] -> 2 cpd00011[c0] + cpd02246[c0] + cpd02817[c0] + Fdred*1[c0]');

% Add charges for added metabolites (specific ferredoxins)
%    'Fdox*1[c0]'
model.metCharge(end-3)=6;
%    'Fdred*1[c0]'
model.metCharge(end-2)=4;

% Allow specific ferredoxin to substitute for general in a pinch
model = addReaction(model,'Oxidized_Fd_promiscuity',...
    'Fdox*1[c0] -> cpd11621[c0]');
model = addReaction(model,'Reduced_Fd_promiscuity',...
    'Fdred*1[c0] -> cpd11620[c0]');

% Add no free energy to either one
[~,idx] = intersect(model.rxns,'Oxidized_Fd_promiscuity');
model.freeEnergy(idx) = 0;
[~,idx] = intersect(model.rxns,'Reduced_Fd_promiscuity');
model.freeEnergy(idx) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Another option is to add a second specific ferredoxin to another set of
% reactions. This is currently not being enforced
% other reactions together
%%Add EhB, indolepyruvate oxidoreductase with new specific ferredoxin; both are dead for now	
%model = addReaction(model,'Ehb',...
%    'Fdox*2[c0] + 2.000000 Na_e0 + H2[c0] <=>	2.000000 H[c0] + Fdred*2[c0] + 2.000000 Na[c0]');
%[~,idx] = intersect(model.rxns,'rxn10561[c0]');
%name = model.rxnNames{idx};
%model = addReaction(model,{'rxn10561[c0]',name},...
%    'Indole-3-pyruvate[c0] + Fdox*2[c0] + CoA[c0] <=> S-2-(indol-3-yl)acetyl-CoA[c0] + CO2[c0] + Fdred*2[c0] + H[c0]');
%%Give the CODH, rxn05938, and rxn05939 the same ferredoxin
%[~,idx] = intersect(model.rxns,'CODH_ACS');
%name = model.rxnNames{idx};	
%model = addReaction(model,{'CODH_ACS',name},...
%    'CO2[c0] + CoA[c0] + 2 H[c0] + Fdred*2[c0] + 5-Methyl-H4MPT[c0] <=> Acetyl-CoA[c0] + Fdox*2[c0] + H2O[c0] + H4MPT[c0]');
%[~,idx] = intersect(model.rxns,'rxn05938[c0]');
%name = model.rxnNames{idx};
%model = addReaction(model,{'rxn05938[c0]',name},...
%    'H[c0] + CO2[c0] + Acetyl-CoA[c0] + Fdred*2[c0] <=>	CoA[c0] + Fdox*2[c0] + Pyruvate[c0]');
%[~,idx] = intersect(model.rxns,'rxn05939[c0]');
%name = model.rxnNames{idx};
%model = addReaction(model,{'rxn05939[c0]',name},...
%    'H[c0] + CO2[c0] + Fdred*2[c0] + Succinyl-CoA[c0] 	<=>	CoA[c0] + 2-Oxoglutarate[c0] + Fdox*2[c0]');
%%%    'Fdox*2[c0]'
%model.metCharge(end-1)=6;
%%%    'Fdred*2[c0]'
%model.metCharge(end)=4;
	


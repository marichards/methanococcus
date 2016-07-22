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
[~,idx] = intersect(model.mets,'Fdox*1[c0]');
model.metCharge(idx)=6;
%    'Fdred*1[c0]'
[~,idx] = intersect(model.mets,'Fdred*1[c0]');
model.metCharge(idx)=4;

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

% 11/23/2015
% Enforce these constraints as well
%%Add EhB, indolepyruvate oxidoreductase with new specific ferredoxin; both are dead for now	
model = addReaction(model,{'Ehb','Ehb'},...
    'Fdox*2[c0] + 2.000000 cpd00971[e0] + cpd11640[c0] <=> 2.000000 cpd00067[c0] + Fdred*2[c0] + 2.000000 cpd00971[c0]');
[~,idx] = intersect(model.rxns,'Ehb');
model.freeEnergy(idx) = 0;

[~,idx] = intersect(model.rxns,'rxn10561[c0]');
name = model.rxnNames{idx};
model = addReaction(model,{'rxn10561[c0]',name},...
    'cpd00010[c0] + Fdox*2[c0] + cpd00278[c0] 	<=>	cpd00067[c0] + cpd00011[c0] + Fdred*2[c0] + cpd11175[c0]');
%%Give the CODH, rxn05938, and rxn05939 the same ferredoxin
[~,idx] = intersect(model.rxns,'CODH');
name = model.rxnNames{idx};	
model = addReaction(model,{'CODH',name},...
    '3 cpd00067[c0] + cpd00011[c0] + Fdred*2[c0] 	<=>	cpd00001[c0] + Fdox*2[c0] + 0.500000 cpd11640[c0] + cpd00204[c0]');

[~,idx] = intersect(model.rxns,'rxn05938[c0]');
name = model.rxnNames{idx};
model = addReaction(model,{'rxn05938[c0]',name},...
    'cpd00067[c0] + cpd00011[c0] + cpd00022[c0] + Fdred*2[c0] 	<=>	cpd00010[c0] + Fdox*2[c0] + cpd00020[c0]');

[~,idx] = intersect(model.rxns,'rxn05939[c0]');
name = model.rxnNames{idx};
model = addReaction(model,{'rxn05939[c0]',name},...
    'cpd00067[c0] + cpd00011[c0] + Fdred*2[c0] + cpd00078[c0] 	<=>	cpd00010[c0] + cpd00024[c0] + Fdox*2[c0]');

% Add charges for added metabolites (specific ferredoxins)
%    'Fdox*2[c0]'
[~,idx] = intersect(model.mets,'Fdox*2[c0]');
model.metCharge(idx)=6;
%    'Fdred*2[c0]'
[~,idx] = intersect(model.mets,'Fdred*2[c0]');
model.metCharge(idx)=4;
	


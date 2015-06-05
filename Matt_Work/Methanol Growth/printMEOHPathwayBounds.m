function rxns = printMEOHPathwayBounds(model)

% Take in the MeOH-consuming model and print out the core metabolic pathway
% along with the bounds of each reaction


% Pull out the whole pathway
rxns = {...
    % Hydrogen exchange
    'Ex_cpd11640_c0'
    % Methanol exchange
    'Methanol_supply';...
    % MeOH --> Methyl-CoM
    'Methanol_to_MCoM';...
    % Methyl-CoM --> Methane
    'rxn03127_c0';...
    % CoB-S-S-CoM --> CoM + CoB
    'HdrABC';...
    % Methyl-CoM <=> Methyl-H4MPT
    'rxn03020_c0';...
    % Methyl-H4MPT <=> Methylene-H4MPT
    'rxn03085_c0';...
    % Methylene-H4MPT <=> Methenyl-H4MPT (two reactions)
    'rxn06696_c0';'rxn03079_c0';...
    % Methenyl-H4MPT <=> Formyl-H4MPT 
    'rxn02480_c0';...
    % Formyl-H4MPT <=> Formyl-MFR
    'rxn02431_c0';...
    % Formyl-MFR <=> CO2
    'rxn11938_c0';...
    % Methyl-H4MPT <=> Acetyl-CoA
    'ACS';...
    % Eha/Ehb
    'Eha/Ehb';...
    % Fru/Frc
    'rxn06299_c0';...
    % ATP Synthase
    'ATPS';...
    % H+ and Na+ pump
    'rxn05209_c0';...
    };

% Clear space a bit
fprintf('\n');
% Print a header
fprintf('Lower\tUpper\tReaction Formula\n\n')
% Print them out with bounds
for i=1:length(rxns)
    %Find the index
    [~,idx] = intersect(model.rxns,rxns{i});
    % Find the formula
    formula = printRxnFormula(model,rxns{i},false);
    %Print out the bounds and formula for the reaction
    fprintf('%d\t%d\t%s\n',model.lb(idx),model.ub(idx),...
        formula{1});
    
end


    



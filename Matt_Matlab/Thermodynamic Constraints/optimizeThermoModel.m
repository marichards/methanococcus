function [solution,gibbs_flux,model] = optimizeThermoModel(model,substrateRxns,concentrations,T,water_rxn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Version 4: 06/23/2015

%GENERAL METHODOLOGY
%Accepts a model, substrate reactions, and initial concentrations
%Accepts a reaction for water
%Adds a metabolite called "dG[e]"
%Adds an exchange for dG called 'GIBBS' (EX_dG[e])
%Modifies exchanges to include dG[e] by adding a row
%Row values: adjusted for concentration (or biomass amount)
%Simulates the model
%Returns solution and overall dG
%(This is a modified addOverallDG code)

% Inputs
% model: a COBRA model structure
% substrateRxns: set of reactions with specified concentrations
% concentrations: set of concentrations specified for exchange reactions
% in the model simulation (in mM)
% T: temperature in Kelvin
% 

%%%%%%%%%%%%%%%%dG values are from Alberty, at pH=7.0%%%%%%%%%%%%%%%%%%
%Catch concentrations that are 0....
if any(~concentrations)
    solution = optimizeCbModel(model,[],'one');
    gibbs_flux = inf;
    
else
    
    %Error for things not the same size
    if length(substrateRxns) ~= length(concentrations)

        error('substrateRxns and concentrations must be of equal length')
    end
    %Gas constant specification (Maybe move these out?)
    R = 8.314e-6; %kJ/mmol*K

    %Add the new reaction first, which adds the metabolite
    model = addReaction(model,'GIBBS_kJ/GDW','dG <=> '); 

    %Find index of dG
    [~,met_idx] = intersect(model.mets,'dG');

    % Alter the free energy values for things with substrate reactions in
    % the free energy vector itself
    % First grab the index of the exchange reactions in the model
    [rxns,rxn_idx] = intersect(model.rxns,substrateRxns,'stable');
    % Make a dictionary
    dict = containers.Map(rxns,rxn_idx);      
    % For those indices, change the free energy numbers using concentration
    %Loop: put in the correct free energy term for each:
    %dG = dG_0 + RTln(C)
    for i = 1:length(substrateRxns)
       %Change the dG weight for the exchange reaction (Conc in mM)
       model.freeEnergy(dict(substrateRxns{i})) = model.freeEnergy(dict(substrateRxns{i}))...
           +R*T*log(concentrations(i));    
    end
    
    % Add free energy values to S matrix for every one at once
    model.S(met_idx,1:end-1) = model.freeEnergy;

    %%New Part (4/30/2013)
    %Add water contribution, which isn't reflected elsewhere
    [~,rxn_idx] = intersect(model.rxns,water_rxn);
    model.S(met_idx,rxn_idx) = model.freeEnergy(rxn_idx);
    
    %Make the biomass value
    %Biomass Modification: -0.1764 kJ/GDW
    %Instead of using specific biomass ID, find it as the objective
    %rxn_idx = find(model.c~=0);
    %model.S(met_idx,(rxn_idx)) = -0.1764 + R*T*log(biomass);
    
    %%DEBUG CHECK: Print the dG metabolite row
    %model.S(met_idx,:)

    %Last step before simulating: restrict model coming out so that it must be < 0 and
    %there's no lower bound
    %For dynamic for now, don't set an upper bound, let dynamic code do that
    %model = changeRxnBounds(model,'GIBBS_kJ/GDW',0,'u');

    %W/O BIOMASS, MAKE UB -15 KJ/MOL, OR -0.015 KJ/MMOL
    %model = changeRxnBounds(model,'GIBBS_kJ/GDW',-0.015,'u');
    model = changeRxnBounds(model,'GIBBS_kJ/GDW',-inf,'l');

    %Simulate the model
    solution = optimizeCbModel(model,[],'one');

    %Find the gibbs flux
    %If no solution, return that!
    if isempty(solution.x)
        fprintf('\nNO THERMODYNAMICALLY FEASIBLE SOLUTION\n')
        gibbs_flux = [];
    else
        [~,idx] = intersect(model.rxns,'GIBBS_kJ/GDW');
        gibbs_flux = solution.x(idx);
    end

    %Final thing: return warnings for things produced/consumed that have no
    %substrate reaction specified
    %Find exchange reactions
    %exchanges = model.rxns(findExcRxns(model));

    %Find the indices of these exchanges where flux is non-zero
    %nonZero = find(solution.x(findExcRxns(model)));

    %Use to find the exchange reaction IDs for non-zero
    %nonZero_exchanges = exchanges(nonZero);

    %Find the ones in the non-Zero that aren't in substrate reactions
    %excluded = setdiff(nonZero_exchanges,substrateRxns);

    %Pull the GIBBS kJ_GDW out of excluded
    %excluded = setdiff(excluded,'GIBBS_kJ/GDW');

    %if ~isempty(excluded)

    %    warning('The following exchange reactions carry flux but are not accounted for in substrateRxns')
    %    for j = 1:length(excluded)
           %Print out the reactions that are problems
    %        fprintf('%s\n',excluded{j})
    %    end
    %end
end




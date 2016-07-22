function [solution,gibbs_flux] = optimizeThermoModel(model,substrateRxns,concentrations,biomass,T,water_rxn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Version 3: 04/24/2013


%GENERAL METHODOLOGY
%Accepts a model, substrate reactions, and initial concentrations/biomass
%Accepts a reaction for water
%Adds a metabolite called "dG[e]"
%Adds an exchange for dG called 'GIBBS' (EX_dG[e])
%Modifies exchanges to include dG[e] by adding a row
%Row values: adjusted for concentration (or biomass amount)
%Simulates the model
%Returns solution and overall dG
%(This is a modified addOverallDG code)

%%%%%%%%%%%%%%%%dG values are from Alberty, at pH=7.0%%%%%%%%%%%%%%%%%%
%Catch concentrations that are 0....
if any(~concentrations)
    solution = optimizeCbModel(model,[],'one');
    gibbs_flux = inf;
    
else
    
    %Error for things not the same size
    if length(substrateRxns) ~= length(concentrations)

        error('Please specify a concentration for each substrate reaction')
    end
    %Gas constant specification (Maybe move these out?)
    R = 8.314e-6; %kJ/mmol*K

    %Add the new reaction first, which adds the metabolite
    model = addReaction(model,'GIBBS_kJ/GDW','dG[e] <=> '); 

    %Find index of dG
    [~,met_idx] = intersect(model.mets,'dG[e]');


    %Intersect the substrate reactions with the model and make a dictionary
    [rxns,rxn_idx] = intersect(model.rxns,substrateRxns,'stable');
    dict = containers.Map(rxns,rxn_idx);

    %Loop: put in the correct free energy term for each:
    %dG = dG_0 + RTln(C)
    for i = 1:length(substrateRxns)

       %Change the dG weight for the exchange reaction (Conc in mM)
       model.S(met_idx,dict(substrateRxns{i})) = model.freeEnergy(dict(substrateRxns{i}))...
           +R*T*log(concentrations(i)/1000);    
    end
    
    %%New Part (4/30/2013)
    %Add water contribution
    [~,rxn_idx] = intersect(model.rxns,water_rxn);
    model.S(met_idx,rxn_idx) = model.freeEnergy(rxn_idx);
    
    %Make the biomass value
    %Biomass Modification: -0.1764 kJ/GDW
    %Instead of using specific biomass ID, find it as the objective
    biomass_idx = find(model.c~=0);
    model.S(met_idx,(biomass_idx)) = -0.1764 + R*T*log(biomass);
    
    %%DEBUG CHECK: Print the dG metabolite row
    %model.S(met_idx,:)

    %Last step before simulating: restrict model coming out so that it must be < 0 and
    %there's no lower bound
    %For dynamic for now, don't set an upper bound, let dynamic code do that
    %model = changeRxnBounds(model,'GIBBS_kJ/GDW',0,'u');

    %W/O BIOMASS, MAKE UB -15 KJ/MOL, OR -0.015 KJ/MMOL
    %model = changeRxnBounds(model,'GIBBS_kJ/GDW',-0.015,'u');
    model = changeRxnBounds(model,'GIBBS_kJ/GDW',-inf,'l');
    
    %%%%%%%%%%%%%%%%%
    %NEW LINE (THE KEY!)
    %%
    %Change objective to be the thermo reaction
    [~,dG_idx] = intersect(model.rxns,'GIBBS_kJ/GDW');
    model.c(biomass_idx) = 0;
    model.c(dG_idx) = 1;

    %Simulate the model (minimize dG!)
    solution = optimizeCbModel(model,'min');

    %Find the gibbs flux
    %If no solution, return that!
    if isempty(solution.x)
        fprintf('\nNO THERMODYNAMICALLY FEASIBLE SOLUTION\n')
        gibbs_flux = [];
    else
        gibbs_flux = solution.x(end);
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

    %    warning('The following exchange reactions carry flux but are not accounted for in subsrateRxns')
    %    for j = 1:length(excluded)
           %Print out the reactions that are problems
    %        fprintf('%s\n',excluded{j})
    %    end
    %end
end




function [solution,gibbs_flux,model] = optimizeThermoModel(model,substrateRxns,concentrations,T,water_rxn)

%%
% Version 4: 06/23/2015
%
% Accepts a model and necessary parameters for estimating thermodynamics of
% the overall system. Adds a new metabolite "dG" to the system that measures
% the free energy contribution for each exchange metabolite and a new
% reaction "GIBBS_kJ/GDW" that sums overall free energy for the system
%
% INPUT:
% model: a COBRA Toolbox model structure with a freeEnergy field
% substrateRxns: a set of exchange reactions for metabolites with specified
% concentrations
% concentrations: a set of concentrations corresponding to the specified
% substrateRxns (in mM)
% T: temperature (in Kelvin) for simulating growth
% water_rxn: identity of the water exchange reaction in the model. Water is
% treated separately from the aqueous metabolites and must be specified
% here.
%
% OUTPUT:
% solution: an FBA flux distribution that optimizes the supplied model
% gibbs_flux: an estimation of free energy produced by the model in the
% specified flux distribution (in kJ/gDCW/h)
% model: the supplied model with the addition of an overall reaction,
% GIBBS_kJ/GDW, that sums free energy for exchange reactions flowing in and
% out of the model
%
% Matthew Richards, 10/06/2015


% dG values are at pH=7.0 and ionic strength of 0.1 M
% Catch concentrations that are 0
if any(~concentrations)
    solution = optimizeCbModel(model,[],'one');
    gibbs_flux = inf;  
else
    
    % Give an error for things not the same size
    if length(substrateRxns) ~= length(concentrations)
        error('substrateRxns and concentrations must be of equal length')
    end
    
    % Gas constant specification
    R = 8.314e-6; %kJ/mmol*K

    % Add the new reaction first, which adds the metabolite
    model = addReaction(model,'GIBBS_kJ/GDW','dG <=> '); 

    % Find index of dG
    [~,met_idx] = intersect(model.mets,'dG');

    % Alter the free energy values for things with substrate reactions in
    % the free energy vector itself
    % First grab the index of the exchange reactions in the model
    [rxns,rxn_idx] = intersect(model.rxns,substrateRxns,'stable');
    % Make a dictionary
    dict = containers.Map(rxns,rxn_idx);      
    % For those indices, change the free energy numbers using concentration
    % Loop: put in the correct free energy term for each:
    % Basis: dG = dG_0 + RTln(C)
    for i = 1:length(substrateRxns)
       % Change the dG weight for the exchange reaction (Conc in mM)
       model.freeEnergy(dict(substrateRxns{i})) = model.freeEnergy(dict(substrateRxns{i}))...
           +R*T*log(concentrations(i));    
    end
    
    % Add free energy values to S matrix for every one at once
    model.S(met_idx,1:end-1) = model.freeEnergy;

    % New Part (4/30/2013)
    % Add water contribution, which isn't reflected elsewhere
    [~,rxn_idx] = intersect(model.rxns,water_rxn);
    model.S(met_idx,rxn_idx) = model.freeEnergy(rxn_idx);

    % Let the free energy be as low as it desires
    model = changeRxnBounds(model,'GIBBS_kJ/GDW',-inf,'l');

    % Simulate the model with minimization of overall flux
    solution = optimizeCbModel(model,[],'one');

    % Find the gibbs flux
    if isempty(solution.x)
        % If no solution, return that
        fprintf('\nNO THERMODYNAMICALLY FEASIBLE SOLUTION\n')
        gibbs_flux = [];
    else
        [~,idx] = intersect(model.rxns,'GIBBS_kJ/GDW');
        gibbs_flux = solution.x(idx);
    end
end
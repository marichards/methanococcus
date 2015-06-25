function [solution,gibbs_flux,thermo_model] = simMMGrowthThermo(model)

% Use addGD2MM to add thermodynamics to the model
thermo_model = addDG2MM(model);

% Simulate it with the optimizeThermoModel code

% Specify the exchanges we'll be using
substrateRxns = {...
    % Methane
    'Ex_cpd01024_c0';...
    % Hydrogen
    'Ex_cpd11640_c0';...
    % CO2(aq), not CO2(total)
    'EX_cpd00011_e0';...
    };

concentrations = [...
    % Methane
    1;...
    % Hydrogen
    1;...
    % CO2(total)
    1;...
    ];

[solution,gibbs_flux,thermo_model] = optimizeThermoModel(thermo_model,substrateRxns,concentrations,310,'EX_cpd00001_e0');




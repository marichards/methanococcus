function thermo_model = addDG2MM(model)

% Add the dG array specifically to the M. mariapludis model using values
% from the Equilibrator site

% Specify the exchanges we'll be using
exchanges = {...
    % Methane
    'EX_cpd01024_e0';...
    % Hydrogen
    'EX_cpd11640_e0';...
    % Water
    'EX_cpd00001_e0';...
    % CO2(total), not CO2(aq)
    'EX_cpd00011_e0';...
    };

% Add the free energy of formation (1 mM, pH=7, ionic strength =0.1 M) from
% the Equilibrator DB (in kJ/mmol)
dGs = [...
    % Methane
    0.1107;...
    % Hydrogen
    0.0816;...
    % Water
    -0.1576;...
    % CO2
    -0.4031;...
    ];

% Add error for dG just in case
dG_error = [...
    % Methane
    0.0058;...
    % Hydrogen
    0.0058;...
    % Water
    0.0016;...
    % CO2(aq)
    0.0058;...
    ];

% Add free energy to the model and create the new model
thermo_model = addDG(model,exchanges,dGs);




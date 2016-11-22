## iMR539, a Genome-Scale Metabolic Model of *Methanococcus maripaludis*

This repository contains the genome scale metabolic model of *Methanococcus maripaludis* and a collection of scripts for helping to simulate growth, perturb the model, and query various aspects of the model. The major directories are as follows:

### Model

Contains the most recent iteration of the *M. maripaludis* model, marked with the date it was created. It also contains the exact iMR539 model from publication in SBML and Matlab formats, the original automated model created using tools in Kbase, and a directory containing all versions of the model that came in between. 

### Matlab Scripts

Contains nearly all scripts created for manipulating and simulating growth with the *M. maripaludis* model. 

## Tutorial: Simulating Steady-State Growth

Once you have cloned this repository to your local machine, we recommend adding the ["Matlab Scripts"](./Matlab Scripts) directory and its sub-directories to your Matlab path, such that the appropriate scripts can be run from any directory rather than navigating back and forth. In the brief procedure that follows, we will assume that the ["Matlab Scripts"](./Matlab Scripts) are on your Matlab path. 

1) Navigate to the ["Model"](./Model) directory and load the model using the following command:

`load('iMR539.mat')`

2) Initialize the COBRA toolbox solver using the following command:

`initCobraToolbox`

Note that you should expect a pause of a few seconds as your machine initializes the solver

3) Simulate maximum steady state growth on H_{2} + CO_{2} using the following command:

`solution = maxGrowthOnH2(model);`

This should give you output that looks something like the following:

` Warning: Metabolite dG not in model - added to the model`
`> In addReaction at 213`
`  In optimizeThermoModel at 55`
`  In maxGrowthOnH2 at 51 `
`Warning: Metabolite name for dG set to dG `
`> In addReaction at 219`
`  In optimizeThermoModel at 55`
`  In maxGrowthOnH2 at 51 `
`Warning: Metabolite formula for dG set to '' `
`> In addReaction at 224`
`  In optimizeThermoModel at 55`
`  In maxGrowthOnH2 at 51 `
`GIBBS_kJ/GDW	dG 	<=>        `

`Biomass flux: 0.097315`

`Formate flux: 0.000000`
`CO2 flux: -45.257937`
`H2 flux: -190.527078`
`H2O flux: 93.716251`
`CH4 flux: 50.000000`
`NH3 flux: -0.756397`
`PO4 flux: -0.022829`
`Acetate flux: -3.666756`

`Overall reaction:`
`CO2 + 4 H2 --> 2 H2O + CH4`

`Model overall reaction (per mole CH4)`
`0.91 CO2 + 3.81 H2 --> 1.87 H2O + CH4`

`Predicted Yield Coefficient: 2.81 gDCW/mol CH4`

`Expected ATP/CH4 Yield: 0.5`
`Predicted ATP/CH4 Yield: 0.550`

`Warning: All external metabolite concentrations set to 1 mM `
`> In maxGrowthOnH2 at 99`
`Predicted Free Energy Generation: -5.592458 kJ/gDCW`

Don't be too alarmed by all the warnings generated here; they all have to do with additional free energy calculations. 

4) Print out the central methanogenesis fluxes using the following command:

`printMethanogenesisFluxes(model,solution,false);`

Here, we are passing the model, the solution we generated, and the parameter *false*, which tells the script not to generate a map of this sub-network. The map-drawing function relies upon the <a href="http://biosystems.lv/index.php/software/paint4net" target="_blank">Paint4Net toolbox</a>, a small set of functions that you can optionally add to your distribution of Matlab. We are not assuming that you have this toolbox installed; however, if you do then changing the *false* to *true* will generate a small metabolic map showing the fluxes through this sub-network. 

(**NOTE: The map drawing functionality relies upon flux variability analysis, which can often cause a numerical error for the first attempt. If you get an error beginning with `Index exceeds matrix dimensions` we recommend retrying the same command). 


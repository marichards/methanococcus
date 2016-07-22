## Matlab Scripts

This directory contains  all scripts created for manipulating and simulating growth with the *M. maripaludis* model. This directory is further divided into the following categories of scripts:

### Data Compilation for Paper

Contains scripts used to create tables and figures for the *M. maripaludis* model manuscript. These scripts are intended to provide a straightforward means to replicate data presented in the paper.

### Diagnostic Non-Map

Contains various scripts used to query metabolic models, not only the *M. maripaludis* model but also models of other species. Some of these scripts can optionally utilize map-drawing functions from the Paint4Net toolbox, but none of them absolutely require it. 

### Map Drawing

Contains several scripts that wrap around functions from the Paint4Net toolbox to quickly visualize flux maps from a metabolic model. As such, these functions explicitly depend on the Paint4Net toolbox found [here](http://biosystems.lv/index.php/software/paint4net).

### MM Growth Sims

Contains scripts used to simulate growth of the *M. maripaludis* model under a variety of different carbon and nitrogen sources. These scripts depend upon those found in the [Model Changing](./Model Changing) directory.

### Model Changing

Contains scripts that alter key parameters of the *M. maripaludis* model, namely the main carbon or nitrogen sources supplied in the *in silico* media formulation and the methane secretion rate that constrains the model's biomass production. 

### Strain Design

Contains a collection of scripts that simulate different metabolic engineering strategies for altering native *M. maripaludis* metabolism. These are primarily centered on the goal of inserting methanol consumption into *M. maripaludis* via insertion of a methanol methyl-CoM transferase, then reversing the pathway to achieve conversion of methane to methanol. 

### Thermodynamic Constraints

Contains scripts that add standard free energies of formation to a metabolic model and add free energy calculations to growth simulations with maximal biomass production. The "optimizeThermoModel.m" script is the central focus, at it builds around the "optimizeCbModel.m" script from the COBRA Toolbox by adding the free energy calculation. 

### Write to Other Formats

In addition to Matlab, there are multiple other platforms for working with constraint based models. This directory contains scripts that help convert the *M. maripaludis* model into alternative formats, principally including the ability to create tab-delimited files for upload to Kbase.

### Other

Contains a miscellaneous group of scripts that are unlikely to be used for standard simulations, but could potentially be useful in some way. 

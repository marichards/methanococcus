function [rxns,fluxes,directions] = metMassBalance(model,met,solution)

%Take in a model, a metabolite, and an FBA solution.  Find all the major
%sinks and sources of that metabolite by 

%Inputs
%model - a COBRA model structure
%met - a metabolite of interest
%solution - an FBA solution of the model

%Outputs
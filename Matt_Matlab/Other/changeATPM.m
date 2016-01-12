function model = changeATPM(model,gam,ngam)

% Given values for growth and non-growth associated ATP maintenance, change
% those values in the M. maripaludis model and return it in this altered
% state.
%
% INPUT
% model
% gam
% ngam
%
% OUTPUT
% model
%
% Matthew Richards 11/23/2015

% NGAM reaction is rxn00062[c0]
% Switch the bounds to the ngam value
model = changeRxnBounds(model,'rxn00062[c0]',ngam,'b');

% GAM is in biomass (biomass0 in this model)
% Find the index of the biomass reaction
[~,bio_idx] = intersect(model.rxns,'biomass0');

% Find the index of ATP in the biomass reaction
[~,atp_idx] = intersect(model.mets,'cpd00002[c0]');

% Change the index to the gam
model.S(atp_idx,bio_idx) = -gam;

% Also find the ADP, H+, water, and Phosphate indices
[~,adp_idx] = intersect(model.mets,'cpd00008[c0]');
[~,h_idx] = intersect(model.mets,'cpd00067[c0]');
[~,h2o_idx] = intersect(model.mets,'cpd00001[c0]');
[~,p_idx] = intersect(model.mets,'cpd00009[c0]');

% Change these indices to the gam too
model.S(adp_idx,bio_idx) = gam;
model.S(h_idx,bio_idx) = gam;
model.S(p_idx,bio_idx) = gam;
model.S(h2o_idx,bio_idx) = -gam;

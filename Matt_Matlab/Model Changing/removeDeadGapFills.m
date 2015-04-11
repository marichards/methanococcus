function [model,removed_rxns] = removeDeadGapFills(model)

% Take in a model, evaluate reactions without genes (gapfills) and dead end
% reactions to find reactions that fall under both headings, and remove
% those reactions from the model. 
%
% Inputs
% model: a COBRA model structure
%
% Outputs
% model: the input model with dead end, non-gene associated reactions
% removed
% removed_rxns: dead end, non-gene associated reactions that were removed
% from the model

% Step 1: Grab dead end reactions from the model as a subset
[~,~,rxns] = removeDeadEnds(model);

% Other Step 1: Run the reaction essentiality to find whether the reactions
% are singularly essential

% Step 2: Find non-gene associated reactions
non_gene_rxns = findRxnsWOGenes(model);

% Step 3: Intersect the two lists
removed_rxns = intersect(rxns,non_gene_rxns);

% Step 4: Weed out the transport and exchange reactions
trans_rxns = findTransRxnsMOD(model,false);
removed_rxns = setdiff(removed_rxns,trans_rxns);

% Step 5: Remove the reactions and return the model
model = removeRxns(model,removed_rxns);


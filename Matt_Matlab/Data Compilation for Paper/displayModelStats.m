function displayModelStats(model)

% Take in the M. maripaludis model and display the statistics I want for
% the paper

% First, take number of genes
num_genes = length(model.genes);

% Then find ORF coverage by dividing by 1722
orf_coverage = 100*num_genes/1722;

% Then find number of external metabolites
ext_mets = 0;
for i=1:length(model.mets)
if regexp(model.mets{i},'_e0')
ext_mets = ext_mets + 1;
end
end

% Use that to find internal mets
int_mets = length(model.mets)-ext_mets;

% Find dead end metabolites and reactions
tic
[~,removedMets,removedRxns] = removeDeadEnds(model);
toc
dead_mets = length(removedMets);
dead_rxns = length(removedRxns);

% Find exchange, transport, and internal reactions
exc_rxns = length(model.rxns(findExcRxns(model)));

trans_rxns = length(findTransRxnsMOD(model));

int_rxns = length(model.rxns) - exc_rxns - trans_rxns;

% Find reactions with genes
rxns_w_genes = length(model.rxns) - length(findRxnsWOGenes(model));

% Find internal reactions w/o genes
rxns_wo_genes = length(setdiff(findRxnsWOGenes(model),findTransRxnsMOD(model)));

%Print out all the results
fprintf('\n\nMethanococcus maripaludis S2 Model Statistics\n')
fprintf('---------------------------------------------------\n')
fprintf('Protein Coding Genes: %d\n',num_genes);
fprintf('%% ORF Coverage : %0.2f\n',orf_coverage);
fprintf('Intra/Extracellular Metabolites: %d/%d\n',int_mets,ext_mets);
fprintf('Dead End Metabolites: %d\n',dead_mets);
fprintf('Internal Reactions: %d\n',int_rxns);
fprintf('Transport/Exchange Reactions: %d/%d\n',trans_rxns,exc_rxns);
fprintf('Gene-Associated Reactions: %d\n',rxns_w_genes);
fprintf('Dead End Reactions: %d\n',dead_rxns);
fprintf('Gapfilled Internal Reactions: %d\n',rxns_wo_genes)
fprintf('%% Reactions with Genes: %0.0f\n',100*rxns_w_genes/(trans_rxns+int_rxns))
function [stats,rxns] = classifyGapRxns(model)

% For all of the gapfilling reactions, pull out the subsystems associated
% with them and display the results

% First, pull out the set of reaction without genes
rxns = setdiff(setdiff(findRxnsWOGenes(model),findTransRxnsMOD(model)),...
    model.rxns(findExcRxns(model)));

% Now find those reactions in the model
[rxns,idx] = intersect(model.rxns,rxns);

% Pull out the subsystems
subs = model.subSystems(idx);


stats = tabulate(subs);
%figure(1)
%bar(cell2mat(tbl(:,2)))
%set(gca,'XTickLabel',tbl(:,1))


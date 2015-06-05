function rxn_list = nadRxns(model)
%Outputs a text file of reactions, 
%plus a cell with the reactions(rxn_list)

%List out all the cofactors I want
metList={'NADP_c0','NAD_c0','NADPH_c0','NADH_c0'};


%This function will find all the rxns with nad/nadp
[rxnList, ~] = findRxnsFromMets(model,metList);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is the perfect place to remove any reactions that are exactly
%identical




%Pull out the reactions in the S-matrix form
%Find the indices in the matrix for each reaction

[~,inds] = intersect(model.rxns,rxnList);

%First remove the cofactors
no_cof = removeMetabolites(model,metList);

%Select the piece of S with my indices
%Transpose for duplicate-finding later
my_S = no_cof.S(:,inds)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Find the duplicates!
%Use inds again...because I can!
[~,inds] = unique(my_S,'rows');
dups=[];
for i = 1:size(my_S,1)
    if ~ismember(i,inds)
        dups = [dups,i];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Have duplicates....want to ouput 
%Initialize cell for the rxns
rxn_list = repmat({''},2,1);
counter=1;
%For each of the indices of the duplicates...
for i=1:length(dups)
    %Pull out the row (rxn)
    my_row = my_S(dups(i),:);
    %Find all occurances of that row in the my_S matrix
    %These are the reaction numbers
    all_inds = find(ismember(my_S,my_row,'rows'));
    %all_inds has indices for all occurances of that row
    %Write out the rxnFormulaList entry for each!
    for j = 1:length(all_inds)
        %Add these to the list of reactions
        rxn_list(counter)=rxnList(all_inds(j));
        counter = counter+1;
    end
    
end

%Shold do it with intersect instead
%all_rxns = intersect(

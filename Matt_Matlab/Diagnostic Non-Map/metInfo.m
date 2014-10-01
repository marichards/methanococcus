function locationInS=metInfo(cbModel,indices)
% Given metabolite indices in vector form, return the name of the 
% metabolite and the names of the rxns it is involved in.
% Inputs are COBRA models and a vector of metabolite indices (according to the S-matrix). 

% For each of the user-inputted metabolite indices, output the name of the 
% metabolite and the names of reactions it is involved in.

%Originally by HA
%modified by bh 4/11/11 to display model.met  and model.rxn as well

jnew=0; 
for i=1:length(indices)
    fprintf('********************** For metabolite index %u **********************',indices(i));
    metName = cbModel.metNames(indices(i));
    fprintf('\nmetName: %s \n', char(metName{:}));
    met = cbModel.mets(indices(i));
    fprintf('met: %s \n', char(met{:}));
    
    
    fprintf('\nName of reactions that it is involved in: \n');
    involvedRxns = find(cbModel.S(indices(i),:));
    for j=1:length(involvedRxns)
%        rxnName = cbModel.rxnNames(involvedRxns(j));
%        fprintf('%u: %s (rxn index %u)\n', j, char(rxnName{:}), involvedRxns(j));
         rxnName = cbModel.rxnNames(involvedRxns(j));
         rxn = cbModel.rxns(involvedRxns(j));
         fprintf('%u: %s (rxn %u, %s)\n', j, ... 
             char(rxnName{:}), involvedRxns(j), char(rxn{:}));
       
        %Use code below if you want to get the met index and rxn index as a
        %ordered pair i.e. [met index, rxn index]
        locationInS(jnew+j,:)=[indices(i),involvedRxns(j)];
    end
    jnew=jnew+j;
end
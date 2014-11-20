function confidences = assembleRxnLikelihoods(model,reactions,probability1)

%Run this script in the Reaction Likelihood Folder!

%Compare the model reactions to the likelihood reactions:
[rxns,model_idx,conf_idx] = intersect(model.rxns,reactions);

%Now create an array for these reactions. First create an empty array the
%length of model reactions
confidences = zeros(length(model.rxns),1);

%Run through each reaction 
for i=1:length(rxns)
    %Assign to the confidenc eof that reaction the value in the probability
    %array for that reaction
    confidences(model_idx(i))=probability1(conf_idx(i));
    
end




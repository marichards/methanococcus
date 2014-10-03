function model = addGenesToModel(model)

%Add the genes and gprs to the model we're comparing to
load('gprs_no523_no525.mat')

%Loop through the reactions, add the appropriate grRules for each
for i = 1:length(model.rxns)
    
    %Add the gene association rule and make it lower case
    model = changeGeneAssociation(model,model.rxns{i},lower(gprs{i}));
end

%Add the other reactions too
model = addReaction(model,'R523',...
    'Fe3_e + atp + h2o <=> fe3 + adp + pi');
model = addReaction(model,'R525',...
    'Fe3_e <=> fe3');

%Add their genes
model = changeGeneAssociation(model,'R523','(MMP0198 or MMP1183)');
model = changeGeneAssociation(model,'R525','MMP0197');

    
function compareGeneEIs(model)

%Compare our model prediction to what they have

%Load the EIs
load('2014_11_04_EIs.mat')

%Pull out the number of genes in common, this is all we need to
%compare to
common_genes = intersect(model.genes,EIs.all);

%Then run gene essentiality
grRatio = singleGeneDeletion(model,'FBA',common_genes);
%Find the essential ones
essentials = common_genes(grRatio<0.1);
non_essentials = setdiff(common_genes,essentials);


%Then compare each set of essentials, non-essentials, and then MCC it

%MCC = [(TP*TN)-(FP*FN)]/sqrt[(TP+FP)(TP+FN)(TN+FP)(TN+FN)]

%TP: Is non-essential, Predicted non-essential (length(intersect(
%TN: Is essential, Predicted essential
%FP: Is essential, Predicted non-essential
%FN: Is non-essential, Predicted essential (

fprintf('\nStatistics for All Cases\n\n')
fprintf('Cases\tTP\tTN\tFP\tFN\tMCC\t\n')
%Compare the all 4
TP = length(intersect(non_essentials,setdiff(common_genes,EIs.four)));
TN = length(intersect(essentials,EIs.four));
FP = length(intersect(non_essentials,EIs.four));
FN = length(intersect(essentials,setdiff(common_genes,EIs.four)));

MCC = (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
fprintf('4\t\t%d\t%d\t%d\t%d\t%0.3f\n',TP,TN,FP,FN,MCC);

%Compare the 3
TP = length(intersect(non_essentials,setdiff(common_genes,EIs.three)));
TN = length(intersect(essentials,EIs.three));
FP = length(intersect(non_essentials,EIs.three));
FN = length(intersect(essentials,setdiff(common_genes,EIs.three)));

MCC = (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
fprintf('3\t\t%d\t%d\t%d\t%d\t%0.3f\n',TP,TN,FP,FN,MCC);

%Compare the 2
TP = length(intersect(non_essentials,setdiff(common_genes,EIs.two)));
TN = length(intersect(essentials,EIs.two));
FP = length(intersect(non_essentials,EIs.two));
FN = length(intersect(essentials,setdiff(common_genes,EIs.two)));

MCC = (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
fprintf('2\t\t%d\t%d\t%d\t%d\t%0.3f\n',TP,TN,FP,FN,MCC);

%Compare the 1
TP = length(intersect(non_essentials,setdiff(common_genes,EIs.one)));
TN = length(intersect(essentials,EIs.one));
FP = length(intersect(non_essentials,EIs.one));
FN = length(intersect(essentials,setdiff(common_genes,EIs.one)));

MCC = (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
fprintf('1\t\t%d\t%d\t%d\t%d\t%0.3f\n',TP,TN,FP,FN,MCC);


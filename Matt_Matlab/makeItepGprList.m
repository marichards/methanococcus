function makeItepGprList(model,filename)
%
%Take in a COBRA model file and print out a text file with two columns, the
%first listing every reaction ID, the second listing the GPR rule
%
%Open the text file
file_id = fopen(sprintf('%s.txt',filename),'w');

%%Loop through the model reactions
for i =1:length(model.rxns)
   
    %If there are genes
    if ~isempty(model.grRules{i})
        
        %Then print them to the file with the reaction
        fprintf(file_id,'%s\t%s\n',model.rxns{i},model.grRules{i});
    end
    

end
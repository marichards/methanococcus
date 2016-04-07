function model = changeGenesToKbase(model)

% First change the genes themselves
for i = 1:length(model.genes)
    %Catch 3-zero cases
    if regexp(model.genes{i},'mmp000')
        model.genes{i} = regexprep(model.genes{i},'mmp000','kb\|g.575.peg.');
    end
    %Catch 2-zero cases
    if regexp(model.genes{i},'mmp00')
        model.genes{i} = regexprep(model.genes{i},'mmp00','kb\|g.575.peg.');
    end
    %Catch 1-zero cases
    if regexp(model.genes{i},'mmp0')
        model.genes{i} = regexprep(model.genes{i},'mmp0','kb\|g.575.peg.');
    end    
          
	%Catch no-zero cases
	model.genes{i} = regexprep(model.genes{i},'mmp','kb\|g.575.peg.');
end

% Do the same for grRules
for i = 1:length(model.grRules)
    %Catch 3-zero cases
    if regexp(model.grRules{i},'mmp000')
        model.grRules{i} = regexprep(model.grRules{i},'mmp000','kb\|g.575.peg.');
    end
    %Catch 2-zero cases
    if regexp(model.grRules{i},'mmp00')
        model.grRules{i} = regexprep(model.grRules{i},'mmp00','kb\|g.575.peg.');
    end
    %Catch 1-zero cases
    if regexp(model.grRules{i},'mmp0')
        model.grRules{i} = regexprep(model.grRules{i},'mmp0','kb\|g.575.peg.');
    end              
	%Catch no-zero cases
	model.grRules{i} = regexprep(model.grRules{i},'mmp','kb\|g.575.peg.');
end
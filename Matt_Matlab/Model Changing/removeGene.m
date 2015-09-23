function model = removeGene(model,gene)

%Remove gene from the model

%First find the index
[~,idx] = intersect(model.genes,gene);

%Remove the column from the gene-reaction matrix
model.rxnGeneMat(:,idx)=[];

%Remove it from the rules...format is x(idx)
%There are | and & characters.  Must remove the gene first, then remove
%useless characters.  Else we'll remove too many

%Loop it
for i=1:length(model.rules)

    %First remove the gene reference itself
    model.rules{i} = regexprep(model.rules{i},sprintf(' *x\\(%d\\) *',idx),'');

    %Now we may have several things: ||, &&, |&, &|, and leading/trailing of
    %either.  Plus, there can be parentheses before or after any of
    %those...
    
    %Remove doubles (and more i suppose)
    model.rules{i} = regexprep(model.rules{i},'\|+','\|');
    model.rules{i} = regexprep(model.rules{i},'\&+','\&');
    
    %Remove mixed doubles with &...and A|B & C|D will still be an "and"
    model.rules{i} = regexprep(model.rules{i},'\&\|','\&');
    model.rules{i} = regexprep(model.rules{i},'\|\&','\&');
    
    %Remove starting or trailing | or & symbols
    model.rules{i} = regexprep(model.rules{i},'^\|','');
    model.rules{i} = regexprep(model.rules{i},'\|$','');
    model.rules{i} = regexprep(model.rules{i},'^\&','');
    model.rules{i} = regexprep(model.rules{i},'\&$','');
    
    %Include possibility of a parentheses    
    model.rules{i} = regexprep(model.rules{i},'\(\|','\(');
    model.rules{i} = regexprep(model.rules{i},'\|\)','\)');
    model.rules{i} = regexprep(model.rules{i},'\(\&','\(');
    model.rules{i} = regexprep(model.rules{i},'\&\)','\)');    
    
    %Finally, remove empty parentheses
    model.rules{i} = regexprep(model.rules{i},'\(\)','');


end

%Do the same removal, but for the gene from the grRules
%Must use the actual gene name there
%Have "and" and "or" instead of escape characters

%Loop it
for i=1:length(model.grRules)
    %First remove the gene itself
    model.grRules{i} = regexprep(model.grRules{i},sprintf(' *%s *',escapeString(gene)),'');
    
    %Remove multiples with just one
    model.grRules{i} = regexprep(model.grRules{i},'(and)+','and');
    model.grRules{i} = regexprep(model.grRules{i},'(or)+','or');
    model.grRules{i} = regexprep(model.grRules{i},'\|+','\|');
    model.grRules{i} = regexprep(model.grRules{i},'\&+','\&');
    
    %Remove mixed doubles with 'and'
    model.grRules{i} = regexprep(model.grRules{i},'andor','and');
    model.grRules{i} = regexprep(model.grRules{i},'orand','and');
    model.grRules{i} = regexprep(model.grRules{i},'\&or','and');
    model.grRules{i} = regexprep(model.grRules{i},'or\&','and');
    model.grRules{i} = regexprep(model.grRules{i},'and\|','and ');
    model.grRules{i} = regexprep(model.grRules{i},'\|and',' and');
    model.grRules{i} = regexprep(model.grRules{i},'\&\|',' and ');
    model.grRules{i} = regexprep(model.grRules{i},'\|\&',' and ');
    model.grRules{i} = regexprep(model.grRules{i},'and\&','and ');
    model.grRules{i} = regexprep(model.grRules{i},'\&and',' and');
    model.grRules{i} = regexprep(model.grRules{i},'\|or',' or');
    model.grRules{i} = regexprep(model.grRules{i},'or\|','or ');
    
    
    %Remove starting or trailing qualifiers
    model.grRules{i} = regexprep(model.grRules{i},'^(and) ?','');
    model.grRules{i} = regexprep(model.grRules{i},' ?(and)$','');
    model.grRules{i} = regexprep(model.grRules{i},'^(or) ?','');
    model.grRules{i} = regexprep(model.grRules{i},' ?(or)$',''); 
    model.grRules{i} = regexprep(model.grRules{i},'^\|','');
    model.grRules{i} = regexprep(model.grRules{i},'\|$','');
    model.grRules{i} = regexprep(model.grRules{i},'^\&','');
    model.grRules{i} = regexprep(model.grRules{i},'\&$','');
    
    %Include possibility of a parentheses  
    model.grRules{i} = regexprep(model.grRules{i},'\((and) ?','\(');
    model.grRules{i} = regexprep(model.grRules{i},' ?(and)\)','\)');
    model.grRules{i} = regexprep(model.grRules{i},'\((or) ?','\(');
    model.grRules{i} = regexprep(model.grRules{i},' ?(or)\)','\)'); 
    model.grRules{i} = regexprep(model.grRules{i},'\(\|','\(');
    model.grRules{i} = regexprep(model.grRules{i},'\|\)','\)');
    model.grRules{i} = regexprep(model.grRules{i},'\(\&','\(');
    model.grRules{i} = regexprep(model.grRules{i},'\&\)','\)');    
    
    %Finally, remove empty parentheses
    model.grRules{i} = regexprep(model.grRules{i},'\(\)','');
end


%Lastly, remove gene from model.genes
model.genes(idx)=[];

%% New addition on 9/22/2015: fix the indexing on model.rules that is caused by removing the gene
% Must go through all reactions and make the rule match the grRule

% Go over each grRule (through each reaction)
for j = 1:length(model.grRules)
    % Pull out the genes and rule encodes by grRules
    [genes,rule] = parseBoolean(model.grRules{j});
    % For each gene in the rule
    for i = 1:length(genes)
        % Find the gene index in the model
        geneID = find(strcmp(model.genes,genes{i}));
        % Set the index to 1 to show it's part of the rxn (index of j)
        model.rxnGeneMat(j,geneID) = 1;
        % Convert the rule using the indices of the genes
        rule = strrep(rule,['x(' num2str(i) ')'],['x(' num2str(geneID) ')']);;
    end
    model.rules{j} = rule;
end
    
    
    

function createKbaseModel(model,model_filename,media_filename,compound_filename)
%%
% Takes in a COBRA model and converts it to a TSV file compatible with the
% Kbase Narrative import function. Outputs a model (reactions) file, a
% media file, and a compounds file that can be separately imported. 
%
% INPUT
% model: a COBRA Toolbox model structure
% model_filename: a supplied name for the model (or reactions) file created
% for upload to Kbase
% media_filename: a supplied name for the media file created for upload to
% Kbase
% compound_filename: a supplied name for the compound file created for
% upload to Kbase
%
%Written by Matthew Richards 01/21/2015
%Edited by MR 03/31/2015 to fix media file and create compound file

%%
% Create 3 files:
% model_file: a tab-delimited file the contains all reactions, their
% directions, and their formulae. The biomass equation comes first.
% Anything that doesn't fit a category gets a blank space. Format is as
% follows: 
% id name direction gpr pathway_enzyme equation reference
%
% media_file: a tab-delimited file that contains all media compounds, their
% concentrations (0.01), their lb, and their ub
% compid conc lb ub
%
% compound_file: a tab-delimited file that contains all model compounds and
% their information encoded in the model with specified format:
% id	name	aliases	charge	formula
% Aliases: KEGG:C00001;BiGG:h2o

%%
%Create the model file
model_file_id = fopen(sprintf('%s.txt',model_filename),'w');
%Print the headers
fprintf(model_file_id,'id\tname\tdirection\tgpr\tpathway\tenzyme\tequation\treference\n');
%Create a full index of reactions
rxn_idx = (1:length(model.rxns))';
%Find the biomass equation index using the objective
bio_idx = find(model.c);
%Take out the biomass index
rxn_idx = setdiff(rxn_idx,bio_idx,'stable');
%Add the biomass index back to the top
rxn_idx = [bio_idx;rxn_idx];

%Now we have the correct index order; create the correct fields
%Loop through the index
for i=1:length(rxn_idx)
    %Make sure to use the number in the index for grabbing each reaction
    %First grab the reaction formula
    formula = printRxnFormula(model,model.rxns{rxn_idx(i)},false);
    %Set a direction based on the reversibility
    %Reversible
    if model.ub(rxn_idx(i))>0 && model.lb(rxn_idx(i))<0
        direction = '=';
    %Forward
    elseif model.ub(rxn_idx(i))>0 && model.lb(rxn_idx(i))>=0
        direction = '>';
    %Reverse        
    elseif model.ub(rxn_idx(i))<=0 && model.lb(rxn_idx(i))<0
        direction = '<';
    %Only remaining case should be both zeros
    else
        direction = '=';
    end
    
    %Now we have formula and direction, still need a gene properly formatted
    %Remove "mmp" and replace it with the 'kb|g.575.peg.'
    genes = model.grRules{rxn_idx(i)};
    %Catch 3-zero cases
    if regexp(genes,'mmp000')
        genes = regexprep(genes,'mmp000','kb\|g.575.peg.');
    end
    %Catch 2-zero cases
    if regexp(genes,'mmp00')
        genes = regexprep(genes,'mmp00','kb\|g.575.peg.');
    end
    %Catch 1-zero cases
    if regexp(genes,'mmp0')
        genes = regexprep(genes,'mmp0','kb\|g.575.peg.');
    end    
          
	%Catch no-zero cases
	genes = regexprep(genes,'mmp','kb\|g.575.peg.');
    fprintf(model_file_id,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t \n',...
        model.rxns{rxn_idx(i)},model.rxnNames{rxn_idx(i)},direction,genes,model.subSystems{rxn_idx(i)},model.rxnECNumbers{rxn_idx(i)},formula{1});
end


%Loop through the reactions and write them to the file
%Change from previous: switch UB and LB (upper bound first) and take negative of each [as per Kbase format, where "upper bound" is how much the model can take up] 

% Create the media file, which is simpler
media_file_id = fopen(sprintf('%s.txt',media_filename),'w');
% Grab all exchange reactions
exc_rxns = model.rxns(findExcRxns(model));
% Check if dG is in the model
if ismember('cpd15000[c0]',model.mets)
    % If so, then grab reactions with it
    rxns = findRxnsFromMets(model,'cpd15000[c0]');
    % And add it to the exc_rxns
    exc_rxns = unique([exc_rxns;rxns]);
end
% Intersect them to get their indices
[exc_rxns,idx] = intersect(model.rxns,exc_rxns);
% Find the gibbs index
[~,gibbs_idx] = intersect(exc_rxns,'GIBBS_kJ/GDW');
% Loop through them
for i=1:length(exc_rxns)
    %Find the metabolites in that reaction and remove dG if necessary
    % But don't remove it if it's GIBBS_KJ
    if i == gibbs_idx
        met = {'cpd15000[c0]'};
    else
        met = setdiff(findMetsFromRxns(model,exc_rxns{i}),'cpd15000[c0]');       
    end
    %For each exchange reaction, print the 4 things needed
    fprintf(media_file_id,'%s\t0.01\t%f\t%f\n',met{1},-model.ub(idx(i)),-model.lb(idx(i)));   
end

%Create the compounds file, with specified format
%id	name	aliases	charge	formula
%Aliases: KEGG:C00001;BiGG:h2o
compound_file_id = fopen(sprintf('%s.txt',compound_filename),'w');
%Write the header
fprintf(compound_file_id,'id\tname\taliases\tcharge\tformula\n');
%Loop through the compounds and write them to the file accordingly
for i=1:length(model.mets)
    %Print out each part of the compounds
    fprintf(compound_file_id,'%s\t%s\tKEGG:%s\t%d\t%s\n',...
       model.mets{i},model.metNames{i},model.metKEGGID{i},model.metCharge(i),model.metFormulas{i}); 
end
function createGAMSModel(model,model_filename)
%%
% Takes in COBRA model, and writes it to a GAMS-readable format, including
% the S-matrix plus all metabolites, reactions, and genes. 
%
% INPUT
% model: a COBRA Toolbox model structure
% model_filename: a filename for the resulting GAMS-compatible model file
%
% Matthew Richards, 09/29/2015


%% GENERAL PLAN
% Write metabolites and reactions as comma-delimited lists

% Add max and min (Vmax, Vmin) fluxes for each reaction

% Tie them together using the S-matrix, new-line delimited
% Format: reaction.metabolite\tcoefficient
%%

% First, create the file and label the top
file_id = fopen(sprintf('%s.txt',model_filename),'w');

fprintf(file_id,'%s for use in SimOptStrain\n\n',model_filename);

% Print metabolites set
fprintf(file_id,'Sets\ni metabolites in S (m)\n/\n');

% 4/21/2015
% Alter mets so it prints the index, not the ID
for i=1:length(model.mets)
    fprintf(file_id,'CPD%0.3d\n',i);
end

% More formatting
fprintf(file_id,'/\nj reactions in S (n)\n/\n');

% Print reactions set
for j=1:length(model.rxns)
    fprintf(file_id,'%s\n',model.rxns{j});
end

% More formatting; No upper/lower limits yet
fprintf(file_id,'/\nParameters\n\n');
% Do upper and lower bounds how Vangelis specified
% Print each reaction and its upper bound
fprintf(file_id,'UpperLimits(j) maximum value flux can take\n/\n');
for j=1:length(model.rxns)
    fprintf(file_id,'%s\t%0.2f\n',model.rxns{j},model.ub(j));
end
% Print each reaction and its lower bound
fprintf(file_id,'/;\n\nLowerLimits(j) minimum value flux can take\n/\n');
for j=1:length(model.rxns)
    fprintf(file_id,'%s\t%0.2f\n',model.rxns{j},model.lb(j));
end
fprintf(file_id,'/;\n\nS(i,j) contains the S matrix\n/\n');

% Make non-sparse matrix A
A = full(model.S);
% Use S matrix to print out everything for each reaction
for i = 1:length(model.mets)
    % Grab index of mets
    idx = find(model.S(i,:));
    % For each thing in the index,
    for j=1:length(idx)
        %Print out the rxn.met and then the coefficient
        % 4/21/2015
        % Alter mets so it prints the index, not the ID
        fprintf(file_id,'CPD%0.3d.%s\t%f\n',i,model.rxns{idx(j)},A(i,idx(j)));
    end
end

fprintf(file_id,'/;\n\nGene-Reaction Relationships\n/\n');

% Print out GPRs in reaction\tgpr format:
for i =1:length(model.rxns)
fprintf(file_id,'%s\t%s\n',model.rxns{i},model.grRules{i});
end
fprintf(file_id,'/;');
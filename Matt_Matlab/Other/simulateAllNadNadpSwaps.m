function [no,yes,yes_gw] = simulateAllNadNadpSwaps(cobra)
%
% function simulateAllNadNadpSwaps(cobra)
%
% Simulate what happens if we pick one of NAD and NADP
% for each of the ambiguous rxns and then run FBA

ambiguousPairs = ...
    {'rxn00248_c0', 'rxn00249_c0'; ...
    'rxn01268_c0', 'rxn01269_c0'; ...
    'rxn01301_c0', 'rxn01302_c0'; ...
    'rxn01313_c0', 'rxn01314_c0'; ...
    'rxn05117_c0', 'rxn05119_c0'...
    };

% Turn off ATPM for this study
%cobra = changeRxnBounds(cobra, {'ATPM'}, 0, 'l');



% Benchmark
%fprintf(1, 'Wild type + ACALD + ALCD2x + F4NPR (no NADNADP)\n');
%simGrowth(cobra);

% Set up a list of indexes to use (nifty trick - thanks Andrew)
x = 0:(2^(size(ambiguousPairs, 1))-1);
x = dec2bin(x);

no = {};
yes = {};
yes_gw = [];
tic
for ii=1:size(x,1)
    newmodel = cobra;
    la = binStringToLogicArray(x(ii,:));
    % Array of reactions to turn off
    arr = [ ambiguousPairs(la, 1); ambiguousPairs(~la, 2) ];
    [~, modelidx] = intersect(newmodel.rxns, arr);
    newmodel.lb(modelidx) = 0;
    newmodel.ub(modelidx) = 0;
    %fprintf(1, 'String: %s\n', x(ii,:));
   % fprintf(1, 'NO NADNADP:\n')
    flag = simGrowth(newmodel);
    if flag > 0
        yes = [yes ; x(ii,:)];
        yes_gw = [yes_gw; flag];
    else
        no = [no ; x(ii,:)];
    end
    continue;
end
toc
%%%%%%%%
function growth_flag = simGrowth(cobra)
%Store the string

%methaneloc = findRxnIDs(cobra, 'EX_ch4(e)');
cobra = changeObjective(cobra, 'biomass0');
f = optimizeCbModel(cobra,[],'one');
%fprintf(1, 'Growth rate: %1.5f\n', f.f);
if f.f < 0.01
    growth_flag = 0;
else
   growth_flag = f.f;
end
%cobra = changeObjective(cobra, 'EX_ch4(e)');
%f = optimizeCbModel(cobra);
%fprintf(1, 'Max methane potential: %1.5f\n', f.f);

%%%%%%
function r = binStringToLogicArray(binstring)
% I have Andrew to thank for this one too... can I do anything myself?
% Don't answer that.
r = logical(binstring - '0');
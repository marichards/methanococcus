function quantifyNADSwapGR(model)

%Take in the M.maripaludis model, swap out the NAD/NADP reactions to see
%all cases, then create a plot showing which ones have the highest average
%growth rates

% Simulate what happens if we pick one of NAD and NADP
% for each of the ambiguous rxns and then run FBA

ambiguousPairs = ...
    {'rxn00248_c0', 'rxn00249_c0'; ...
    'rxn01268_c0', 'rxn01269_c0'; ...
    'rxn01301_c0', 'rxn01302_c0'; ...
    'rxn01313_c0', 'rxn01314_c0'; ...
    'rxn05117_c0', 'rxn05119_c0'...
    };


% Benchmark
%fprintf(1, 'Wild type + ACALD + ALCD2x + F4NPR (no NADNADP)\n');
wt_yield = simGrowth(model);

% Set up a list of indexes to use (nifty trick - thanks Andrew)
x = 0:(2^(size(ambiguousPairs, 1))-1);
x = dec2bin(x);

no = {};
yes = {};
yes_gw = [];
%Make a set of reactions that are turned off each time
off_rxns = cell(size(x));
tic
for ii=1:size(x,1)
    newmodel = model;
    la = binStringToLogicArray(x(ii,:));
    % Array of reactions to turn off
    arr = [ ambiguousPairs(la, 1); ambiguousPairs(~la, 2) ];
    %Save them to the turned-off reactions
    off_rxns(ii,:)=arr';
    [~, modelidx] = intersect(newmodel.rxns, arr);
    newmodel.lb(modelidx) = 0;
    newmodel.ub(modelidx) = 0;
    %fprintf(1, 'String: %s\n', x(ii,:));
    growth_yield = simGrowth(newmodel);
    if growth_yield > 0
        yes = [yes ; x(ii,:)];
        yes_gw = [yes_gw; growth_yield];
    else
        no = [no ; x(ii,:)];
    end
    continue;
end
toc

%I'm left with "yes", an array of growing things with the binary
%combination, and "yes_gw", an array of their growth yields

%Now I also have off_rxns, a 32x5 cell array where each row contains the
%reactions turned off there.

%For each of the 10 reactions, go through each row and if it's there, then
%pull out the growth rate

%First create something to store rates in. It should be 16/10 (half the
%total conditions, twice the total offs)
%Grab the size of the "off_rxns", which is 32 x 5
[m,n]=size(off_rxns);
on_yields = zeros(m/2,2*n);

%Now for each thing in ambiguous pairs, make a corresponding row
%First make an array that's the ambiguous pairs in a row
all_pairs = [ambiguousPairs(:,1);ambiguousPairs(:,2)];
%Now for each one, loop through, find it, and pull things out
%Loop through the reactions
for j = 1:length(all_pairs)
    %Start at row 1 of on_yields
    row_num = 1;
    %Loop through the off_rxns
    for k = 1:m
        %If it's not there
        if ~ismember(all_pairs{j},off_rxns(k,:))
            %That means it was on; save its growth yield and add to the
            %counter
            on_yields(row_num,j)=yes_gw(k);
            row_num=row_num+1;
        end
    end
end

%So now I should have a 16x10 array called "on_yields" with all the growth
%yiels and "all_pairs", a cell array that labels what reaction corresponds
%to each column

%Take the average of each column (can take standard dev too...)
avg_yield = mean(on_yields,1);
std_yield = std(on_yields,1);

figure(1)
bar(avg_yield,'r')
set(gca,'XTickLabel',all_pairs)
ylabel('Average Yield When On')
title ('Growth Yields For NAD/NADP')

end

%%%%%%%%
function growth_yield = simGrowth(model)
%Store the string
%methaneloc = findRxnIDs(model, 'EX_ch4(e)');
model = changeObjective(model, 'biomass0');
f = optimizeCbModel(model,[],'one');
%fprintf(1, 'Growth rate: %1.5f\
%Find the growth yield
[~,ch4_idx] = intersect(model.rxns,'Ex_cpd01024_c0');
growth_yield = f.f*1000/f.x(ch4_idx);
if growth_yield < 0.01
    growth_yield = 0;
else
   growth_yield = growth_yield;
end
%model = changeObjective(model, 'EX_ch4(e)');
%f = optimizeCbModel(model);
%fprintf(1, 'Max methane potential: %1.5f\n', f.f);

end
%%%%%%
function r = binStringToLogicArray(binstring)
% I have Andrew to thank for this one too... can I do anything myself?
% Don't answer that.
r = logical(binstring - '0');
end
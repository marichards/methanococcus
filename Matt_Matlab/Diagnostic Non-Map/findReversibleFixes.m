function fixRxns = findReversibleFixes(model,rxn)

%Given an essential reversible reaction for a model, go through the irreverisble
%internal reactions of the model and find all those that, when they are
%made reversible, allow model growth when the essential reaction is
%constrained
%
%Inputs
%
%Outputs
tic
%First, simulate the model to ensure it grows
solution = optimizeCbModel(model,[],'one');

if solution.f<=1e-8
    error('Please supply a model that grows')
end

%Now simulate it with the reaction turned forward
fwd_model = changeRxnBounds(model,rxn,0,'l');
fwd_solution = optimizeCbModel(fwd_model,[],'one');

%Find the solution that matches KO by testing forward and make the model
%the INCORRECT direction going forward
if fwd_solution.f <= 1e-8 
    model = fwd_model;
elseif fwd_solution.f == solution.f
    model = changeRxnBounds(model,rxn,0,'u');
else
    error('Making the supplied reaction irreversible is not lethal')
end

%Now we have a model that shouldn't grow
%Get a list of all irreversible reactions
rxns = model.rxns(model.rev==0);
 
%Find exchanges and setdiff with them
rxns = setdiff(rxns,model.rxns(findExcRxns(model)));

%Create an array to store reaction indices
fixRxns = {};
%Cycle through the reactions and try turning each reversible
for i=1:length(rxns)
    %For each of these reactions, turn it reversible (-1000-1000)
    test = changeRxnBounds(model,rxns{i},-1000,'l');
    test = changeRxnBounds(test,rxns{i},1000,'u');
    
    %Simulate for growth
    sol = optimizeCbModel(test,[],'one');
    %Use 10% of growth as a barometer
    if sol.f >= 0.1*solution.f
        % It restored growth; store it
        fixRxns =[rxns{i};fixRxns];
    end
end


toc
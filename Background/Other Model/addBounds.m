function model = addBounds(model)

%For all reactions that are reversible, make the bounds -1000 and 1000
%For all irreversible reactions, make the bounds 0 and 1000

%Loop through the reactions
for i = 1:length(model.rxns)
    
    %Make the UB 1000
    model.ub(i) = 1000;
    
    %Look at the model.rev value
    %If it's 1, then make the lower bound -1000
    %Otherwise, make it 0
    if model.rev(i) == 1
        %Then make it -1000
        model.lb(i)=-1000;
    else
        %Make it 0
        model.lb(i)=0;
    end
end

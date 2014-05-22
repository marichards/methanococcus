function [missingMets,directions] = diagnoseMissingMet(model,rxn)

%Takes in a reaction that is essential for growth, turns off that reaction,
%and creates a source/sink for each metabolite.  Turns one of those
%source/sinks off at a time to diagnose which of the metabolites in the
%reaction are essential

%Inputs
%model - a COBRA model structure
%rxn - an essential reaction

%Outputs
%missingMets - the list of metabolites in rxn that require this reaction as a
%source or a sink
%directions - either 'source' or 'sink' to indicate the function of the
%reaction for the corresponding metabolite

%Matthew Richards, 05/21/2014


%Step 1: Turn off the reaction
model = changeRxnBounds(model,rxn,0,'b');

%Step 2: Find all the metabolites in the reaction
mets = findMetsFromRxns(model,rxn);

%Step 3: For each metabolite, create a sink/source reaction 
for i=1:length(mets)
    %Create an exchange reaction and call it a "sink"
    model = addReaction(model,sprintf('%s sink',mets{i}),...
        sprintf('%s <=> ',mets{i}));
end

%Step 4: Simulate growth, return error if it doesn't grow now
solution = optimizeCbModel(model,[],'one');
%The new model SHOULD grow if the given one did
if solution.f == 0
    error('Please check to ensure that your model grows')
end

%Step 5: Use solution to determine source/sink for each metabolite
source_sink = cell(length(mets),1);
for i=1:length(mets)
    %Find the index of the reaction for that metabolite
    [~,idx] = intersect(model.rxns,sprintf('%s sink',mets{i}));
    %Use the index to pull out the flux for that value  
    %If flux is greater than 0, it's a sink
    if solution.x(idx)>0
        source_sink{i}='sink';
    %If flux is less than 0, it's a source
    elseif solution.x(idx)<0
        source_sink{i}='source';
    %If flux is 0, I think there's an error, but I'll print 'zero'
    else
        source_sink{i}='zero';
    end
    
end

%Step 6: Knockout one source/sink at a time. If KO stops growth, record it
%in missingMets
missingMets={};
directions={};
for i=1:length(mets)
    %KO the sink for that metabolite
    ko_model = changeRxnBounds(model,sprintf('%s sink',mets{i}),0,'b');
    %Simulate growth
    ko_solution = optimizeCbModel(ko_model,[],'one');
    %If growth is 0, save it to missingMets
    if ko_solution.f==0
        missingMets=[missingMets; mets{i}];
        %Also add the source/sink designation to 'directions'
        directions=[directions; source_sink{i}];
    end
end

end

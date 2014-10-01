function [unbalancedRxns,unbalancedCharges]=checkChargeBalance(model)

%Take in a COBRA model, check charge balance for every reaction, return
%reactions that are imbalanced and a vector of their overall charges

%Inputs
%model: a COBRA model structure

%Outputs
%unbalancedRxns: a list of unbalanced reactions in 'model'
%unbalancedCharges: overall charges corresponding to unbalancedRxns

%Initiate the arrays for reactions and charges
unbalancedRxns = {};
unbalancedCharges = [];

%Loop through the reactions in the model
for i=1:length(model.rxns)
    
    %Find indices
    indices = find(model.S(:,i));
    %Pull out the coefficient of each metabolite
    coeffs = model.S(indices,i);
    
    %Pull out the charge of each metabolite and convert to a double
    charges = double(model.metCharge(indices));

    %Multiply the array of coefficients by the array of charges to get
    %overall charge
    %Each is nx1; multiply coeffs' by charges
    overall_charge = coeffs'*charges;
    
    %Check charge
    if overall_charge~=0
        %If it's not zero, add it to the lists
        unbalancedRxns=[unbalancedRxns;model.rxns{i}];
        unbalancedCharges=[unbalancedCharges;overall_charge];
    end
end

function [predicted_gr,full_growth_rates] = runATPMLOOCV(model,printFlag)%,growth_rates,ch4_rates)

% Check for a print flag
if nargin<2
    printFlag = true;
end

% Input measured values for growth rates and ch4 rates
full_growth_rates = [0.0814,0.0902,0.0892,0.0465,0.0705,0.0458,0.0602];

full_ch4_rates = [51.73,48.34,44.11,28.40,41.13,28.12,34.87];

% Create an array to store predicted growth rates
predicted_gr = zeros(1,length(full_growth_rates));

% Print header if not told no
if printFlag
    fprintf('\nActual Growth\tPredicted Growth\n')
end
% Loop through the dataset
for i=1:length(full_growth_rates)
    
    % Trim the ith element off both arrays
    growth_rates = setdiff(full_growth_rates,full_growth_rates(i));
    ch4_rates = setdiff(full_ch4_rates,full_ch4_rates(i));
    
    % Fit the model ATPM values
    [gam,ngam] = determineATPM(model,growth_rates,ch4_rates);
    fitted_model = changeATPM(model,gam,ngam);
    
    % Set the fitted model's secretion rate of methane
    fitted_model = changeRxnBounds(...
        fitted_model,'EX_cpd01024[e0]',full_ch4_rates(i),'b');
    
    % Simulate growth by optimizing biomass
    solution = optimizeCbModel(fitted_model);
    
    % Add the rate to predicted_gr
    predicted_gr(i) = solution.f; 
    % Print if not told no
    if printFlag
        fprintf('%0.4f\t\t\t%0.4f\n',full_growth_rates(i),solution.f)
    end
    
end





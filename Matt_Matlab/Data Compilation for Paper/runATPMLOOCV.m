function [predicted_yields,measured_yields] = runATPMLOOCV(model,printFlag)%,growth_rates,ch4_rates)

% Check for a print flag
if nargin<2
    printFlag = true;
end

% Input measured values for growth rates and ch4 rates
full_growth_rates = [0.0902,0.0892,0.0465,0.0705,0.0458,0.0602,0.0587,0.1297,0.1301];

full_ch4_rates = [48.34,44.11,28.40,41.13,28.12,34.87,35.17,64.97,61.58];

% Create an array to store predicted growth rates
predicted_gr = zeros(1,length(full_growth_rates));
pred_error = zeros(1,length(full_growth_rates));
% Print header if not told no
if printFlag
    fprintf('\nActual GR\tPredicted GR\t %% Error\n')
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
    fitted_model = setMethaneSecretion(fitted_model,full_ch4_rates(i));
    
    % Simulate growth by optimizing biomass
    solution = optimizeCbModel(fitted_model);
    
    % Add the rate to predicted_gr
    predicted_gr(i) = solution.f; 
    pred_error(i) = 100*(abs(solution.f-full_growth_rates(i))/full_growth_rates(i));
    % Print if not told no
    if printFlag
        fprintf('%0.4f\t\t\t%0.4f\t\t\t%0.1f\n',full_growth_rates(i),solution.f,pred_error(i))
    end
    
end

% Convert Growth rates into yields
measured_yields = full_growth_rates./full_ch4_rates*1000/log(2);
predicted_yields = predicted_gr./full_ch4_rates*1000/log(2);

% Add yield error I calculated
error_95 = [0.249,0.383,0.163,0.175,0.191,0.207,0.142,0.115,0.161];

figure(1)
%
%h = plot(full_ch4_rates,predicted_gr,'bo',...
% Plot error bar graph instead for the second point
h2 = errorbar(full_ch4_rates,measured_yields,error_95,'rx','MarkerSize',11);
set(h2,'LineWidth',2.5);
hold on
h1 = plot(full_ch4_rates,predicted_yields,'bo','MarkerSize',11);
set(h1,'LineWidth',2.5);
%data = [predicted_gr',full_growth_rates'];
%hb = bar(data,);
%set(hb(1),'FaceColor','b')
%set(hb(2),'FaceColor','y')
%set(hb,'XTickLabel',full_ch4_rates)
legend('Measured Growth Yield','Predicted Growth Yield','Location','northwest')
xlabel('Methane Evolution Rate ($$\frac{mmol}{gDCW \cdot h}$$)'...
    ,'Interpreter','latex','FontSize',14,'FontWeight','bold')
ylabel('Yield ($$\frac{gDCW}{mol Methane}$$)'...
    ,'Interpreter','latex','FontSize',14,'FontWeight','bold')
axis([25,70,0,4])
% Increase size of labels
set(gca,'FontSize',14);
% Take off the hold
hold off





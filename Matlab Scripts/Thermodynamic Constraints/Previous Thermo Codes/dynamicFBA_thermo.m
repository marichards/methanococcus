%Added dGVec to the outputs
function [concentrationMatrix,excRxnNames,timeVec,biomassVec,dG_overall] = ...
    dynamicFBA_thermo(model,substrateRxns,initConcentrations,initBiomass,timeStep,nSteps,plotRxns,exclUptakeRxns,water_rxn,pH,p_co2)
%dynamicFBA Perform dynamic FBA simulation using the static optimization
%approach
%
% [concentrationMatrix,excRxnNames,timeVec,biomassVec,dGVec]
% dynamicFBA_thermo(model,substrateRxns,initConcentrations,initBiomass,timeStep,nSteps,plotRxns,exclUptakeRxns)
%
%INPUTS
% model                 COBRA model structure
% substrateRxns         List of exchange reaction names for substrates
%                       initially in the media that may change (e.g. not
%                       h2o or co2)
% initConcentrations    Initial concentrations of substrates (in the same
%                       structure as substrateRxns)
% initBiomass           Initial biomass (must be non zero)
% timeStep              Time step size
% nSteps                Maximum number of time steps
%
%OPTIONAL INPUTS
% plotRxns              Reactions to be plotted (Default =
%                       {'EX_glc(e)','EX_ac(e)','EX_for(e)'})
% exclUptakeRxns        List of uptake reactions whose substrate
%                       concentrations do not change (Default =
%                       {'EX_co2(e)','EX_o2(e)','EX_h2o(e)','EX_h(e)'})
% 
%OUTPUTS
% concentrationMatrix   Matrix of extracellular metabolite concentrations
% excRxnNames           Names of exchange reactions for the EC metabolites
% timeVec               Vector of time points
% biomassVec            Vector of biomass values
%
% If no initial concentration is given for a substrate that has an open
% uptake in the model (i.e. model.lb < 0) the concentration is assumed to
% be high enough to not be limiting. If the uptake rate for a nutrient is
% calculated to exceed the maximum uptake rate for that nutrient specified
% in the model and the max uptake rate specified is > 0, the maximum uptake 
% rate specified in the model is used instead of the calculated uptake
% rate.
%
% NOTE: The dynamic FBA method implemented in this function is essentially 
% the same as the method described in
% [Varma, A., and B. O. Palsson. Appl. Environ. Microbiol. 60:3724 (1994)].
% This function does not implement the dynamic FBA using dynamic optimization approach
% described in [Mahadevan, R. et al. Biophys J, 83:3431-3440 (2003)].
%
% Markus Herrgard 8/22/06

% Modified: Matthew Richards 2/21/2013

%4/23/2013: Took out "originalBound > 0" on aboveOriginal, so now things
%that are irreversible forward can't change to uptake things instead

%Modified for methanogens:
%ch4, co2
if (nargin < 7)
    plotRxns = {'EX_ch4(e)','EX_co2(e)'};
end

% Uptake reactions whose substrate concentrations do not change
%%Need to modify this%%
%Took out O2, which isn't here
%Not sure what else to do with these...yeesh!
%Took out all of them for now...
if (nargin < 8)
    exclUptakeRxns = {};
end

%If water isn't given, won't be included in EQ
if (nargin < 10)
    water_rxn = {};
end

%If pH not given, assume neutral
if (nargin < 10)
    pH = 7;
end

%If pressure of co2 not given, assume atmospheric
if (nargin < 11)
    p_co2 = 0.00035;
end

%%%%%%%%%%%%%%%%%%%
%INITIAL CO2 EQ
%Fix CO2 and HCO3 concentrations using initial equilibrium
%*Be able to adjust for if others exist...?
%%%%%%%%%%%%%%%%%%%
%pH->c_h...adjust by factor of 1000
c_h=10^(3-pH);

%Henry's Law
c_co2 = p_co2*33.6;
%CO2->H2CO3
c_h2co3 = c_co2*1.7e-3;
%H2CO3->HCO3
c_hco3 = 0.25*c_h2co3/c_h;
%HCO3->CO3
c_co3 = 4.69e-8*c_hco3/c_h;

%ADD THEM TO THE SUBSTRATES/CONCENTRATIONS

%Methanogen-specific + pH (at least for now)
substrateRxns = [substrateRxns,'EX_h(e)','EX_co2(e)','EX_hco3(e)'];
initConcentrations = [initConcentrations,c_h,c_co2,c_hco3];

%%%%%%%%%%%%%%%%%%

% Find exchange rxns
excInd = findExcRxns(model,false);
excInd = excInd & ~ismember(model.rxns,exclUptakeRxns);
excRxnNames = model.rxns(excInd);
%find(excInd);
%length(excRxnNames);
% Figure out if substrate reactions are correct
%%%%%
%Comment this out (7 lines) for now, etoh problem
%missingInd = find(~ismember(substrateRxns,excRxnNames));
%if (~isempty(missingInd))
%    for i = 1:length(missingInd)
%        fprintf('%s\n',substrateRxns{missingInd(i)});
%    end
%    error('Invalid substrate uptake reaction!');
%end

%%%%%%%%%%%%%%%%%%%%
%MODIFY TO USE INTERSECT!
%%%%%%%%%%%%%%%%%%%%%%%
%%THEIR METHOD%%%%
% Initialize concentrations
%substrateMatchInd = ismember(excRxnNames,substrateRxns);
%concentrations = zeros(length(excRxnNames),1);
%concentrations(substrateMatchInd) = initConcentrations;
%%%%%%%%%%%%%%%%%%%%

%MY METHOD (04/03/2013)
%Initialize concentrations
%Should fix indexing issues
%sub_idx is index in supplied vector, exc_idx is index in model exchanges
[~,sub_idx,exc_idx] = intersect(substrateRxns,excRxnNames,'stable'); 
concentrations = zeros(length(excRxnNames),1);
%Concentrations will be indexed by exc_idx, must transfer over
concentrations(exc_idx)=initConcentrations(sub_idx);

%Pull out the indices for carbon species 
%Methanogen-specific
[~,co2_idx] = intersect(excRxnNames,'EX_co2(e)');
[~,hco3_idx] = intersect(excRxnNames,'EX_hco3(e)');
%pH
[~,h_idx] = intersect(excRxnNames,'EX_h(e)');

% Deal with reactions for which there are no initial concentrations
%%These are exchanges, so this should only be a problem if I don't specify
%%a starting concentration for something I need

originalBound = -model.lb(excInd);
noInitConcentration = (concentrations == 0 & originalBound > 0);
concentrations(noInitConcentration) = 1000;

biomass = initBiomass;

%THESE BREAK IT

%H2 Debug
%fprintf('Initial H2 Bounds:\n')
% Initialize bounds
uptakeBound =  concentrations/(biomass*timeStep);

%Debug with H2
%fprintf('Initial Concentration: %f\n',concentrations(34))
%fprintf('Initial Uptake Bound: %f\n',uptakeBound(34))
%Make sure bounds are not higher than what are specified in the model
aboveOriginal = (uptakeBound > originalBound); %& (originalBound > 0);

%Debug with H2
%fprintf('Is above original bound: %d\n',aboveOriginal(34))

uptakeBound(aboveOriginal) = originalBound(aboveOriginal);

%Debug with H2
%fprintf('Actual Uptake Bound: %f\n',uptakeBound(34))
model.lb(excInd) = -uptakeBound;
%%%%%%%%%%%%%%%%%%

concentrationMatrix = sparse(concentrations);
biomassVec = biomass;
timeVec(1) = 0;
%%%%Initiate the dG vector
dG_overall = [];

fprintf('Step number\tBiomass\n');
h = waitbar(0,'Dynamic FBA analysis in progress ...');

%DEBUG FOR THERMO
%substrateRxns
%concentrations(substrateMatchInd)

for stepNo = 1:nSteps
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %NEW PART
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    % Run Thermodynamic FBA
    %Pass only the substrate reaction concentrations
    [sol,gibbs] = optimizeThermoModel(model,substrateRxns,concentrations(exc_idx),biomass,310,water_rxn);
    
    %Debug the thermo inputs
    my_set = concentrations(exc_idx);
    for i = 1:length(substrateRxns)
        fprintf('%s\t%f\n',substrateRxns{i},my_set(i))
    end
        
    %printFluxVector(model,sol.x,'True','True')
    mu = sol.f;

    %Check if dG is valid
    if gibbs >=0
        fprintf('No feasible solution - Thermodynamics Violation\n');
        break;
    end
    
    %%This part needs to be commented out eventually
    %Check if nutrients exhausted
    if (sol.stat ~= 1 || mu == 0)
        fprintf('No feasible solution - nutrients exhausted\n');
        break;    
    end
    
    %Grab what's happening to exchange fluxes
    uptakeFlux = sol.x(excInd);
    
    %Grab the possible new concentration and biomasses
    new_biomass = biomass*exp(mu*timeStep);
    
    new_conc = concentrations - uptakeFlux/mu*biomass*(1-exp(mu*timeStep));
    %%%%%%%%%%%%%%%%%%%%
    %IF LOOP STARTS HERE
    if (condition)
    %If no solution, change the timestep to timeStep/10
    timeStep = timeStep/10;
    
    %Copy the regular uptake bound adjustment using the new timestep
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Update bounds for uptake reactions
    uptakeBound =  concentrations/(biomass*timeStep);
    
    % This is to avoid any numerical issues
    uptakeBound(uptakeBound > 1000) = 1000;
    
    % Figure out if the computed bounds were above the original bounds
    aboveOriginal = (uptakeBound > originalBound); %& (originalBound > 0);
    
    % Revert to original bounds if the rate was too high
    uptakeBound(aboveOriginal) = originalBound(aboveOriginal);
    uptakeBound(abs(uptakeBound) < 1e-9) = 0;

    model.lb(excInd) = -uptakeBound; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Otherwise, do what you regularly do
    else
    %Index the dG
    dG_overall(end+1) = gibbs;
    
    %DEBUG
    %model.rxns(excInd);
    biomass = new_biomass;
    %biomass = biomass*(1+mu*timeStep);
    biomassVec(end+1) = biomass;
    
    % Update concentrations
    concentrations = concentrations-takeaway;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Adjust carbon species concentrations
    %
    %Assign the new intial concentrations
    c_co2 = concentrations(co2_idx);
    c_hco3 = concentrations(hco3_idx);
    c_h = concentrations(h_idx);
    %Others are same as before (no output)...balance them with co2EQ code
    [c_co2,c_h2co3,c_hco3,c_co3,c_h,p_co2] = ...
        co2EQ(c_co2,c_h2co3,c_hco3,c_co3,c_h,p_co2);
    
    %Put in the 3 things we want
    concentrations(co2_idx) = c_co2;
    concentrations(hco3_idx) = c_hco3;
    concentrations(h_idx) = c_h;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %concentrations = concentrations + uptakeFlux*biomass*timeStep;
    concentrations(concentrations <= 0) = 0;
    concentrationMatrix(:,end+1) = sparse(concentrations);
            
    %%%%%%%%%%%%%%%%%
    %THESE BREAK IT
    % Update bounds for uptake reactions
    uptakeBound =  concentrations/(biomass*timeStep);
    
    %Debug with H2
    %%fprintf('New H2 Bounds:\n')
    %%fprintf('New Concentration: %f\n',concentrations(34))
    %%fprintf('New Uptake Bound: %f\n',uptakeBound(34))
    
    % This is to avoid any numerical issues
    uptakeBound(uptakeBound > 1000) = 1000;
    
    % Figure out if the computed bounds were above the original bounds
    aboveOriginal = (uptakeBound > originalBound); %& (originalBound > 0);
    
    %Debug with H2
    %%fprintf('Is above original bound: %d\n',aboveOriginal(34))
    %%fprintf('Original Bound: %f\n',originalBound(34))
    % Revert to original bounds if the rate was too high
    uptakeBound(aboveOriginal) = originalBound(aboveOriginal);
    uptakeBound(abs(uptakeBound) < 1e-9) = 0;
    
    %Debug with H2
    %%fprintf('Actual Uptake Bound: %f\n',uptakeBound(34))
    model.lb(excInd) = -uptakeBound;  
    %%%%%%%%%%%%%%%%%
    
    fprintf('%d\t%f\n',stepNo,biomass);
    waitbar(stepNo/nSteps,h); %Makes the progress bar go
    timeVec(stepNo+1) = stepNo*timeStep;
    end
end

%%%
%This ends the simulation part
%%%

if ( regexp( version, 'R20') )
        close(h);
end

selNonZero = any(concentrationMatrix>0,2);
concentrationMatrix = concentrationMatrix(selNonZero,:);
excRxnNames = excRxnNames(selNonZero);
selPlot = ismember(excRxnNames,plotRxns);

% Plot concentrations as a function of time
clf
figure(1)
subplot(1,2,1);
plot(timeVec,biomassVec,'LineWidth',1.5);
axis tight
xlabel('Time (h)','FontSize',10)
ylabel('Dry Cell Mass (g)','FontSize',10)
title('Biomass','FontSize',10);
set(gca,'LineWidth',1.5,'FontSize',10)
subplot(1,2,2);
plot(timeVec,concentrationMatrix(selPlot,:),'LineWidth',1.5);
axis tight
legend(strrep(excRxnNames(selPlot),'EX_',''),'FontSize',10);
xlabel('Time (h)','FontSize',10)
ylabel('Concentration (mM)','FontSize',10)
title('Concentration Profile','FontSize',10)
set(gca,'LineWidth',1.5,'FontSize',10)

%Plot the dG separately
figure(2)
plot(timeVec(1:end-1),dG_overall,'LineWidth',1.5);
xlabel('Time (h)','FontSize',10)
ylabel('\DeltaG Rate (kJ\cdotGDW^{-1}h^{-1})','FontSize',10)
%Previous implementation
%plot(timeVec(1:end-1),dG_overall);
set(gca,'LineWidth',1.5,'FontSize',10)

%Plot the kJ/hr
%You have to .* them
kj_hr = dG_overall.*biomassVec(1:end-1);
figure(3)
plot(timeVec(1:end-1),kj_hr,'LineWidth',1.5);
xlabel('Time (h)','FontSize',10)
ylabel('\DeltaG Rate (kJ/h)','FontSize',10)
set(gca,'LineWidth',1.5,'FontSize',10)

%ADDITIONAL THERMO MEASURES...
%Integrate to find total kJ/GDW
total_dG_gdw = trapz(timeVec(1:end-1),dG_overall);
fprintf('\nTotal Specific Free Energy: %f kJ/GDW\n',total_dG_gdw)

total_dG = trapz(timeVec(1:end-1),kj_hr);
fprintf('\nTotal Free Energy: %f kJ\n',total_dG)

end

function [c_co2,c_h2co3,c_hco3,c_co3,c_h,p_co2] = co2EQ(c0_co2,c0_h2co3,c0_hco3,c0_co3,c0_h,p0_co2)

%Calculate the concentrations of carbon species given initial conditions
%Make the guesses equal to the initial values!
guess = [0,0,0,0,c0_co2,c0_h2co3,c0_hco3,c0_co3,c0_h,p0_co2];

%Solve using system of equations
results=fsolve(@eqns,guess,[],c0_co2,c0_h2co3,c0_hco3,c0_co3,c0_h,p0_co2);

%Pull the things I want out of the results...
c_co2 = results(5);
c_h2co3 = results(6);
c_hco3 = results(7);
c_co3 = results(8);
c_h = results(9);
p_co2 = results(10);

end

%Supply set of functions, matlab solves it...define functions:
function fcns = eqns(z,c0_co2,c0_h2co3,c0_hco3,c0_co3,c0_h,p0_co2)

%**z contains values for all unknowns
%Length of equations/unkonwns
x0 = z(1);
x1 = z(2);
x2 = z(3);
x3 = z(4);
c_co2 = z(5);
c_h2co3 = z(6);
c_hco3 = z(7);
c_co3 = z(8);
c_h = z(9);
p_co2 = z(10);

%Equilibrium data (Matt B article)
%H = 2; %mM/atm, Henry's Constant
%k1 = 1.3e-3; %unitless
%k2 = 0.2; %mM
%k3 = 4.45e-4; %mM

%Equilibrium data (Wikipedia)
H = 33.6; %mM/atm, Henry's Constant
k1 = 1.7e-3; %unitless
k2 = 0.25; %mM
k3 = 4.45e-8; %mM

%Equilibrium relationships
fcns(1) = c_co2 - p_co2*H;
fcns(2) = c_h2co3 - k1*c_co2;
fcns(3) = c_hco3 .* c_h - k2*c_h2co3;
fcns(4) = c_co3 .* c_h - k3*c_hco3;

%Mass Balances
%In solution
fcns(5) = c0_h + x1 + x3 - c_h;
fcns(6) = c0_co2 + x0 - x1 - c_co2;
fcns(7) = c0_h2co3 + x1 - x2 - c_h2co3;
fcns(8) = c0_hco3 + x2 - x3 - c_hco3;
fcns(9) = c0_co3 + x3 - c_co3;

%Gaseous
fcns(10) = p0_co2 - x0/H - p_co2;
end

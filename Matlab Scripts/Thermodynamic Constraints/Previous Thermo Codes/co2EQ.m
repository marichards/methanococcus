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


function [x,y] = testSolve(x0,y0)

guess=[x0,y0];

results = fsolve(@eqns,guess,[],x0,y0);

x=results(1);
y=results(2);

end


%Note: x is 4, y is 2
function fcns = eqns(z,x0,y0)

x = z(1);
y = z(2);

fcns(1) = x.^2-y-x0;
fcns(2) = x.*y-y0;

end
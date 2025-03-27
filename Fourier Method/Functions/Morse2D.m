function [v2_12] = Morse2D(epsilon, alpha , r_e, B , n)
% This function is used to evaluate the minimum point for cut off process
% then it find the potential according to the distance given with respect
% to the smoothing process in the part which goes to Inf.
% 
% epsilon = 0.3492; % is 0.3429 in thesis (D). (eV Units) 
% alpha = 1.357; % (Units of 1/Angstrom)
% r_e = 2.866; % The equilibrium separation for copper. (Angstrom Units)
% B = distance 12 
% epsilon = constant value of the epsilon in Morse potential
% alpha =  constant value of the alpha in Morse potential
% r_e = equilibrium bond distance
% U3 = @(epsilon,alpha,r,r_e)(epsilon .* (exp(-2 * alpha .*(r - r_e))  - (2 .* exp(-alpha .*(r - r_e)))));
% evaluation of r_mid
% syms y
% eqn2 =  (exp( -2 * alpha .*(y - r_e))  - 2 .* exp( -alpha .*( y- r_e))) == 25; 
% solr2 = solve(eqn2,y);
% solr2 = double(solr2);
% r_mid_morse = solr2(1);

A = 50;
r_mid_morse = r_e - log(1+sqrt(1+A))/alpha;
% smooth part constants ( V_morse = a * r^2 + b * r + c )
tmp = exp(-alpha * (r_mid_morse - r_e));
a = (epsilon * alpha ) ./ r_mid_morse * tmp * (1. - tmp);
b = 0;
c = epsilon * tmp * (tmp - 2.) - a * r_mid_morse^2. - b * r_mid_morse;

v2_12 = epsilon * (exp(-2. * alpha * (B - r_e)) - 2. * exp(-alpha * (B - r_e)));
filter = (B < r_mid_morse);
v2_12(filter) = a * B(filter).^2. + b * B(filter) + c;
end
% solr(abs(imag(solr))<eps & real(solr)>0)
% not smoothed version Morse:
% %  v2_12 = epsilon * (( 1 - exp( -alpha * ( distance12 - 1))).^2 - 1);
% B = distance12;
% v2_12 =  epsilon .*(exp( -2 * alpha .*(B - r_e))  - 2.* exp( -alpha .*( B - r_e)));

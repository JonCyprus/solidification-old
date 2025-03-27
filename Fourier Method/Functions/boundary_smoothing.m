function [damping_factor] = boundary_smoothing(r)
% Function is intended to decrease interactions farther from the origin to
% prevent boundary conditions from interfering.
% The function is a 5th degree polynomial that smoothly decreases the
% interactions from the origin to 0 starting at distance ro to r1 from the
% origin
% Polynomial's General Form: g(r) = a + br + cr^2 + dr^3 + er^4 + fr^5
% Constants were calculated using r0 = 4.4r_e and r1 = 4.9r_e
% r_e is the equilibrium bond distance, grid in twobodydistribution.m is
% 5r_e in each direction.
% r_e = 2.866

r_e = 2.866;
r_inner = 4.4 * r_e;
r_outer = 4.9 * r_e;

damping_factor = zeros(size(r));
damping_factor(r < r_inner) = 1.0;

a = 413409.1743022878;
b = -155701.4152797109;
c = 23434.223422034785;
d = -1761.8134897417644;
e = 66.16364790843001;
f = -0.9929338091898702;

filter = (r >= r_inner & r < r_outer);
damping_factor(filter) = a + b*r(filter) + c*r(filter).^2 + d*r(filter).^3+ e*r(filter).^4 + f*r(filter).^5;

end


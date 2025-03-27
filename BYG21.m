% This script finds the two-body distribution function for a uniform liquid
% in two dimensions.

% Overall parameters
max_t = 1e-4;   % increment of time
kbT = 1.2;       % temperature, regulates diffusion term
G = 1.;         % overall mobility constant
red = 16;       % factor that reduces the set of frequencies

% Two dimensional parameters
L = 40;         % length of simulation cell edge
N = 1700;        % number of particles

% Uniform liquid parameters
n = 4096;       % number of lattice sites in one dimension
R = 20.;        % maximum distance of interacting particles

%%%%%%%%%%%%%%%%%%%%

% Quantities for the normalization of the distribution functions
one_body = N / L^2;
two_body = N * ( N - 1 ) / L^4;

% Quantities for the Fourier-Bessel series
% Most efficient when ( k_max * bond length ) is around forty
z = bessel_root( 1, n )';
k = z( 2:n/red ) / R;
r = R / z(n) * z;
rdr = integrate( r );

% Transformation matrices
c = sqrt( 2. ) ./ ( R * besselj( 0, z( 2:n/red ) ) );
T0 = besselj( 0, k * r' );
T1 = besselj( 1, k * r' );

% Potentials: exp-6, 10-6, Morse, hard sphere
% v = 1. / ( 1. - 6. / 12. ) * ( ( 6. / 12. ) * exp( 12. * ( 1 - r / 1. ) ) - ( 1. ./ r ).^6 );
v = 1. / ( 10. - 6. ) * ( 6. * ( 1. ./ r ).^10 - 10. * ( 1. ./ r ).^6 );
% v = ( 1. - exp( -2.4 * ( r - 1 ) ) ).^2 - 1.;
% v = 15. * ( tanh( 15. * ( 1. -  r ) ) + 1. );

% Modifies the potential - should change with kbT and density
% v = plateau_a( v, r, 2. );
% v = plateau_b( v, r, 2. );
v( 1:find( v > 25., 1, 'last' ) ) = 25.;
v = cutoff( v, r, 6., 6.5 );

[ dv, lv ] = taylor( v, r );

% Initial two-body distribution function
p = ones( n, 1 ) * two_body;
p = peak( r, rdr, p, 0., 0.75, 1., 0.75, 1.25 );
lr = lrange( r, 12., 12.5 );

% Reference for forward and reverse transform
% P0 = c .* ( T0 * ( p0 .* rdr ) );
% p0 = T0' * ( c .* P0 );

time = 0;
for a = 0:1000000
    [ dp, lp ] = taylor( p, r );
    P0 = c .* ( T0 * ( ( lr .* ( p - two_body) + two_body ) .* rdr ) );
    
    % Calculates change
    fst = kbT * lp;
    snd = dv .* dp + lv .* p;
    trd = p .* ( T0' * ( P0 .* c .* ( T0 * ( snd .* rdr ) ) ) );
    fth = dp .* ( T1' * ( P0 .* c .* ( T1 * ( dv .* p .* rdr ) ) ) );
    change = 2. * G * ( fst + snd + 2. * pi * ( trd + fth ) / one_body^3 );
    
    % Enforces positivity
    tmp = -p ./ change;
    dt = min( [ min( tmp( tmp > 0 ) ) / 2., max_t ] );
    time = time + dt;
    p = p + dt * change;
    
    lim = n;
    % Visualization
    if mod( a, 1000 ) == 0
        disp(a);
        disp(dt);
        figure(1);
        plot( r(1:lim), p(1:lim) / two_body );
        figure(2), clf, hold on;
        plot( r(1:lim), fst(1:lim), 'r' );
        plot( r(1:lim), snd(1:lim), 'g' );
        plot( r(1:lim), trd(1:lim), 'b' );
        plot( r(1:lim), fth(1:lim), 'c' );
        figure(3);
        plot( r(1:lim), change(1:lim) / two_body, 'k' );
    end
    
    % Stopping criterion
    if max( abs( change ) ) / max( p ) < 0.1;
        break;
    end
end


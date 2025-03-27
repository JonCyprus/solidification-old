% Evolution Equation for the two-body distribution function
% Morse potential constants
epsilon = 0.3492; % epsilon = constant value of the epsilon in Morse potential
alpha = 1.357; % alpha =  constant value of the alpha in Morse potential
r_e = 2.866; % r_e = equilibrium bond distance

% Constants
G = 1; % Gamma
n = 2000; % number of grid points along the edge
L = 10 * r_e; %length of simulation cell edge
N = 80; % number of particles
% Expect in MD to have density of 0.14 atoms/Angstrom^2 @ 1350K
% This corresponds to 114 atoms when L = 10 * r_e
T = 1200; %Temperature of the system in Kelvin [K]
step =  1e-6; % length of time step
total_step = 2000000; % total number of steps
plotting_step =200; % every plotting step it plots p0_12

delta = L / n; %the spacing of the grid
dA = delta * delta;
cell_area = L^2;
uniform_density = N / cell_area; 
p_ref = N * (N-1) / cell_area^2;

Kb = 8.617333262e-5; %Boltzmann constant Kb = 8.617333262...e-5 [eV/K]
KbT = Kb * T;
Beta = 1 / KbT;

% save('old.mat');
% save('new.mat');

% Definition of Particles
% First Particle
x1 = 0;
y1 = 0;

% Second Particle
x2 = linspace(-L/2,L/2,n+1)'* ones(1,n+1);
x2 = x2(1:end-1, 1:end-1);
y2 = ones(1,n+1)' * linspace (-L/2,L/2,n+1);
y2 = y2(1:end-1, 1:end-1);

% Evaluation of distance
distance12 = sqrt(x2.^2 + y2.^2);
smoothing = boundary_smoothing(distance12);

% Definition of Potential:
% Morse Potential
v2_12 = Morse2D(epsilon, alpha, r_e, distance12, n); % Changed to Morse2D

%To Speed up Fast Fourier Transforms
fftw('dwisdom', []);
fftw('planner', 'patient');
ifft2(fft2(v2_12), 'symmetric');
fftinfo = fftw('dwisdom');
fftw('dwisdom', fftinfo);

% Third Particle
x3 = x2;
y3 = y2;

% Evaluation of distance & potential v_2
% distance13 = distance12;
 
% SECOND PARTICLE 1---2
lap_v2_12 = laplacian3(v2_12, delta);
grad_v2_12_x = nablax(v2_12,delta);
grad_v2_12_y = nablay(v2_12,delta);

% THIRD PARTICLE 1---3
lap_v2_13 = lap_v2_12;
grad_v2_13_x = grad_v2_12_x;
grad_v2_13_y = grad_v2_12_y;  

% SECOND and THIRD 2---3
lap_v2_23 = circshift(lap_v2_12, [n/2 , n/2]);
grad_v2_23_x = circshift (grad_v2_12_x, [n/2 , n/2]);
grad_v2_23_y = circshift(grad_v2_12_y, [n/2 , n/2]);

% values relates to P_0 that do not change as it iterates:j
% SECOND PARTICLE 1---2
% p0_12 = exp(-Beta.*v2_12);
% p0_12 = p0_12 * (N * (N-1)/cell_area) / (sum(p0_12(:)) * dA);
% 
% tmp = 0.;
% fig_num = 1;
% abs_delta = [];
% min_step = 1000;

start = 780000;
for s = start:total_step
    grad_p0_12_x = nablax(p0_12, delta);
    grad_p0_12_y = nablay(p0_12, delta);
    %THIRD PARTICLE 1---3
    % Second and Third 2--3
    p0_23 = circshift(p0_12, [n/2, n/2]);
    % Evaluation of formula:
    %%% 1st term
    first = KbT .* laplacian3(p0_12, delta); % lap_p0_12 = laplacian3(p0_12, delta);
    %%% 2nd term
    second = lap_v2_12 .* p0_12 + ...
        grad_p0_12_x .* grad_v2_12_x + grad_p0_12_y .* grad_v2_12_y;
    %%% 3rd term
    % 1st Integral
    wx1 = grad_v2_13_x .* grad_p0_12_x; % grad_p0_13_x = grad_p0_12_x;
    wy1 = grad_v2_13_y .* grad_p0_12_y; % grad_p0_13_y = grad_p0_12_y;
    fw1 = fft2(lap_v2_13 .* p0_12 + wx1 + wy1); % p0_13 = p0_12;
    fp1 = fft2(p0_23);
    int1 = ifft2(fw1 .* fp1, 'symmetric');
    % 2nd Integral
    wx2 = grad_v2_23_x .* circshift(grad_p0_12_x, [n/2, n/2]); % grad_p0_23_x = circshift(grad_p0_12_x, [n/2, n/2]);
    wy2 = grad_v2_23_y .* circshift(grad_p0_12_y, [n/2, n/2]); % grad_p0_23_y = circshift(grad_p0_12_y, [n/2, n/2]);
    fw2 = fft2(lap_v2_23 .* p0_23 + wx2 + wy2);
    fp2 = fft2(p0_12); % p0_13 = p0_12;
    int2 = ifft2 (fw2 .* fp2, 'symmetric');
    % 3rd Integral
    fwx3 = fft2(p0_12 .* grad_v2_13_x); % p0_13 = p0_12;
    fwy3 = fft2(p0_12 .* grad_v2_13_y); % p0_13 = p0_12;
    intx3 = ifft2(fwx3 .* fp1, 'symmetric'); % fp3 = fp1;
    inty3 = ifft2(fwy3 .* fp1, 'symmetric'); % fp3 = fp1;
    % 4th Integral
    fwx4 = fft2(p0_23 .* grad_v2_23_x);
    fwy4 = fft2(p0_23 .* grad_v2_23_y);
    intx4 = ifft2(fwx4 .* fp2, 'symmetric'); % fp4 = fp2;
    inty4 = ifft2(fwy4 .* fp2, 'symmetric'); % fp4 = fp2;
    % Total Contributions:
    integral1 = int1 .* p0_12 * dA;
    integral2 = int2 .* p0_12 * dA;
    integral3 = (intx3 .* grad_p0_12_x + inty3 .* grad_p0_12_y) * dA;
    integral4 = (intx4 .* grad_p0_12_x + inty4 .* grad_p0_12_y) * dA;
    third = (integral1 + integral2 + integral3 + integral4) ./ (uniform_density^3);
    % Final Step
    change = 2 * G * (first + second + third); % dP_12
    p0_23 = p0_12; % Previous value
    p0_12 = p0_12 + change * step;

    % Smoothing to 0 from r0 to r1
    p0_12 = p0_12 + abs(min(min(p0_12(:)), 0.0));
    p0_12 = p0_12 .* smoothing + p_ref * (1. - smoothing);
    p0_12 = p0_12 * (N * (N-1) / cell_area) / (sum(p0_12(:)) * dA);
    
    if mod(s,50000) ==0
        save(num2str(s));
    end
    % if mod(s, 100) == 0
    %     delete('old.mat');
    %     movefile('new.mat', 'old.mat');
    %     disp(['old.mat corresponds to step ', num2str(s - 100)]);
    %     save('new.mat');
    %     disp(['new.mat corresponds to step ', num2str(s)]);
    % end
    if any(abs(p0_12(:)) > 1.e3)
        disp(['Limit exceeded on step ', num2str(s)]);
        save('final.mat');
        break;
    end
   
    tmp = sum(sum(abs(p0_12 - p0_23))) * dA / step;
    disp([num2str(s), ', ', num2str(tmp), ', ', num2str(min(p0_12(:)))]);

    abs_delta(end + 1) = tmp;
    if numel(abs_delta) > 2
        tmp = mean(abs_delta((end-2):end));
    end
    if s > min_step && tmp < 1.
        if step < 5e-7
            disp(['Cut here: ', num2str(s)])
            break;
        else
            step = step / 10.;
            disp(['Reduce step size: ', num2str(step)])
            min_step = s + 1000;
        end
    end
    if mod(s, plotting_step) == 0
        figure(fig_num);
        hh = surf(y2, x2, p0_12, 'edgealpha',0);
        colorbar
        az = 90;
        el = 90;
        view([az,el])
        axis tight
        title(['Total steps = ', num2str(fig_num * s), '   n =',num2str(n),'   K_bT =', num2str(KbT),'   N=', num2str(N)]);
        saveas(hh, sprintf('figure_Zdirect_%d.jpg',fig_num));
        az = 90;
        el = 0;
        view([az,el])
        axis tight
        title (['Total steps = ',num2str(fig_num * s),'   n =',num2str(n),'   K_bT =',num2str(KbT),'   N =', num2str(N)]); 
        saveas(hh, sprintf('figure_Xdirect_%d.jpg',fig_num));
        if rem(fig_num,3) == 0
            close all
        end
    end
end
% Plot:
close all
save('p12_0.1163_450k.mat', 'p0_12', 'v2_12', 'dA', 'grad_p0_12_x', 'grad_p0_12_y', 'grad_v2_12_x', 'grad_v2_12_y', 'lap_v2_12', 'y2', 'x2')
% Plot absolute value part
fig_num = fig_num + 1;
figure(fig_num)
ps = size(abs_delta);
plot_abs = linspace(1,ps(2),ps(2))'; 
h2 = plot(plot_abs, abs_delta);
title(['maximum absolute','   ', 'Total steps = ', num2str(fig_num * plotting_step),'   ','KbT =',num2str(KbT),'   ','N =', num2str(N)]); 
saveas(h2, 'maximum_absolute.jpg')

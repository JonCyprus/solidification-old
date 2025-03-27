% Evolution Equation for the two-body distribution function

% Add functions
addpath('Functions');

% Morse potential constants
D = 0.3492; % constant value of the epsilon in Morse potential
alpha = 1.357; % constant value of the alpha in Morse potential
Re = 2.866; % equilibrium bond distance
Rc = 15;  

% Constants
G = 1; % Gamma
n = 2000; % number of grid points along the edge
N = 400; % number of particles
uniform_density = 0.1362; % uni-form density of a liquid from MD
L = sqrt(N / uniform_density); %length of simulation cell edge done in terms of density, N, and r_e to preserve density
% Expect in MD to have density of 0.1362 atoms/Angstrom^2 @ 1350K
T = 1200; %Temperature of the system in Kelvin [K]
step =  1e-4; % length of time step
total_step = 2000000; % total number of steps
plotting_step = 500; % every plotting step it plots p0_12
save_step = 100000;

delta = L / n; %the spacing of the grid
dA = delta * delta;
cell_area = L^2;
%uniform_density = N / cell_area; put in constants to keep the density same
p_ref = N * (N-1) / cell_area^2;

Kb = 8.617333262e-5; %Boltzmann constant Kb = 8.617...e-5 [eV/K]
kbT = Kb * T;
Beta = 1 / kbT;

% save('old.mat');
% save('new.mat');

% New Gamma Terms from non-constant mobility
r = linspace(0, L, n);
g_xy = gamma_g(G, r, 4 * Re, 7 * Re);
g_xy = interp_data(L, n, Re .* 8, 2, r, g_xy, 0.);

% Morse Potential
f1 = poly_solver([[0 8 0], [0 0 1], [0 -2 2]], 'r'); % polynomial at origin
v2_12 = morse_modified(r, f1, 0.1 * Re, 0.75 * Re);
[ dv, lv ] = taylor( v2_12, r );
v2_12 = interp_data(L, n, Rc, 2, r, v2_12, 0.);

% Definition of Particles
% First Particle
% x1 = 0;
% y1 = 0;

% Second Particle
x2 = linspace(-L/2, L/2, n);
[x2, y2] = meshgrid(x2, x2);

% Old incorrect way
% x2 = linspace(-L/2,L/2,n+1)'* ones(1,n+1);
% x2 = x2(1:end-1, 1:end-1);
% y2 = ones(1,n+1)' * linspace (-L/2,L/2,n+1);
% y2 = y2(1:end-1, 1:end-1);

% Evaluation of distance
distance12 = sqrt(x2.^2 + y2.^2);
smoothing = g_xy; %boundary_smoothing(distance12);

% To Speed up Fast Fourier Transforms
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
 
% PARTICLE 1--2 POTENTIAL GRAD/LAPLACIAN
% New method to have a smoother laplacian by interpolating from 1-D
lap_v2_12 = interp_data(L, n, Rc + 1, 2, r, lv, 0.); 
dv_dr = interp_data(L, n, Rc+1, 2, r, dv, 0.);
r_2d = sqrt(x2.^2 + y2.^2);
grad_v2_12_x = dv_dr .* x2 ./ r_2d;
grad_v2_12_y = dv_dr .* y2 ./ r_2d;

% Handling singularity division by 0 
grad_v2_12_x(r_2d == 0) = 0;
grad_v2_12_y(r_2d == 0) = 0;

% Old finite difference methods
%lap_v2_12 = laplacian3(v2_12, delta);
%grad_v2_12_x = nablax(v2_12,delta);
%grad_v2_12_y = nablay(v2_12,delta);

% THIRD PARTICLE 1---3
lap_v2_13 = lap_v2_12;
grad_v2_13_x = grad_v2_12_x;
grad_v2_13_y = grad_v2_12_y;  

% SECOND and THIRD 2---3
lap_v2_23 = circshift(lap_v2_12, [n/2 , n/2]);
grad_v2_23_x = circshift (grad_v2_12_x, [n/2 , n/2]);
grad_v2_23_y = circshift(grad_v2_12_y, [n/2 , n/2]);

% GRAD OF g(x,y) FROM NEW TREATMENT OF GAMMA: 
grad_g_x = nablax(g_xy, delta);
grad_g_y = nablay(g_xy, delta);

% SECOND PARTICLE 1---2 Initial conditions
p0_12 = exp(-Beta .* v2_12);
p0_12 = p0_12 * (N * (N-1)/cell_area) / (sum(p0_12(:)) * dA);
p0_23 = circshift(p0_12, [n/2, n/2]); % for gpuArray

% Initialize loop parameters
tmp = 0.;
fig_num = 1;
abs_delta = [];
min_step = 1000;

% Before loop, move data to the GPU for Parallel Computing Toolbox
p0_12 = gpuArray(p0_12);
p0_23 = gpuArray(p0_23);
v2_12 = gpuArray(v2_12);
smoothing = gpuArray(smoothing);
g_xy = gpuArray(g_xy);
grad_g_x = gpuArray(grad_g_x);
grad_g_y = gpuArray(grad_g_y);
lap_v2_12 = gpuArray(lap_v2_12);
grad_v2_12_x = gpuArray(grad_v2_12_x);
grad_v2_12_y = gpuArray(grad_v2_12_y);
lap_v2_13 = gpuArray(lap_v2_13);
grad_v2_13_x = gpuArray(grad_v2_13_x);
grad_v2_13_y = gpuArray(grad_v2_13_y);
lap_v2_23 = gpuArray(lap_v2_23);
grad_v2_23_x = gpuArray(grad_v2_23_x);
grad_v2_23_y = gpuArray(grad_v2_23_y);

start = 0;
for s = start:total_step

    % Recalculate grads for p0_12
    grad_p0_12_x = nablax(p0_12, delta);
    grad_p0_12_y = nablay(p0_12, delta);

    % Second and Third 2--3
    p0_23 = circshift(p0_12, [n/2, n/2]);
    grad_p0_23_x = circshift(grad_p0_12_x, [n/2, n/2]);
    grad_p0_23_y = circshift(grad_p0_12_y, [n/2, n/2]);

    % w and z terms
    w_12 = lap_v2_12 .* p0_12 + (grad_v2_12_x .* grad_p0_12_x) + ...
                                (grad_v2_12_y .* grad_p0_12_y);

    w_23 = lap_v2_23 .* p0_23 + (grad_v2_23_x .* grad_p0_23_x) + ...
                                (grad_v2_23_y .* grad_p0_23_y);

    z_23_x = - grad_v2_23_x .* p0_23;
    z_23_y = - grad_v2_23_y .* p0_23;

    % Evaluation of formula:
    %%% 1st term
    first = 2 * G .* g_xy .* (kbT * laplacian3(p0_12, delta) + w_12);

    %%% 2nd term
    second_x = (kbT * grad_p0_12_x + grad_v2_12_x .* p0_12);
    second_y = (kbT * grad_p0_12_y + grad_v2_12_y .* p0_12);
    second = 2 * G * (grad_g_x .* second_x + grad_g_y .* second_y);

    %%% 3rd term
    third_factor = 2 * G .* g_xy .* (p0_12 / (uniform_density ^ 3));
    fw_23 = fft2(w_23);
    fp0_13 = fft2(p0_12); % is p0_13 in derivation
    int1 = ifft2(fw_23 .* fp0_13, 'symmetric') .* dA;
    int1 = (int1 + circshift(int1, [1, 0])) / 2.;
    int1 = (int1 + circshift(int1, [0, 1])) / 2.;
    % int1 = real(ifft2(fw_23 .* fp0_13)) .* dA;
    third = third_factor .* int1;

    %%% 4th term
    fourth_x = (2 * g_xy .* grad_p0_12_x + 2 * grad_g_x .* p0_12);
    fourth_y = (2 * g_xy .* grad_p0_12_y + 2 * grad_g_y .* p0_12);
    int2_x = ifft2(fft2(z_23_x) .* fp0_13, 'symmetric') .* dA;
    int2_x = (int2_x + circshift(int2_x, [1, 0])) / 2.;
    int2_x = (int2_x + circshift(int2_x, [0, 1])) / 2.;
    int2_y = ifft2(fft2(z_23_y) .* fp0_13, 'symmetric') .* dA;
    int2_y = (int2_y + circshift(int2_y, [1, 0])) / 2.;
    int2_y = (int2_y + circshift(int2_y, [0, 1])) / 2.;
    % int2_x = real(ifft2(fft2(z_23_x) .* fp0_13)) .* dA;
    % int2_y = real(ifft2(fft2(z_23_y) .* fp0_13)) .* dA;
    fourth = (fourth_x .* int2_x) + (fourth_y .* int2_y);
    fourth = (G / (uniform_density ^3)) * fourth;
    
    % Final Step
    change = first + second + third + fourth; % dP_12
    p0_23 = p0_12; % store for use in tmp
    p0_12 = p0_12 + change * step;

    % Smoothing to 0 from r0 to r1
    p0_12 = p0_12 + abs(min(min(p0_12(:)), 0.0));
    p0_12 = p0_12 .* smoothing + p_ref * (1. - smoothing);
    p0_12 = p0_12 * (N * (N-1) / cell_area) / (sum(p0_12(:)) * dA);
    
    % Saving criterion
    if mod(s,save_step) == 0
        if s > 0
            save(num2str(s));
        end
    end
    % if mod(s, 100) == 0
    %     delete('old.mat');
    %     movefile('new.mat', 'old.mat');
    %     disp(['old.mat corresponds to step ', num2str(s - 100)]);
    %     save('new.mat');
    %     disp(['new.mat corresponds to step ', num2str(s)]);
    % end

    % End criterion for p0_12 being too large
    if any(abs(p0_12(:)) > 1.e3)
        disp(['Limit exceeded on step ', num2str(s)]);
        save('final.mat');
        break;
    end
   
    % Tracking of total mass
    total_mass = sum(p0_12(:)) * dA;  

    tmp = sum(sum(abs(p0_12 - p0_23))) * dA / step;
    disp(['Step: ', num2str(s), ...
            ', tmp: ', num2str(tmp), ', min(p0_12): ', ... 
            num2str(min(p0_12(:))), ... 
            ', Total Mass: ', num2str(total_mass)]);
    
    abs_delta(end + 1) = tmp;
    if numel(abs_delta) > 2
        tmp = mean(abs_delta((end-2):end));
    end

    % Dynamically change timestep based on tmp
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
    

    % Visualization
    if mod(s, plotting_step) == 0
        figure(fig_num);
        hh = surf(x2, y2, p0_12, 'edgealpha',0);
        colorbar
        az = 90;
        el = 90;
        view([az,el])
        axis tight
        title(['Total steps = ', num2str(s), '   n =',num2str(n),'   K_bT =', num2str(kbT),'   N=', num2str(N)]);
        xlabel('x');
        ylabel('y');
        saveas(hh, sprintf('figure_Zdirect_%d.jpg',fig_num));
        filename1= fullfile('Figure 1', [ num2str(s/plotting_step), '.png']);
        saveas(gcf,filename1);
        az = 90;
        el = 0;
        view([az,el])
        axis tight
        title (['Total steps = ',num2str(s),'   n =',num2str(n),'   K_bT =',num2str(kbT),'   N =', num2str(N)]); 
        %saveas(hh, sprintf('figure_Xdirect_%d.jpg',fig_num));
        
        % Graph the first change term
        figure(10);
        hh = surf(x2, y2, first, 'edgealpha',0);
        colorbar
        az = 90;
        el = 90;
        view([az,el])
        axis tight
        title (['First Change Term: Step - ', num2str(s)])
        xlabel('x');
        ylabel('y');
        
        % Graph the second change term
        figure(11);
        hh = surf(x2, y2, second, 'edgealpha',0);
        colorbar
        az = 90;
        el = 90;
        view([az,el])
        axis tight
        title(['Second Change Term: Step - ', num2str(s)])
        xlabel('x');
        ylabel('y');
        
        % Graph the third change term
        figure(12);
        hh = surf(x2, y2, third, 'edgealpha',0);
        colorbar
        az = 90;
        el = 90;
        view([az,el])
        axis tight
        title(['Third Change Term: Step - ', num2str(s)])
        xlabel('x');
        ylabel('y');
        
        % Graph the fourth change term
        figure(13);
        hh = surf(x2, y2, fourth, 'edgealpha',0);
        colorbar
        az = 90;
        el = 90;
        view([az,el])
        axis tight
        title(['Fourth Change Term: Step - ', num2str(s)])
        xlabel('x');
        ylabel('y');

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
title(['maximum absolute','   ', 'Total steps = ', num2str(fig_num * plotting_step),'   ','KbT =',num2str(kbT),'   ','N =', num2str(N)]); 
saveas(h2, 'maximum_absolute.jpg')

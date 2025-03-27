% Evoluation equation for the two-body distirbution function 
%----------------------------------------------------------
% FOURIER CODE
%----------------------------------------------------------
% Constants
G = 1; % Gama
n = 2048;     % number of grid points along the edge
L = 13.52 * 2.866;     % length of simulation cell edge
A = L.^2;   % Area of the system
% Density = (L^2) / ((sqrt(3)/2));
N = 200;    % number of particles 
KbT = 0.125;      %1450 K  
% T = 1450;
beta = 1/KbT;
% The Boltzmann constant : k_b =  8.617332478e-5 [ eV/K ] 
timestep = 1e-6;   % length of time step
total_step = 100000; % total number of steps
ploting_step = 5000; % every ploting step it plots p0_12
% MORSE potential constants
epsilon = 0.3492;
alpha = 1.357;
r_e = 2.866;
constantP = N/(L^2); 
delta = L/n; % the spacing of the grid
% Strategy 1:
% Neg_parameter = 0.8; 
% gamma = 0.1;   % 0 < Gamma < 0.25 smoothing parameter
%----------------------------------------------------------
% % Definition of Particles:
% % FIRST PARTICLE
x1 = 0;
y1 = 0;
% SECOND PARTICLE
x2 = linspace(-L/2,L/2,n+1)' * ones(1,n+1);
x2 = x2(1:end-1, 1:end-1);  
y2 = ones(1,n+1)' * linspace(-L/2,L/2,n+1);
y2 = y2(1:end-1, 1:end-1);
filter = sqrt(x2.^2 + y2.^2) < 3.5 .* r_e;
% evaluation of distance
distance12 = sqrt ( x2.^2 + y2.^2);
%----------------------------------------------------------
% % Definition of Potential:
%----------------------------------------------------------
% Morse Potential 
v2_12 = Morse2D( epsilon, alpha , r_e, distance12 , n );
% figure(1);
% surf(x2,y2,v2_12,'edgealpha',0);
%----------------------------------------------------------
% Lennard Jones Potential:
% epsilon2 = 0.415;
% sigma = 2.27;
% v2_12 = lennard_jones( sigma, epsilon2, distance12 , n );
% figure(3);
% surf(x2,y2,v2_12,'edgealpha',0);
%----------------------------------------------------------
% THIRD PARTICLE
x3 = x2;
y3 = y2;
% evaluation of distance & Potential v_2
distance13 = distance12;
v2_13 = v2_12;
%----------------------------------------------------------
% SECOND PARTICLE     1---2
lap_v2_12 = laplacian3( v2_12, delta);
grad_v2_12_x  = nablax(v2_12,delta);
grad_v2_12_y  = nablay(v2_12,delta); 
% THIRD PARTICLE      1---3
lap_v2_13 = lap_v2_12;
grad_v2_13_x = grad_v2_12_x;
grad_v2_13_y = grad_v2_12_y;
% Second and Third    2---3
lap_v2_23 = circshift(lap_v2_12,[n/2 , n/2]);
grad_v2_23_x = circshift(grad_v2_12_x,[n/2 , n/2]);
grad_v2_23_y = circshift(grad_v2_12_y,[n/2 , n/2]);
%----------------------------------------------------------
diffx = x2(2,1)- x2(1,1);
diffy = y2(1,2)- y2(1,1);
dA = diffx * diffy;
%----------------------------------------------------------
% values relates to P_0 that do not change as we iterate to a solution:
% SECOND PARTICLE     1---2
p0_12 = exp(-beta.* v2_12);
% p0_12 = p0_12 ./ p0_12;
% test = (N * (N-1) / A) / (sum(sum(p0_12)) * dA)
p0_12 = p0_12 * (N * (N-1) / A) / (sum(sum(p0_12)) * dA);
Neg = zeros(total_step,1);
% val = zeros(total_step,1);
var_plot = zeros(n,n,total_step/ploting_step);
% var_plot2 = zeros(n,n,total_step/ploting_step);
% plot_noise = zeros(n,n,20);
% max_delta = 0;
% abs_delta = [];
% ii = 1;
iii =1;
cut = 0;
% p0_12 = p0_12_S;
start = 1;
% %----------------------------------------------------------
for s = start:total_step
    lap_p0_12 = laplacian3(p0_12,delta);
    grad_p0_12_x = nablax(p0_12,delta);
    grad_p0_12_y = nablay(p0_12,delta);
    p0_23 = circshift(p0_12,[n/2, n/2]);
    % Evaluation of formula :
    %----------------------------------------------------------
    %%% 1st term
    first = KbT .* lap_p0_12;
%     figure(5); surf(x2,y2,first,'edgealpha',0);
%     %----------------------------------------------------------
    %%% 2nd t erm 
    second = lap_v2_12 .* p0_12 + ...
    grad_p0_12_x.* grad_v2_12_x+ grad_p0_12_y.* grad_v2_12_y; 
%     figure(6);surf(x2,y2,second,'edgealpha',0);
    %----------------------------------------------------------
    %%% 3rd term 
    % 1st Integral
    fw13 = fft2(second);
    fp1 = fft2( p0_23 );
%     test = fw13 .* fp1;
%     test = ifft2 ( test , 'symmetric');
    int1 = ifft2 ( fw13 .* fp1 );
    int1 = int1 * (dA );
    integral1 = int1 .* p0_12 ;
    trd = integral1./(constantP^3);
%     figure(7);surf(x2,y2,trd,'edgealpha',0); 
    % 3rd Integral
    p0_13 = p0_12;
    fwx3 = fft2 ( p0_13 .* grad_v2_13_x);
    fwy3 = fft2 ( p0_13 .* grad_v2_13_y);
    fp3 = fp1;
    intx3 = ifft2( fwx3 .* fp3 );
    inty3 = ifft2( fwy3 .* fp3 );
    integral3 = (intx3 .* grad_p0_12_x + inty3 .* grad_p0_12_y)* dA;
    fth = integral3./(constantP^3);
%     figure(8);surf(x2,y2,fth,'edgealpha',0);
%----------------------------------------------------------
    change = 2 * G * ( first + second + trd + fth); % dP_12
    p0_12 = p0_12 + change * timestep;
    disp(s)
    % Data for ploting of p0_12 :
    if rem(s,ploting_step)== 0
    var_plot(:,:,iii) = p0_12;
    iii = iii+1;
    end
    %     % Negative part 
    Neg(s) = min(min(p0_12));
%     for xx=1:n
%         for yy=1:n
%             if p0_12(xx,yy) < 0
%             p0_12(xx,yy)= Neg_parameter .* p0_12(xx,yy);
%             end
%          end  
%     end
% disp([num2str(s),'   ', 'Min: ', num2str(min(min(p0_12))),'   ','Max: ',num2str(max(max(p0_12)))]);
% %----------------------------------------------------------
%     figure(1);surf(x2,y2,p0_12,'edgealpha',0);
%     az = 90;
%     el = 0;
%     view([az,el])
%     axis tight
%     title(['Step = ',num2str(s),'  ','dt =',num2str(timestep),'  ','T =',num2str(T)]);    
end
T = 1450;
hh = zeros(n,n,total_step/ploting_step);
for fig_num = 1: total_step/ploting_step
    figure(fig_num);
    hh(fig_num) = surf(y2,x2,var_plot(:,:,fig_num),'edgealpha',0);
    colorbar
    az = 90;
    el = 90;
    view([az,el])
    axis tight
    title(['Step = ',num2str(fig_num * ploting_step),'  ','dt =',num2str(timestep),'  ','T =',num2str(T),'  ','N =',num2str(Neg(fig_num * ploting_step))]);    
    saveas(hh(fig_num),sprintf('figure_Zdirect_%d.jpg',fig_num));
    az = 90;
    el = 0;
    view([az,el])
    axis tight
    title(['Step = ',num2str(fig_num * ploting_step),'  ','dt =',num2str(timestep),'  ','T =',num2str(T),'  ','N =',num2str(Neg(fig_num * ploting_step))]);    
    saveas(hh(fig_num),sprintf('figure_Xdirect_%d.jpg',fig_num));
%     if ( fig_num * ploting_step ) > cut
%         break
%     end
    if rem(fig_num,3) == 0 
        close all
    end
end
save ('Twobody_0.125_100k.mat','p0_12','v2_12','dA','grad_p0_12_x','grad_p0_12_y','grad_v2_12_x','grad_v2_12_y','lap_v2_12','y2','x2')
%----------------------------------------------------------                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
    % 2nd Integral
%     wx23 = grad_v2_23_x .* grad_p0_23_x;
%     wy23 = grad_v2_23_y .* grad_p0_23_y;
%     fw23 = fft2( lap_v2_23 .* p0_23 + wx23 + wy23 );
%     fp2 = fft2(p0_13);
%     int2 = ifft2 ( fw23 .* fp2 );
    %----------------------------------------------------------
%     %----------------------------------------------------------
%     % 4th Integral
% %     fwx4 = fft2 ( p0_23 .* grad_v2_23_x );
% %     fwy4 = fft2 ( p0_23 .* grad_v2_23_y );
% %     fp4 = fp2;
% %     intx4 = ifft2( fwx4 .* fp4 );
% %     inty4 = ifft2( fwy4 .* fp4 );            
%     %----------------------------------------------------------
%     integral1 = int1 .* p0_12 * dA;
% %     integral2 = int2 .* p0_12 * dA;
%     integral3 = (intx3 .* grad_p0_12_x + inty3 .* grad_p0_12_y) * dA;
% %     integral4 = (intx4 .* grad_p0_12_x + inty4 .* grad_p0_12_y) * dA;
%     third = (integral1 + integral3 )./(constantP^3);
%     %----------------------------------------------------------
%     change = 2 * G * ( first + second + third); % dP_12
%     p0_12 = p0_12 + change * timestep;
%     p0_12 = p0_12 * (N * (N-1) / A) / (sum(sum(p0_12)) * dA);
%     tt = (N * (N-1) / A) / (sum(sum(p0_12)) * dA);
%     disp(tt)
%     %----------------------------------------------------------
%     % Negative part 
%     Neg(s) = min(min(p0_12));
%     for xx=1:n
%         for yy=1:n
%             if p0_12(xx,yy) < 0
%             p0_12(xx,yy)= Neg_parameter .* p0_12(xx,yy);
%             end
%          end  
%     end
% %     disp(num2str(s))
% disp([num2str(s),'   ', 'Min: ', num2str(min(min(p0_12))),'   ','Max: ',num2str(max(max(p0_12)))]);
% %----------------------------------------------------------
%     % Data for ploting of p0_12 :
%     if rem(s,ploting_step)== 0
%     var_plot(:,:,iii) = p0_12;
%     iii = iii+1;
%     end
% %----------------------------------------------------------
%     val = sum(sum(abs(change(filter)/total_step)));
%     abs_delta(end+1) = val;
%     if val > max_delta
% % %         max_delta = val;
%     elseif s > 10000 && (val < max_delta/2) && mean(abs_delta((end-9):(end-5)))-mean(abs_delta(((end-4):end)))< 1e-4
%             disp(['Cut Here: Step = ', num2str(s)])
%             cut = s;
%             disp('Stopping criteria worked')
%             break
%     end
%     if(min(min(p0_12)) < -5e-3) && (min(min(p0_12) < 0))
%         disp('Break due to the negative value increase')
%         cut = s;
%         break
%     end
% %     runcheck = 1000;
% %     if rem(s,runcheck)== 0;
% %         figure(1);
% %     surf(y2,x2,p0_12,'edgealpha',0);
% %         axis tight
% %     end
% % end
% % Plot :
% % % PLot all code :  
% % save ('TB_0.1077_400k.mat','p0_12','v2_12','dA','grad_p0_12_x','grad_p0_12_y','grad_v2_12_x','grad_v2_12_y','lap_v2_12','y2','x2')
% % close all
% % fig_num = 1;
% % hh = zeros(n,n,total_step/ploting_step);
% % for fig_num = 1: total_step/ploting_step
% % for fig_num = 1: 10
% %     figure(fig_num);
% %     hh(fig_num) = surf(y2,x2,var_plot(:,:,fig_num),'edgealpha',0);
% %     colorbar
% %     az = 90;
% %     el = 90;
% %     view([az,el])
% %     axis tight
% %     title(['Step = ',num2str(fig_num * ploting_step),'  ','dt =',num2str(timestep),'  ','T =',num2str(T),'  ','N =',num2str(Neg(fig_num * ploting_step))]);    
% %     saveas(hh(fig_num),sprintf('figure_Zdirect_%d.jpg',fig_num));
% %     az = 90;
% %     el = 0;
% %     view([az,el])
% %     axis tight
% %     title(['Step = ',num2str(fig_num * ploting_step),'  ','dt =',num2str(timestep),'  ','T =',num2str(T),'  ','N =',num2str(Neg(fig_num * ploting_step))]);    
% %     saveas(hh(fig_num),sprintf('figure_Xdirect_%d.jpg',fig_num));
% % %     if ( fig_num * ploting_step ) > cut
% % %         break
% % %     end
% %     if rem(fig_num,3) == 0 
% %         close all
% %     end
% % end
% % close all
% % save ('save_plot_450k.mat','var_plot','-v7.3');
% %----------------------------------------------------------
% % for xx = 1:n
% %     for yy = 1:n
% %         if p0_12(xx,yy)<0
% %             p0_12(xx,yy)=0;
% %         end
% %     end
% % end    
% % gamma = 0.1;
% % p0_12_S = smooth ( p0_12 , gamma, 5 );
% %----------------------------------------------------------
% % Plot :
% PLot all code :  
% hh = zeros(n,n,total_step/ploting_step);
% for fig_num = 1: total_step/ploting_step
%     figure(fig_num);
%     hh(fig_num) = surf(y2,x2,var_plot(:,:,fig_num),'edgealpha',0);
%     colorbar
%     az = 90;
%     el = 90;
%     view([az,el])
%     axis tight
%     title(['Total steps = ',num2str(fig_num * ploting_step),'   ','n =',num2str(n),'   ','K_bT =',num2str(KbT),'   ','N =',num2str(N)]);
%     saveas(hh(fig_num),sprintf('figure_Zdirect_%d.jpg',fig_num));
%     az = 90;
%     el = 0;
%     view([az,el])
%     axis tight
%     title(['Total steps = ',num2str(fig_num * ploting_step),'   ','n =',num2str(n),'   ','K_bT =',num2str(KbT),'   ','N =',num2str(N)]);
%     saveas(hh(fig_num),sprintf('figure_Xdirect_%d.jpg',fig_num));
% %     if ( fig_num * ploting_step ) > cut
% %         break
% %     end
%     if rem(fig_num,3) == 0 
%         close all
%     end
% end
% % close all
% % % %----------------------------------------------------------
% % Plot Negative Part
% % fig_num = fig_num + 1;
% % figure(fig_num)
% % ss = linspace(1,total_step,total_step); 
% % h1 = plot(ss,Neg);
% % title('Negative Part ');
% % saveas(h1,'negative.jpg')
% % %----------------------------------------------------------
% % % Plot absolute value Partabsolute change for steady state
% % fig_num = fig_num + 1;
% % figure(fig_num)
% % ps = size(abs_delta);
% % plot_abs = linspace(1,ps(2),ps(2))';
% % h2 = plot(plot_abs,abs_delta);
% % title(['maximum absolute','   ', 'Total steps = ',num2str(fig_num * ploting_step),'   ','K_bT =',num2str(KbT),'   ','N =',num2str(N)]);
% % saveas(h2,'maximum_absolute.jpg')
% % %----------------------------------------------------------
% % % Saving for 1-body distribution :
% % % save smooth p0_12 v2_12 dA grad_p0_12_x grad_p0_12_y grad_v2_12_x grad_v2_12_y lap_v2_12 y2 x2
